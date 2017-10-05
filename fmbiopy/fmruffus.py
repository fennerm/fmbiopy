""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""

import logging
import os
from ruffus.proxy_logger import make_shared_logger_and_proxy
from ruffus.proxy_logger import setup_std_shared_logger
import typing

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmsystem as fmsystem
from fmbiopy.fmtype import StringOrSequence

"""Default logging format"""
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"


class RuffusLog(object):
    """A logger/mutex pair used for logging in Ruffus pipelines

    Usage
    -----
        ruflog = RuffusLog("foo", "bar.log")
        with ruflog.mutex:
            ruflog.log.info("Logged message")

    Parameters
    ----------
    name
        The ID used to refer to the logging instance
    location
        The path to where the logfile should be stored
    formatter
        Specifies logging format
    delay
        If True, file opening deferred until first log action
    level
        Logging level, see logging module docs

    Attributes
    ----------
    config: Dict[str, str, bool, int]
        Arguments passed to logger setup
    log: proxy_logger.LoggerProxy
        The logging proxy instance
    mutex: multiprocessing.managers.AcquirerProxy
        The mutex lock used for synchronous access to logfile
    """
    def __init__(
            self,
            name: str,
            location: str,
            formatter: str = DEFAULT_LOG_FORMAT,
            delay: bool = True,
            level: int = logging.DEBUG) -> None:

        logdir = os.path.dirname(location)
        if logdir:
            if not os.path.exists(logdir):
                raise ValueError("Logfile directory does not exist")

        self.config = {'file_name' : location,
                       'formatter' : formatter,
                       'delay' : delay,
                       'level' : level}

        self.log, self.mutex = make_shared_logger_and_proxy(
                setup_std_shared_logger, name, self.config)

    def write(self, message: StringOrSequence) -> None:
        """Log a message with mutex lock

        Parameters
        ----------
        message
            Message to be logged to the Ruffus logfile. List is converted to
            string and spaces are added. If already a string just log directly
        header, optional
            Header string. Will precede message and be wrapped in dashes to
            stand out from other text.

        """
        with self.mutex:
            if not isinstance(message, str):
                message = ' '.join(message)
            self.log.info(message)

    def _divider(self):
        """Write a text divider to logfile"""
        self.write('-' * 30)

    def write_header(self, message: StringOrSequence) -> None:
        """Write a header to the logfile"""
        self._divider()
        self.write(message)
        self._divider()


"""The `RuffusLog` shared across all `RuffusTasks`"""
ROOT_LOGGER = RuffusLog("", "pipeline.log")


class RuffusTask(object):
    """A superclass representing a task to be run with Ruffus

    Each subclass is expected to define its own method for constructing the
    Bash command to be run. The superclass handles the setup, logging and
    running of the command.

    Parameters
    ----------
    input_files
        One or more input file names
    output_files
        One or more output file names
    log_results : Optional
        If True, stdout and stderr are logged to the pipeline's logfile
    param : Optional
        A list of bash parameters e.g ['--foo', 'bar', '-x', '1']

    Other Parameters
    ----------------
    run_on_init : Optional
        If True, the associated bash command is run upon initialization. If
        False, the command is constructed but not run.

    Attributes
    ----------
    exit_code : int
        The exit code returned by the task's bash process
    stdout : str
        The standard output produced by the task
    stderr : str
        The standard error produced by the task
    input_typ str
        (Class attribute) Expected input file extension. If None then any file
        type is acceptable input.
    output_type : str
        (Class attribute) Expected output file extension. If '' then the
        extension of the input is stripped in the output.
    """
    input_type : typing.List[typing.Iterable] = ['']
    output_type = ['']
    _logger = ROOT_LOGGER

    def __init__(self,
                 input_files: typing.Sequence[str],
                 output_files: typing.Sequence[str],
                 log_results: bool = True,
                 param: typing.Sequence[str] = [],
                 run_on_init: bool = True) -> None:
        # Store parameters
        self._input_files = input_files
        self._output_files = output_files
        self._log_results = log_results
        self._param = param

        # Init attributes

        self.exit_code : int = None
        self.stdout = ''
        self.stderr = ''

        # Bash command as list
        self._command : typing.List = []

        # If True, command run directly in shell
        try:
            self._shell
        except AttributeError:
            self._shell = False

        # If True, the command is a list of subcommands, which should be run
        # successively. This is essentially a mini-pipeline.
        try:
            self._multi_stage
        except AttributeError:
            self._multi_stage = False

        # If True, the input files are deleted once the output files are
        # produced.
        try:
            self._inplace
        except AttributeError:
            self._inplace = False

        # Extra output files produced which are not specified by the user.
        self._extra_outputs : typing.List[str] = []
        self._add_extra_outputs()
        self._construct_command()
        self._run_command()

    def _construct_command(self)-> None:
        """Subclasses construct their own tasks"""
        pass

    def _add_extra_outputs(self) -> None:
        """Add additional non-specified outputs to the output file list"""
        pass

    def _run_command(self)-> None:
        """Run a command set with `_construct_command`"""
        try:
            if not self._multi_stage:
                # Convert to list of lists so we can loop through anyway
                self._command = [self._command]

            for subcommand in self._command:

                if self._log_results:
                    self._logger.write_header(['Running:'] + subcommand)

                (self.exit_code, self.stdout, self.stderr) = (
                    fmsystem.run_command(
                        command=subcommand,
                        log_stdout=self._log_results,
                        log_stderr=self._log_results,
                        mutex_log=self._logger,
                        shell=self._shell))

            if self._inplace:
                if fmcheck.all_filesize_nonzero(self._output_files):
                    fmsystem.remove_all(self._input_files, silent=True)

        except Exception:
            self._cleanup()

    def _cleanup(self)-> None:
        fmsystem.remove_all(
                self._output_files + self._extra_outputs,
                silent=True)


class SamtoolsIndexFasta(RuffusTask):
    """Index a .fasta file using samtools faidx"""

    input_type = ['fasta']
    output_type = ['fai']

    def __init__(self, *args, **kwargs)-> None:
        super().__init__(*args, **kwargs)

    def _construct_command(self) -> None:
        """Construct the bash command"""
        self._command = fmlist.flatten([
            'samtools', 'faidx', self._input_files])


class Gunzip(RuffusTask):
    """Unzip a file with gunzip"""
    input_type = ['gz']
    output_type = ['']

    def __init__(self, *args, **kwargs)-> None:
        self._shell = True
        self._inplace = True
        super().__init__(*args, **kwargs)

    def _construct_command(self) -> None:
        """Construct the bash command"""
        self._command = fmlist.flatten([
            'gunzip', '-c', self._param, self._input_files, '>',
            self._output_files])

    def _cleanup(self)-> None:
        """Gzip the file back up"""
        Gzip(self._output_files[0], self._input_files[0])


class Gzip(RuffusTask):
    """Compress a file with g zip"""
    input_type = ['any']
    output_type = ['gz']

    def __init__(self, *args, **kwargs)-> None:
        self._shell = True
        self._inplace = True
        super().__init__(*args, **kwargs)

    def _construct_command(self) -> None:
        """Construct the bash command"""
        self._command = fmlist.flatten([
            'gzip', '-c', self._param, self._input_files, '>',
            self._output_files])

    def _cleanup(self) -> None:
        """Unzip the files again"""
        Gunzip(self._output_files[0], self._input_files[0])


class PairedBowtie2Align(RuffusTask):
    """Align a pair of fastq files to a fasta file using bowtie2"""
    input_type = [('fastq', 'fastq'), 'fasta']
    output_type = ['bam']

    def __init__(self, *args, **kwargs)-> None:
        self._multi_stage = True
        self._shell = True
        super().__init__(*args, **kwargs)

    def _add_extra_outputs(self) -> None:
        """Add the .bt2 and .bai files to the output file list"""
        prefix = fmpaths.remove_suffix(self._input_files[2])
        bowtie2_indices = fmpaths.get_bowtie2_indices(prefix[0])
        bai_file = fmpaths.add_suffix(self._output_files[0], '.bai')

        self._extra_outputs = fmlist.flatten([
                self._output_files, bowtie2_indices, bai_file])

    def _construct_command(self) -> None:
        """Construct the bash command"""
        fwd_fastq = self._input_files[0]
        rev_fastq = self._input_files[1]
        fasta = self._input_files[2]
        output_bam = self._output_files[0]

        # Index fasta first
        bowtie2_index = fmpaths.remove_suffix(fasta)
        self._command.append(fmlist.flatten([
            'bowtie2-build', fasta, bowtie2_index, '>', '/dev/null', '2>&1']))

        #  Run Bowtie2 and pipe the output to a sorted bam file
        self._command.append(fmlist.flatten([
            'bowtie2', '-1', fwd_fastq, '-2', rev_fastq, '-x', bowtie2_index,
            '|', 'samtools', 'view', '-bS', '-', '|', 'samtools', 'sort',
            '-f', '-', output_bam]))

        # Index the bam file
        self._command.append(fmlist.flatten([
            'samtools', 'index', output_bam, '/dev/null', '2>&1']))
