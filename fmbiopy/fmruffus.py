""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""

import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmsam as fmsam
import fmbiopy.fmsystem as fmsystem
from fmbiopy.fmtype import StringOrSequence
import logging
import os
from ruffus.proxy_logger import make_shared_logger_and_proxy
from ruffus.proxy_logger import setup_std_shared_logger
import typing

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
    log_results
        If True, stdout and stderr are logged to the pipeline's logfile
    param
        A list of bash parameters e.g ['--foo', 'bar', '-x', '1']

    Attributes
    ----------
    exit_code : int
       The exit code returned by the task's bash process
    """
    _logger = ROOT_LOGGER

    def __init__(self,
                 input_files: typing.Sequence[str],
                 output_files: typing.Sequence[str],
                 log_results: bool = True,
                 param: typing.Sequence[str] = None) -> None:
        # Store parameters
        self._input_files = input_files
        self._output_files = output_files
        self._log_results = log_results
        self._param = param

        # Init attributes
        self._command = None
        self.exit_code = None
        self._shell = False

    def _construct_command(self)-> None:
        """Subclasses construct their own tasks"""
        pass

    def _run_command(self)-> None:
        """Run a command set with `_construct_command`"""
        if self._command:
            if self._log_results:
                _logger.write_header(['Running:'] + self._command)

            self.exit_code = fmsystem.run_command(
                    command=self._command,
                    log_stdout=self._log_results,
                    log_stderr=self._log_results,
                    mutex_log=_logger,
                    shell=self._shell)

        return self.exit_code

class Bowtie2IndexFasta(RuffusTask):
    """Index a .fasta file using bowtie2-build

    Parameters
    ----------
    input_files
        Name of input .fasta file
    output_files
        Name of output bowtie2 index prefix
    log_results
        If True, stdout and stderr are logged to the pipeline's logfile
    param
        A list of bash parameters e.g ['--foo', 'bar', '-x', '1']
    """
    def __init__(self,
                 input_files: typing.Sequence[str],
                 output_files: typing.Sequence[str],
                 log_results: bool = True,
                 param: typing.Sequence[str] = None) -> None:
        super().__init__(input_files, output_files, log_results, param)
        self._construct_command()
        self._run_command()

    def _construct_command(self) -> typing.List[str]:
        """Construct the bash command"""
        self._command = ['bowtie2-build'] + param + [input_files, output_files]


def samtools_index_fasta(input_fasta: str,
        output_index: str,
        logger: RuffusLog = None,
        log_results: bool = False
        ) -> int:
    """Index a list of .fasta files using samtools faidx

    Parameters
    ----------
    input_fasta
        Path to fasta file to be indexed
    output_index
        Path to the output .fai file
    logger
        The Ruffus logging instance
    log_results
        If False, results will not be logged to file
    Returns
    -------
    The error code of the process

    """

    # Construct the samtools command
    command = ['samtools', 'faidx', input_fasta]

    # Run the command
    exit_code = _run_ruffus_command(command, logger, log_results)
    return exit_code


def gunzip(
        input_file: str,
        output_file: str,
        param: typing.List[str] = []) -> int:
    """Gunzip a file and keep the original

    Parameters
    ----------
    input_file
        Path to input file
    output_file
        Path to the gunzipped output
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    """
    command = ['gunzip', '-c'] + param + [input_file]
    exit_code, stdout, _ = fmsystem.run_command(
            command, log_stdout=False, log_stderr=False)
    with open(output_file, "w") as f:
        f.write(stdout)
    return exit_code


def gzip(
        input_file: str,
        output_file: str,
        param: typing.List[str] = []) -> int:
    """Gzip a file

    Parameters
    ----------
    input_file
        Path to input file
    output_file
        Path to the gzipped output
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    """
    command = ['gzip'] + param + [input_file]
    exit_code = fmsystem.run_command(
            command, log_stdout=False, log_stderr=False)[0]
    return exit_code


def paired_bowtie2_align(
        input_files: typing.Tuple[str, str, str],
        output_bam: str,
        param: typing.List[str] = [],
        logger: RuffusLog = None,
        log_results: bool = False) -> int:
    """Align a pair of fastq files to a fasta file using bowtie2

    Parameters
    ----------
    input_files
        A Tuple of the form (Forward reads, Reverse reads, Bowtie2 index)
    output_bam
        Path to the output .bam file
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    logger
        The Ruffus logging instance
    log_results
        If False, results will not be logged to file
    Returns
    -------
    The exit code of the process
    """

    fwd_reads = input_files[0]
    rev_reads = input_files[1]
    bowtie2_index = input_files[2]

    # Construct command
    output_sam = fmpaths.replace_suffix(output_bam, '.bam', '.sam')
    command = ['bowtie2', '-1', fwd_reads, '-2', rev_reads, '-x',
            bowtie2_index, '-S', output_sam]

    # Run command
    exit_code = _run_ruffus_command(command, logger, log_results)

    # Convert to a sorted Bam
    fmsam.sam_to_bam(output_sam, output_bam)
