""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""

import logging
from pathlib import Path
from typing import Any
from typing import Callable
from typing import List
from typing import Sequence
from typing import Type

from ruffus.proxy_logger import make_shared_logger_and_proxy
from ruffus.proxy_logger import setup_std_shared_logger

import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmsystem as fmsystem

# Default logging format"""
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
        The ID of the logger
    location
        The path to where the logfile should be stored
    formatter
        Specifies logging format
    delay
        If True, file opening deferred until first log action
    level
        Logging level, see logging module docs
    buffered
        If True, the `write` statements are appended to the buffer, but nothing
        is logged until the buffer is `flush`ed.

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
            location: Path,
            formatter: str = DEFAULT_LOG_FORMAT,
            delay: bool = True,
            level: int = logging.DEBUG,
            buffered: bool = True) -> None:

        self.config = {'file_name': str(location),
                       'formatter': formatter,
                       'delay': delay,
                       'level': level,
                       'buffered': buffered}

        if self.config['buffered']:
            self._buffer = ''

        if not location.parent.exists():
            raise ValueError('Parent directory of logfile doesnt exist')

        self.log, self.mutex = make_shared_logger_and_proxy(
                setup_std_shared_logger, name, self.config)

    def write(self, msg: str) -> None:
        """Log a message with mutex lock

        Parameters
        ----------
        msg
            Message to be logged to the Ruffus logfile. List is converted to
            string and spaces are added. If already a string just log directly
        """
        with self.mutex:
            if self.config['buffered']:
                self._buffer += '\n' + msg
            else:
                self.log.info(msg)

    def _divider(self, sub: bool):
        """Write a text divider to logfile

        Parameters
        ----------
        sub
            If True, the divider is made up of dashes. Else, the divider is
            made up of equal signs
        """
        if sub:
            self.write('-' * 80)
        else:
            self.write('=' * 80)

    def write_header(self, msg: str, sub=False) -> None:
        """Write a header to the logfile

        Parameters
        ----------
        msg
            A string or list of words to be written as the header text
        sub, optional
            If True, write a subheader, instead of a header
        """
        self._divider(sub=sub)
        self.write(msg=msg)
        self._divider(sub=sub)

    def flush(self):
        """Write the buffer to the logfile and clear it"""
        self.log.info(self._buffer)
        self._buffer = ''


# The `RuffusLog` shared across all `RuffusTasks`"""
ROOT_LOGGER = RuffusLog("", Path("pipeline.log"))


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
    input_type: List[str] = ['']
    output_type: List[str] = ['']
    _logger = ROOT_LOGGER

    def __init__(self,
                 input_files: Sequence[str],
                 output_files: List[str],
                 log_results: bool = True,
                 param: Sequence[str] = [''],
                 run_on_init: bool = True) -> None:
        # Store parameters
        if isinstance(input_files, str):
            self._input_files: Sequence[str] = [input_files]
        else:
            self._input_files = input_files
        if isinstance(output_files, str):
            self._output_files: List[str] = [output_files]
        else:
            self._output_files = output_files
        self._log_results = log_results
        self._param = param

        # Init attributes

        self.exit_code: List[int] = []
        self.stdout: List[str] = []
        self.stderr: List[str] = []

        # Bash command as list
        self._command: List = []

        # If True, command run directly in shell
        if not hasattr(self, '_shell'):
            self._shell: bool = False

        # If True, the command is a list of subcommands, which should be run
        # successively. This is essentially a mini-pipeline.
        if not hasattr(self, '_multi_stage'):
            self._multi_stage: bool = False

        # If True, the input files are deleted once the output files are
        # produced.
        if not hasattr(self, '_inplace'):
            self._inplace: bool = False

        # Extra output files produced which are not specified by the user.
        self._extra_outputs: List[str] = []
        self._add_extra_outputs()
        self._construct_command()

        if run_on_init:
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
            if self._log_results:
                header = 'Task: ' + self.__class__.__name__ + '\n'
                header += 'Input_files: ' + ', '.join(self._input_files) + '\n'
                header += 'Output_files: ' + ', '.join(self._output_files)
                self._logger.write_header(header)
            if not self._multi_stage:
                # Convert to list of lists so we can loop through anyway
                self._command = [self._command]

            for subcommand in self._command:

                if self._log_results:
                    self._logger.write_header(
                            ' '.join(['Running:'] + subcommand),
                            sub=True)

                (exit_code, stdout, stderr) = (
                    fmsystem.run_command(
                        command=subcommand,
                        log_stdout=True,
                        log_stderr=True,
                        mutex_log=self._logger,
                        shell=self._shell))
                self.exit_code.append(exit_code)
                self.stdout.append(stdout)
                self.stderr.append(stderr)

            # Write the task's log
            self._logger.flush()

            if self._inplace:
                paths = fmpaths.as_path(self._input_files)
                fmsystem.remove_all(paths, silent=True)

        except Exception:
            self._cleanup()
            raise

    def _cleanup(self)-> None:
        fmsystem.remove_all(
                fmpaths.as_path(self._output_files + self._extra_outputs),
                silent=True)


class SamtoolsIndexFasta(RuffusTask):
    """Index a .fasta file using samtools faidx"""

    input_type = ['fasta']
    output_type = ['fasta.fai']

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
    input_type = ['ANY']
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
    input_type = ['fasta', 'fwd_fastq', 'rev_fastq']
    output_type = ['bam']

    def __init__(self, *args, **kwargs)-> None:
        self._multi_stage = True
        self._shell = True
        super().__init__(*args, **kwargs)

    def _add_extra_outputs(self) -> None:
        """Add the .bt2 and .bai files to the output file list"""
        prefix = Path(self._input_files[0]).stem
        bowtie2_indices = fmpaths.get_bowtie2_indices(prefix)
        bai_file = '.'.join([self._output_files[0], 'bai'])

        self._extra_outputs = fmlist.flatten([
                self._output_files, bowtie2_indices, bai_file])

    def _construct_command(self) -> None:
        """Construct the bash command"""
        fasta = self._input_files[0]
        fwd_fastq = self._input_files[1]
        rev_fastq = self._input_files[2]
        output_bam = self._output_files[0]

        # Index fasta first
        bowtie2_index = Path(fasta).stem
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


class SymlinkInputs(RuffusTask):
    """Create symlinks of the input arguments"""
    input_type = ['ANY']
    output_type = ['SAME']

    def _construct_command(self) -> None:
        """Construct the bash command"""
        self._command = fmlist.flatten([
            'ln', '-s', self._input_files, self._output_files])


"""Type variable for a function with same inputs and outputs as a RuffusTask"""
TaskFunction = Callable[[Sequence[str], Sequence[str], Any], None]


def apply(task: Type[RuffusTask])-> TaskFunction:
    """Return function which applies a `RuffusTask` to multiple inputs

    All extra parameters are passed to `task`

    Parameters
    ----------
    task
        A RuffusTask class

    Returns
    -------
    A function of two arguments, 1. a sequence of input files, 2. a
    sequence of output files, which maps `task` from input[i] -> output[i]
    """
    def _apply_task(
            input_files: Sequence[str],
            output_files: Sequence[str],
            *args,
            **kwargs,
            )-> None:
        for inp, out in zip(input_files, output_files):
            task([inp], [out], *args, **kwargs)

    return _apply_task
