""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""
from abc import (
        ABC,
        abstractmethod,
        )
from logging import DEBUG
from pathlib import Path
from typing import (
        Any,
        Callable,
        List,
        Sequence,
        Tuple,
        Type,
        )

from ruffus.proxy_logger import (
        make_shared_logger_and_proxy,
        setup_std_shared_logger,
        )

from fmbiopy.fmlist import (
        ensure_list,
        exclude_blank,
        flatten,
        )
from fmbiopy.fmpaths import (
        add_suffix,
        as_paths,
        get_bowtie2_indices,
        )
from fmbiopy.fmsystem import (
        remove_all,
        run_command,
        )


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
    overwrite
        If True, logfile will be overwritten if it already exists

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
            level: int = DEBUG,
            buffered: bool = True,
            overwrite: bool = True) -> None:

        self.config = {'file_name': str(location),
                       'formatter': formatter,
                       'delay': delay,
                       'level': level,
                       'buffered': buffered}

        if self.config['buffered']:
            self._buffer: str = ''

        if not location.parent.exists():
            raise FileNotFoundError('Parent directory of logfile doesnt exist')

        # Overwrite logfile if it already exists
        if overwrite:
            if location.exists():
                location.unlink()

        self.log, self.mutex = make_shared_logger_and_proxy(
                setup_std_shared_logger, name, self.config)

    def flush(self)-> None:
        """Write the buffer to the logfile and clear it"""
        self.log.info(self._buffer)
        self._clear_buffer()

    def write(self, msg: str)-> None:
        """Log a message with mutex lock

        Parameters
        ----------
        msg
            Message to be logged to the Ruffus logfile. List is converted to
            string and spaces are added. If already a string just log directly
        """
        with self.mutex:
            if self.config['buffered']:
                self._new_msg(msg)
            else:
                self.log.info(msg)

    def write_header(self, msg: str, subheader: bool = False)-> None:
        """Write a header to the logfile

        Parameters
        ----------
        msg
            A string or list of words to be written as the header text
        sub, optional
            If True, write a subheader, instead of a header
        """
        self._write_divider(subheader)
        self.write(msg)
        self._write_divider(subheader)

    def _clear_buffer(self)-> None:
        self._buffer = ''

    def _new_msg(self, msg: str)-> None:
        """Add a new message to the buffer"""
        self._buffer += ''.join(['\n', msg])

    def _write_divider(self, subheader: bool)-> None:
        """Write a text divider to logfile

        Parameters
        ----------
        subheader
            If True, the divider is made up of dashes. Else, the divider is
            made up of equal signs
        """
        if subheader:
            self.write('-' * 80)
        else:
            self.write('=' * 80)


"""The `RuffusLog` shared across all `RuffusTransforms`"""
ROOT_LOGGER = RuffusLog("", Path("pipeline.log"))

class RuffusTask(ABC):
    """Abtract class representing a `ruffus` function with output but no input

    The running and task logging of the command is handled here. Extra
    initialization and command construction are handled by subclasses.
    `RuffusFunctions` which inherit from `RuffusTask` will produce data which
    will be input to `RuffusTransform`s. In practice, these functions will be used
    for `RuffusFunction`s called during `ruffus.Pipeline.originate` tasks.

    Subclasses are expected to define their own method for constructing their
    command (`_build`). They may also define extra outputs
    (`_add_extra_outputs`) and steps to be run after the main command
    (`_post_command`).

    Parameters
    ----------
    log_results: optional
        If True, OS stdout and stderr are logged to the pipeline's logfile
    param: optional
        A list of bash parameters e.g ['--foo', 'bar', '-x', '1']
    shell: optional
        If True, command should be run in shell rather than through python
    inplace: optional
        If True, input files are deleted once output files are produced

    Attributes
    ----------
    output_type : str
        (Class attribute) Expected output file extension. If '' then the
        extension of the input is stripped in the output.
    exit_code: int
        The exit code returned by the task's bash process
    stdout: str
        The standard output produced by the task
    stderr: str
        The standard error produced by the task
    """

    _logger = ROOT_LOGGER
    output_type: List[str] = None

    def __init__(
            self,
            param: List[str] = [],
            log_results: bool = True)-> None:

        # Store parameters
        self._log_results = log_results
        self.param = param

        # If True, command will be run in shell
        self._shell = False

        # Initialize attributes
        # ---------------------

        # Bash command as list
        self._command: List = []

        # Returns from command(s)
        self.exit_code: List[int] = []
        self.stdout: List[str] = []
        self.stderr: List[str] = []

        # Output files as strings and paths
        self.output_files: Sequence[str] = None
        self._output_paths: List[Path] = None

        # Output directories
        self._outdirs: List[Path] = None

        # Output file stems
        self._output_stems: Sequence[str] = None

        # Extra outputs produced which are not specified by the user
        self._extra_outputs: Sequence[Path] = None

        # The header to be written to the logfile at the start of the task
        self._task_header: str = ''

    def cleanup(self)-> None:
        """Cleanup"""
        if self._output_paths is not None:
            remove_all(self._output_paths, silent=True)
        if self._extra_outputs is not None:
            remove_all(self._extra_outputs, silent=True)

    def run(self, *args):
        """Run a command

        Positional arguments are passed from the `ruffus.Pipeline`. If the task
        has inputs this will be the first argument. Otherwise the first and
        only argument will be the outputs."""
        output_files = args[0]
        self._define_output_attributes(output_files)
        self._check_run_inputs()
        self._build()
        self._execute()
        self._post_command()

    def _add_command(self, subcommand: Sequence)-> None:
        """Append a command to the command list"""
        self._command.append(exclude_blank(flatten(subcommand)))

    def _add_extra_outputs(self)-> List[Path]:
        """Add additional non-specified outputs to the output file list

        Should be overwritten by subclasses which produce outputs other than
        those defined in output_files parameter

        Returns
        -------
        A list of the extra outputs
        """
        return None

    def _add_taskname_to_header(self)-> None:
        self._add_to_header('Task', self.__class__.__name__)

    def _addoutput_files_to_header(self)-> None:
        self._add_to_header('Output', self.output_files)

    def _add_to_header(self, key: str, value: Sequence[str])-> None:
        if not isinstance(value, str):
            value = ', '.join(value)
        self._task_header += ' '.join([key, ':', value, '\n'])

    @abstractmethod
    def _build(self)-> None:
        """Construct the bash command"""
        pass

    def _check_run_inputs(self)-> None:
        """Check that the inputs to `run` are valid"""
        pass

    def _define_output_attributes(self, output_files: Sequence[str])-> None:
        """Define attributes related to output files"""
        self.output_files = ensure_list(output_files)
        self._output_paths: List[Path] = as_paths(self.output_files)
        self._outdirs: List[Path] = [f.parent for f in self._output_paths]
        self._output_stems = self._get_stems(self._output_paths)

        # Extra output files which are not specified by the user.
        self._extra_outputs: List[Path] = self._add_extra_outputs()

        # Create output directories if necessary
        self._mkdirs()

    def _execute(self)-> None:
        """Run a command set with `_build`"""
        try:
            if self._log_results:
                self._log_task_header()

            for subcommand in self._command:

                if self._log_results:
                    self._logger.write_header(
                            ' '.join(['Running:'] + subcommand),
                            subheader=True)

                results = (
                        run_command(
                            subcommand,
                            log=(True, True),
                            mutex_logger=self._logger,
                            shell=self._shell))

                self._store_results(results)

            # Write the task's log
            self._logger.flush()

        except Exception:
            self.cleanup()
            raise

    def _get_stems(self, paths: List[Path])-> List[str]:
        """Get the filename stems (everything but the extension)"""
        stems = []
        for path in paths:
            if path.suffix == '.gz':
                stems.append(Path(path.stem).stem)
            else:
                stems.append(path.stem)
        return stems

    def _log_task_header(self)-> None:
        """Write the task header to the logfile"""
        self._add_taskname_to_header()
        self._addoutput_files_to_header()
        self._logger.write_header(self._task_header)

    def _mkdirs(self)-> None:
        """Create output directories if necessary"""
        for d in self._outdirs:
            d.mkdir(parents=True, exist_ok=True)

    def _post_command(self)-> None:
        """Commands to be run after `run`"""
        pass

    def _store_results(self, results: Tuple[int, str, str])-> None:
        """Update attributes with the results of `fmsystem.run_command`"""
        self.exit_code.append(results[0])
        self.stdout.append(results[1])
        self.stderr.append(results[2])


class RuffusTransform(RuffusTask, ABC):
    """A superclass representing a `ruffus` function with output and input

    `RuffusFunctions` which inherit from `RuffusTransform` will produce data which
    will be input to further `RuffusTransform`s. The vast majority of functions used
    in `ruffus` pipelines fall under the domain of `RuffusTransform`s.

    Subclasses are expected to define their own method for constructing their
    command (`_build`). They may also define extra outputs
    (`_add_extra_outputs`) and steps to be run after the main command
    (`_post_command`).

    Attributes
    ----------
    input_type: str
        (Class attribute) Expected input file extension. If None then any file
        type is acceptable input.
    """
    input_type: List[str] = None

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        # Initialize attributes
        # ---------------------

        # Input files as strings and paths
        self.input_files: Sequence[str] = None
        self._input_paths: List[Path] = None

        # Input file stems
        self._input_stems: List[str] = None

        # Input directories
        self._indirs: List[Path] = None

        # If True, input files will be deleted upon creation of output files
        self._inplace = False

    def run(self, *args):
        """Run a command

        Positional arguments are passed from the `ruffus.Pipeline`. If the task
        has inputs this will be the first argument. Otherwise the first and
        only argument will be the outputs."""
        input_files = args[0]
        output_files = args[1]
        self._define_input_attributes(input_files)
        super().run(output_files)
        if self._inplace:
            remove_all(self._input_paths, silent=True)

    def _addinput_files_to_header(self)-> None:
        self._add_to_header('Output', self.input_files)

    @abstractmethod
    def _build(self)-> None:
        pass

    def _check_run_inputs(self)-> None:
        """Check that input and output files are valid"""
        for inp in self._input_paths:
            if inp in self._output_paths:
                raise FileExistsError('Attempted to overwrite input file(s)')

    def _define_input_attributes(self, input_files)-> None:
        self.input_files: Sequence[str] = ensure_list(input_files)
        self._input_paths: List[Path] = as_paths(self.input_files)
        self._input_stems = self._get_stems(self._input_paths)
        self._indirs: List[Path] = [f.parent for f in self._input_paths]

    def _log_task_header(self)-> None:
        """Write the task header to the logfile"""
        self._add_to_header('Task', self.__class__.__name__)
        self._addinput_files_to_header()
        self._addoutput_files_to_header()
        self._logger.write_header(self._task_header)


class SamtoolsIndexFasta(RuffusTransform):
    """Index a .fasta file using samtools faidx"""

    input_type = ['fasta']
    output_type = ['fasta.fai']

    def _build(self) -> None:
        """Construct the bash command"""
        self._add_command(['samtools', 'faidx', self.input_files, self.param])


class Gunzip(RuffusTransform):
    """Unzip a file with gunzip"""
    input_type = ['gz']
    output_type = ['']

    def __init__(self, *args, **kwargs)-> None:
        super().__init__(*args, **kwargs)
        self._inplace = True
        self._shell = True

    def cleanup(self)-> None:
        """Gzip the file back up"""
        Gzip().run(self.output_files[0], self.input_files[0])

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command(['gunzip', '-c', self.param, self.input_files, '>',
            self.output_files])


class Gzip(RuffusTransform):
    """Compress a file with g zip"""
    input_type = ['ANY']
    output_type = ['gz']

    def __init__(self, *args, **kwargs)-> None:
        super().__init__(*args, **kwargs)
        self._inplace = True
        self._shell = True

    def cleanup(self) -> None:
        """Unzip the files again"""
        Gunzip().run(self.output_files[0], self.input_files[0])

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command([
            'gzip', '-c', self.param, self.input_files, '>',
            self.output_files])


class PairedBowtie2Align(RuffusTransform):
    """Align a pair of fastq files to a fasta file using bowtie2"""
    input_type = ['fasta', 'fwd_fastq', 'rev_fastq']
    output_type = ['bam']

    def __init__(self, *args, **kwargs)-> None:
        super().__init__(*args, **kwargs)
        self._shell = True

    def _add_extra_outputs(self) -> List[Path]:
        """Add the .bt2 and .bai files to the output file list"""
        bowtie2_indices = get_bowtie2_indices(self._input_stems[0])
        bai_file: Path = add_suffix(self._output_paths[0], '.bai')
        return flatten([bowtie2_indices, bai_file])

    def _build(self) -> None:
        """Construct the bash command"""
        fasta = self.input_files[0]
        fwd_fastq = self.input_files[1]
        rev_fastq = self.input_files[2]
        output_bam = self.output_files[0]

        # Index fasta first
        bowtie2_index = self._input_stems[0]
        self._add_command([
            'bowtie2-build', fasta, bowtie2_index, '>', '/dev/null', '2>&1'])

        # Move the indices to the assembly directory
        self._add_command(['mv', '*.bt2', str(self._indirs[0])])

        #  Run Bowtie2 and pipe the output to a sorted bam file
        self._add_command([
            'bowtie2', *self.param, '-1', fwd_fastq, '-2', rev_fastq, '-x',
            bowtie2_index, '|', 'samtools', 'view', '-bS', '-', '|',
            'samtools', 'sort', '-f', '-', output_bam])

        # Index the bam file
        self._add_command([
            'samtools', 'index', output_bam, '/dev/null', '2>&1'])


class SymlinkInputs(RuffusTransform):
    """Create symlinks of the input arguments"""
    input_type = ['ANY']
    output_type = ['SAME']

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command([
            'ln', '-sf', self.param, self.input_files, self.output_files])


class BuildCentrifugeDB(RuffusTask):
    """Build a custom centrifuge database"""
    def __init__(self, *args, **kwargs)-> None:
        # Location of the seqid2taxid mapping file
        super().__init__(*args, **kwargs)

        # Location of the seqid2taxid file
        self._seqid2taxid: Path = None

        # Location of the concatenated sequences file
        self._concat_seqs: Path = None

        # Location of the taxonomy directory
        self._taxonomy_dir: Path = None

        # Location of the library directory
        self._library_dir: Path = None
        self._shell = True

    def _add_extra_outputs(self)-> List[Path]:
        """Add the library, taxonomy and SeqID2Taxmap files to extra outputs"""
        self._seqid2taxid = self._outdirs[0].joinpath('seqid2taxid.map')

        # The index files
        centrifuge_idx = [
                ''.join([self.output_files[0], '.', i, '.cf'])
                for i in ['1', '2', '3']]

        return flatten([self._seqid2taxid, centrifuge_idx])

    def _build(self)-> None:
        """Construct the bash command"""
        self._taxonomy_dir = self._outdirs[0].joinpath('taxonomy')
        self._library_dir = self._outdirs[0].joinpath('library')

        # Download the complete NCBI taxonomy
        self._add_command([
            "centrifuge-download", "-o", str(self._taxonomy_dir), "taxonomy"])

        # Download all complete Archaeal, Bacterial, Viral, Fungal, Protozoan
        # and Plant genomes
        self._add_command([
            "centrifuge-download", "-o", str(self._library_dir), "-m", "-d",
            "archaea,bacteria,viral,fungi,protozoa,plant", "refseq", ">",
            str(self._seqid2taxid)])

        # Download human reference sequence
        self._add_command([
            "centrifuge-download", "-o", str(self._library_dir), "-d",
            "vertebrate_mammalian", "-a", "Chromosome", "-t", 9606, "-c",
            "reference_genome", ">>", str(self._seqid2taxid)])

        # Concatenate the downloaded sequences
        self._concat_seqs = self._outdirs[0].joinpath("catseq.fna")
        self._add_command([
            "cat", str(self._library_dir / '*' / '*.fna'), ">",
            str(self._concat_seqs)])

        # Build the database
        self._add_command([
            "centrifuge-build", "--conversion-table", str(self._seqid2taxid),
            "--taxonomy-tree", str(self._taxonomy_dir / 'nodes.dmp'), "--name-table",
            str(self._taxonomy_dir / 'names.dmp'), str(self._concat_seqs),
            self.output_files[0]])

    def _post_command(self)-> None:
        """Cleanup unnecessary large files"""

        # Be as explicit as possible to avoid accidental deletions
        self._concat_seqs.unlink()
        remove_all(self._taxonomy_dir.glob('*.dmp'))
        remove_all(self._library_dir.glob('*/*.fna'))

        combined_contents = list(self._taxonomy_dir.glob('*')) + \
                list(self._library_dir.glob('*'))

        for f in combined_contents:
            try:
                f.rmdir()
            except OSError:
                pass


"""Type variable for a function with same inputs and outputs as a RuffusTransform"""
TaskFunction = Callable[[Sequence[str], Sequence[str], Any], None]


def apply_(task: Type[RuffusTransform])-> TaskFunction:
    """Return function which applies a `RuffusTransform` to multiple inputs

    All extra parameters are passed to `task`

    Parameters
    ----------
    task
        A RuffusTransform class

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
            task(*args, **kwargs).run([inp], [out])

    return _apply_task

class ParameterError(ValueError):
    pass
