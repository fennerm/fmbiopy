""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""
from abc import (
        ABC,
        abstractmethod,
        )
from logging import DEBUG
from multiprocessing.managers import AcquirerProxy  #type: ignore
from pathlib import Path
from typing import (
        Any,
        Callable,
        List,
        Sequence,
        Tuple,
        Type,
        )

from ruffus.ruffus_utility import formatter
from ruffus.proxy_logger import (
        LoggerProxy,
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


"""Default logging format"""
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"

"""Binaries"""
SAMTOOLS = "samtools"
BOWTIE2 = "bowtie2"
CENTRIFUGE_DOWNLOAD = "centrifuge-download"
CENTRIFUGE_BUILD = "centrifuge-build"
CENTRIFUGE = "centrifuge"
BOWTIE2_BUILD = "bowtie2-build"

"""Default output directory"""
OUTPUT_DIR: Path = Path.cwd() / 'pipe'

class RuffusLogger(object):
    """A logger/mutex pair used for logging in Ruffus pipelines

    Examples
    --------
    >>>>ruflog = RuffusLogger("foo", "bar.log")
    >>>>with ruflog.mutex:
            ruflog.log.info("RuffusLoggerged message")

    Attributes
    ----------
    config: Dict[str, str, bool, int]
        Arguments passed to logger setup
    log: LoggerProxy
        The logging proxy instance
    mutex: AcquirerProxy
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
        """Initialization

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
            If True, the `write` statements are appended to the buffer, but
            nothing is logged until the buffer is `flush`ed.
        overwrite
            If True, logfile will be overwritten if it already exists

        Raises
        ------
        FileNotFoundError
            Parent directory of logfile location doesn't exist
        """


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

        logprox: Tuple[LoggerProxy, AcquirerProxy] = \
                make_shared_logger_and_proxy(
                        setup_std_shared_logger, name, self.config)
        self.log, self.mutex = logprox

    def flush(self)-> None:
        """Write the buffer to the logfile and clear it"""
        self.log.info(self._buffer)
        self._clear_buffer()

    def write(self, msg: str)-> None:
        """RuffusLogger a message with mutex lock

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


"""The `RuffusLogger` shared across all `Transform`s"""
ROOT_LOGGER = RuffusLogger('', Path("pipeline.log"))


# """Valid input type for Pipeline.schedule"""
# RuffusTaskInput = Union[RuffusRuffusTask, Sequence[str], Sequence[Sequence[str]]]
#
#
# class Pipeline(object):
#     """A bioinformatics pipeline
#
#     Attributes
#     ----------
#     logger : RuffusLogger
#         Object which handles writing to the pipeline logfile.
#     name : str
#         Name of the pipeline
#     """
#     def __init__(self, name: str)-> None:
#         logfile_location: Path = OUTPUT_DIR / ''.join([name, '.log'])
#         self.logger: RuffusLogger = RuffusLogger('', logfile_location)
#         self.name: str = name
#         self._pipeline: RuffusPipeline = RuffusPipeline(name=name)
#         self._idx: int = 1
#
#     def schedule(
#             self,
#             task: Type[RuffusTask],
#             input_files: RuffusTaskInput,
#             input_suffixes: List[str],
#             output_suffixes: Iterable[str],
#             param: List[str],
#             )-> RuffusRuffusTask:
#         inst: RuffusTask = RuffusTask(param)
#
#         formatter = format_input(input_suffixes)
#
#         output_dir = OUTPUT_DIR / ''.join([classname(task) + self._idx])
#         output_format = []
#         for i, suffix in enumerate(output_suffixes):
#             output_format.append(str(
#                 output_dir / ''.join(['{PREFIX[', str(i), ']},', suffix])))
#
#         if 'Transform' in parent_names(task):
#             self._pipeline.transform(
#                     task_func = inst.run,
#                     filter = output_format(input_suffixes),
#                     input = input_files,
#                     output = output_format)
#
#         elif 'RuffusTask' in parent_names(task):
#             self._pipeline.originate(
#                     task_func = inst.run,
#                     output = )
#
#         self._idx += 1
#
#     def run(self, *args, **kwargs):
#         pass


class RuffusTask(ABC):
    """Abtract class representing a `ruffus` function with output but no input

    The running and task logging of the command is handled here. Extra
    initialization and command construction are handled by subclasses.
    `RuffusFunctions` which inherit from `RuffusTask` will produce data which
    will be input to `Transform`s. In practice, these functions will be used
    for `RuffusFunction`s called during `ruffus.Pipeline.originate` tasks.

    Subclasses are expected to define their own method for constructing their
    command (`_build`). They may also define extra outputs
    (`_add_extra_outputs`) and steps to be run after the main command
    (`_post_command`).


    Attributes
    ----------
    output_type : List[str]
        (Class attribute) Expected output file extensions. If '' then the
        extension of the input is stripped in the output.
    exit_code : List[int]
        The exit codes returned by the task's bash process
    stdout : List[str]
        The standard outputs produced by the task
    stderr : List[str]
        The standard errors produced by the task
    """

    _logger = ROOT_LOGGER
    output_type: List[str] = None

    def __init__(
            self,
            param: List[str] = [],
            log_results: bool = True)-> None:
        """Initialization

        Parameters
        ----------
        log_results: optional
            If True, OS stdout and stderr are logged to the pipeline's logfile
        param: optional
            A list of bash parameters e.g ['--foo', 'bar', '-x', '1']
        """


        # Store parameters
        self._log_results = log_results
        self.param: List[str] = param

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

        # If True, the task takes parameters
        self._parameterized = True

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
        self._add_to_header('RuffusTask', self.__class__.__name__)

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


class Transform(RuffusTask, ABC):
    """A superclass representing a `ruffus` function with output and input

    `RuffusFunctions` which inherit from `Transform` will produce data which
    will be input to further `Transform`s. The vast majority of functions used
    in `ruffus` pipelines fall under the domain of `Transform`s.

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

        if self.param and not self._parameterized:
            raise ParameterError('Cannot add parameters to this task')

    def _define_input_attributes(self, input_files)-> None:
        self.input_files: Sequence[str] = ensure_list(input_files)
        self._input_paths: List[Path] = as_paths(self.input_files)
        self._input_stems = self._get_stems(self._input_paths)
        self._indirs: List[Path] = [f.parent for f in self._input_paths]

    def _log_task_header(self)-> None:
        """Write the task header to the logfile"""
        self._add_to_header('RuffusTask', self.__class__.__name__)
        self._addinput_files_to_header()
        self._addoutput_files_to_header()
        self._logger.write_header(self._task_header)


class SamtoolsIndexFasta(Transform):
    """Index a .fasta file using samtools faidx"""

    input_type = ['fasta']
    output_type = ['fasta.fai']

    def _build(self) -> None:
        """Construct the bash command"""
        # Index the file
        self._add_command([SAMTOOLS, 'faidx', self.input_files, self.param])
        # Move the index to the correct directory
        index_file = add_suffix(self._input_paths[0], '.fai')
        self._add_command(['mv', str(index_file), self.output_files[0]])


class Gunzip(Transform):
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


class Gzip(Transform):
    """Compress a file with gzip"""
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


class PairedBowtie2Align(Transform):
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
            BOWTIE2_BUILD, fasta, bowtie2_index, '>', '/dev/null', '2>&1'])

        # Move the indices to the assembly directory
        self._add_command(['mv', '*.bt2', str(self._indirs[0])])
        bowtie2_index = str(self._indirs[0] / bowtie2_index)

        #  Run Bowtie2 and pipe the output to a sorted bam file
        self._add_command([
            BOWTIE2, '-1', fwd_fastq, '-2', rev_fastq, '-x',
            bowtie2_index, *self.param, '|', SAMTOOLS, 'view', '-bS', '-', '|',
            SAMTOOLS, 'sort', '-o', output_bam, '-'])

        # Index the bam file
        self._add_command([
            SAMTOOLS, 'index', output_bam, '/dev/null', '2>&1'])


class SymlinkInputs(Transform):
    """Create symlinks of the input arguments"""
    input_type = ['ANY']
    output_type = ['SAME']

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command([
            'ln', '-sf', self.param, self.input_files, self.output_files])


class BuildCentrifugeDB(RuffusTask):
    """Build a custom centrifuge database"""
    output_type = ['cf']

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
        self._parameterized = False

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
            CENTRIFUGE_DOWNLOAD, "-o", str(self._taxonomy_dir), "taxonomy"])

        # Download all complete Archaeal, Bacterial, Viral, Fungal, Protozoan
        # and Plant genomes
        self._add_command([
            CENTRIFUGE_DOWNLOAD, "-o", str(self._library_dir), "-m", "-d",
            "archaea,bacteria,viral,fungi,protozoa,plant", "refseq", ">",
            str(self._seqid2taxid)])

        # Download human reference sequence
        self._add_command([
            CENTRIFUGE_DOWNLOAD, "-o", str(self._library_dir), "-d",
            "vertebrate_mammalian", "-a", "Chromosome", "-t", 9606, "-c",
            "reference_genome", ">>", str(self._seqid2taxid)])

        # Concatenate the downloaded sequences
        self._concat_seqs = self._outdirs[0].joinpath("catseq.fna")
        self._add_command([
            "cat", str(self._library_dir / '*' / '*.fna'), ">",
            str(self._concat_seqs)])

        # Build the database
        self._add_command([
            CENTRIFUGE_BUILD, "--conversion-table", str(self._seqid2taxid),
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


class Centrifuge(Transform):
    """Run centrifuge"""
    input_type = ['fasta', 'cf']
    output_type = ['tsv']

    def _build(self)-> None:
        """Build the command line arguments"""
        assembly = self.input_files[0]
        centrifuge_idx = self.input_files[1]
        self._add_command([CENTRIFUGE, '-k', '1', '-f', '-x',
            centrifuge_idx, '-U', assembly, '--report-file',
            self.output_files[0], self.param])


class ConvertCentrifugeToHits(Transform):
    """Convert a centrifuge output file to a hits file"""
    input_type = ['centrifuge_output']
    output_type = ['tsv']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parameterized = False

    def _build(self)-> None:
        """Build the command line arguments"""
        self._add_command([
            'centrifuge_to_hits.py', self.input_files[0],
            self.output_files[0]])


#  class BlobtoolsCreate(Transform):
#      """Run centrifuge"""
#      input_type = ['fasta', 'hits', 'fastq', 'fastq']
#      output_type = ['tsv']
#
#      def _build(self)-> None:
#          """Build the command line arguments"""
#          assembly = self.input_files[0]
#          centrifuge_idx = self.input_files[1]
#          self._add_command(['centrifuge', '-k', '1', '-f', '-x',
#              centrifuge_idx, '-U', assembly, '--report-file',
#              self.output_files[0], self.param])


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


"""Type variable for Transform-like function"""
TransformFunction = Callable[[Sequence[str], Sequence[str], Any], None]

def apply_(task: Type[Transform])-> TransformFunction:
    """Return function which applies a `Transform` to multiple inputs

    All extra parameters are passed to `task`

    Parameters
    ----------
    task
        A Transform class

    Returns
    -------
    Callable[[Sequence[str], Sequence[str], Any], None]
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


def output_format(suffixes: List[str], dirname: str = '')-> List[str]:
    """Produce `ruffus` output format strings

    The format specifies that one output file should be placed in `dirname` for
    each suffix in `suffixes`. The return is intended to be passed as the
    output argument to `ruffus.Pipeline.RuffusTask`.

    Parameters
    ----------
    suffixes
        List of suffixes to append to the output files.
    dirname
        The output directory basename

    Returns
    -------
    List[str]
        A list of format strings
    """
    output_format: List[str] = []
    for suffix in suffixes:
        output_format.append(str(
            OUTPUT_DIR / dirname / ''.join(['{PREFIX[0]}', suffix])))

    return output_format


def format_input(suffixes: List[str])-> formatter:
    """Create a ruffus formatter which splits input files into three parts

    The return is intended to be passed as the filter argument to
    `ruffus.Pipeline.RuffusTask`.

    Example
    -------
    format_(['fasta', 'fa'], ['fastq', 'fq'])
        -> formatter(
                ".*/(?P<PREFIX>[^\.]*)\.(?P<MID>.*)\.?(?P<SUFFIX>fasta|fa|),
                ".*/(?P<PREFIX>[^\.]*)\.(?P<MID>.*)\.?(?P<SUFFIX>fastq|fq|))

    Parameters
    ----------
    suffixes
        List of lists of accepted suffixes.

    Returns
    -------
    formatter
        A ruffus formatter with PREFIX, MID, and SUFFIX.
    """
    regexes = []
    regexes.append(r".*/(?P<PREFIX>[^.]*).*") # Capture Prefix

    # We reverse sort so that 'fq.gz' is matched before 'fq'. Without the
    # sort, only the first part of two part suffixes would be matched
    suffixes = sorted(suffixes, reverse=True)
    suff_str = '|'.join(suffixes)
    regexes.append(''.join([r'.*/.*', suff_str]))

    input_formatter = formatter(*regexes)
    return input_formatter


# -----------------------------------------------------------------------------
# Exceptions
# -----------------------------------------------------------------------------


class ParameterError(ValueError):
    """Thrown when unparameterized task is given parameters"""
    pass
