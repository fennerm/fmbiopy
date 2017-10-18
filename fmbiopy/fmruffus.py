""" Classes and functions for use in Ruffus pipelines

        cat library/*/*.fna > input-sequences.fna

Ruffus: http://www.ruffus.org.uk/
"""

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


"""The `RuffusLog` shared across all `RuffusTasks`"""
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
    log_results: optional
        If True, OS stdout and stderr are logged to the pipeline's logfile
    param : Optional
        A list of bash parameters e.g ['--foo', 'bar', '-x', '1']

    Other Parameters
    ----------------
    debug : Optional
        If True, run in debugging mode

    Attributes
    ----------
    input_type: str
        (Class attribute) Expected input file extension. If None then any file
        type is acceptable input.
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
    input_type: List[str] = None
    output_type: List[str] = None
    _logger = ROOT_LOGGER

    def __init__(self,
                 input_files: Sequence[str],
                 output_files: Sequence[str],
                 log_results: bool = True,
                 param: Sequence[str] = None,
                 run_on_init: bool = True) -> None:

        # Store parameters
        self._input_files: Sequence[str] = self._ensure_list(input_files)
        self._input_paths: List[Path] = as_paths(self._input_files)
        self._output_files: Sequence[str] = self._ensure_list(output_files)
        self._output_paths: List[Path] = as_paths(self._output_files)
        self._log_results = log_results
        self._param = param

        # Initialize attributes
        self._init_attributes()

        # Check inputs
        self._check_inputs()

        # Create output directories if necessary
        self._mkdirs()

        # Construct command
        self._build()

        # Run command
        if run_on_init:
            self._run()

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

    def _build(self)-> None:
        """Subclasses construct their own tasks"""
        pass

    def _check_inputs(self)-> None:
        """Check that input and output files are valid"""
        for inp in self._input_paths:
            if inp in self._output_paths:
                raise FileExistsError('Attempted to overwrite input file(s)')

    def _cleanup(self)-> None:
        """Cleanup"""
        remove_all(self._output_paths, silent=True)
        remove_all(self._extra_outputs, silent=True)

    def _ensure_list(self, x: Sequence[str])-> Sequence[str]:
        """If list, return list. If str, return str"""
        if isinstance(x, str):
            return [x]
        return x

    def _get_stems(self, paths: List[Path])-> List[str]:
        """Get the filename stems (everything but the extension)"""
        stems = []
        for path in paths:
            if path.suffix == '.gz':
                stems.append(Path(path.stem).stem)
            else:
                stems.append(path.stem)
        return stems

    def _init_attributes(self):
        """Initialize empty attributes"""
        # Bash command as list
        self._command: List = []

        # Returns from command(s)
        self.exit_code: List[int] = []
        self.stdout: List[str] = []
        self.stderr: List[str] = []

        # Output directory
        self._outdirs: List[Path] = [f.parent for f in self._output_paths]

        # Input directory
        self._indirs: List[Path] = [f.parent for f in self._input_paths]

        # Input and output file stems
        self._input_stems = self._get_stems(self._input_paths)
        self._output_stems = self._get_stems(self._output_paths)

        # Extra output files which are not specified by the user.
        self._extra_outputs: List[Path] = self._add_extra_outputs()

        # These attributes may or may not be set by subclasses:
        # -----------------------------------------------------

        # If True, command run directly in shell
        if not hasattr(self, '_shell'):
            self._shell: bool = False

        # If True, the input files are deleted once the output files are
        # produced.
        if not hasattr(self, '_inplace'):
            self._inplace: bool = False

    def _mkdirs(self)-> None:
        """Create output directories if necessary"""
        for d in self._outdirs:
            d.mkdir(parents=True, exist_ok=True)

    def _post_task(self)-> None:
        """Commands to be run after `_run`"""
        pass

    def _run(self)-> None:
        """Run a command set with `_build`"""
        try:
            if self._log_results:
                header = ''.join([
                    'Task: ', self.__class__.__name__, '\n', 'Input_files: ',
                    ', '.join(self._input_files), '\n', 'Output_files: ',
                    ', '.join(self._output_files)])
                self._logger.write_header(header)

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

            if self._inplace:
                remove_all(self._input_paths, silent=True)

            self._post_task()

        except Exception:
            self._cleanup()
            raise

    def _store_results(self, results: Tuple[int, str, str])-> None:
        """Update attributes with the results of `fmsystem.run_command`"""
        self.exit_code.append(results[0])
        self.stdout.append(results[1])
        self.stderr.append(results[2])


class SamtoolsIndexFasta(RuffusTask):
    """Index a .fasta file using samtools faidx"""

    input_type = ['fasta']
    output_type = ['fasta.fai']

    def _build(self) -> None:
        """Construct the bash command"""
        self._add_command(['samtools', 'faidx', self._input_files])


class Gunzip(RuffusTask):
    """Unzip a file with gunzip"""
    input_type = ['gz']
    output_type = ['']

    def __init__(self, *args, **kwargs)-> None:
        """Initialize"""
        self._shell = True
        self._inplace = True
        super().__init__(*args, **kwargs)

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command(['gunzip', '-c', self._param, self._input_files, '>',
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

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command([
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
        self._shell = True
        super().__init__(*args, **kwargs)


    def _add_extra_outputs(self) -> List[Path]:
        """Add the .bt2 and .bai files to the output file list"""
        bowtie2_indices = get_bowtie2_indices(self._input_stems[0])
        bai_file: Path = add_suffix(self._output_paths[0], '.bai')
        return flatten([bowtie2_indices, bai_file])

    def _build(self) -> None:
        """Construct the bash command"""
        fasta = self._input_files[0]
        fwd_fastq = self._input_files[1]
        rev_fastq = self._input_files[2]
        output_bam = self._output_files[0]

        # Index fasta first
        bowtie2_index = self._input_stems[0]
        self._add_command([
            'bowtie2-build', fasta, bowtie2_index, '>', '/dev/null', '2>&1'])

        # Move the indices to the assembly directory
        self._add_command(['mv', '*.bt2', str(self._indirs[0])])

        #  Run Bowtie2 and pipe the output to a sorted bam file
        self._add_command([
            'bowtie2', '-1', fwd_fastq, '-2', rev_fastq, '-x', bowtie2_index,
            '|', 'samtools', 'view', '-bS', '-', '|', 'samtools', 'sort',
            '-f', '-', output_bam])

        # Index the bam file
        self._add_command([
            'samtools', 'index', output_bam, '/dev/null', '2>&1'])


class SymlinkInputs(RuffusTask):
    """Create symlinks of the input arguments"""
    input_type = ['ANY']
    output_type = ['SAME']

    def _build(self)-> None:
        """Construct the bash command"""
        self._add_command([
            'ln', '-sf', self._input_files, self._output_files])


class BuildCentrifugeDB(RuffusTask):
    """Build a custom centrifuge database"""
    def __init__(self, *args, **kwargs)-> None:
        self._shell = True
        super().__init__(*args, **kwargs)

    def _add_extra_outputs(self)-> List[Path]:
        """Add the library, taxonomy and SeqID2Taxmap files to extra outputs"""
        # The map from seqid2taxid
        self._seqid2taxid = self._outdirs[0].joinpath('seqid2taxid.map')

        # The index files
        centrifuge_idx = [
                ''.join([self._output_files[0], '.', i, '.cf'])
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
            self._output_files[0]])

    def _post_task(self)-> None:
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


"""Type variable for a function with same inputs and outputs as a RuffusTask"""
TaskFunction = Callable[[Sequence[str], Sequence[str], Any], None]


def apply_(task: Type[RuffusTask])-> TaskFunction:
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

