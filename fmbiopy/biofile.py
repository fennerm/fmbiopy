"""
Classes for storing and validating various types of bioinformatics
files.

Most common use case will be for pipeline scripts, which benefit from strict
validation of files between analysis steps.

Usage
-----
fastq = Fastq(<list of .fastq file paths>)
sam = Sam(<list of .sam file paths>)

Notes
-----
Files need not necessarily exist upon initialization, but they must exist upon
the first access attempt. In this way, a file group class can act as a promise
that files will exist later.

Assumes that files are named using a dot based naming scheme with the unique ID
of the file before the first dot:
    <unique_ID>.<processing_step>.<...>.<extension>
E.g FA_SC.trim.map.sam
This could be changed with a little work, but it would probably be easier to
just rename the files.

FileGroups are designed to be initialized and then used without editing. Once
a group has been initialized, it is not recommended to attempt to change the
files it maps to.
"""

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmstring as fmstring
import os
from typing import List
from typing import Sequence
from typing import Union


class BioFileGroup(object):
    """Superclass for storing and validating groups of bioinformatics files.

    Classes of more specific filetypes inherit the majority of their attributes
    and validation methods from it, however there should be no reason to
    initialize a BioFileGroup itself. Rather a new type of filegroup should
    inherit from it.

    Notes
    -----
    Input files to BioFileGroups are automatically sorted alphabetically.

    Parameters
    ----------
    files :
        The list of files to be stored in the group. Files must be the same
        type and non-empty (unless possibly_empty==True).
    gzipped :
        If True, files should be gzipped and have the .gz extension. Detecting
        that files are gzipped is obviously trivial, this flag is just used to
        ensure that the output function knows the gzip state. If this flag is
        unset but the files are gzipped, this indicates that an unintended zip
        step has occured or vice versa.

    Attributes
    ----------
    gzipped : bool
        Same as parameter
    names : List[str]
        The unique IDs of each file.
    """

    def __init__(
            self,
            files: Sequence[str],
            gzipped: bool = False,
            possibly_empty: bool = False
            ) -> None:

        self._paths = files
        self.validated = False
        self._accepted_extensions = None
        self._possibly_empty = possibly_empty
        self.gzipped = gzipped
        self._extensions = self._get_extensions()

        self._prevalidate()

        # These steps are done after prevalidation to prevent unexpected errors
        self._paths = [os.path.abspath(f) for f in files]
        self._paths = sorted(files)
        self.names = self._get_names()

    def validate(self) -> None:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt
        """
        self._check_extensions()
        self._check_files_non_empty()
        self.validated = True

    def _get_extensions(self) -> List[str]:
        """Get the file extensions of the files.

        Returns
        -------
        If gzipped then it returns the two part extension E.g fq.gz. Otherwise
        just returns the final extension.
        """

        # If gzipped we return the two part extension E.g 'fq.gz'
        if self.gzipped:
            return [fmpaths.last_two_suffix(f) for f in self._paths]
        else:
            return [fmpaths.get_final_suffix(f) for f in self._paths]

    def _get_names(self) -> List[str]:
        """Get the names of the files in the group

        Returns
        -------
        The substring of the paths' basenames before the first dot. E.g
        'FA_SC.trim.map.sam' -> 'FA_SC'
        """
        basenames = [os.path.basename(path) for path in self._paths]
        names = [fmpaths.get_prefix(basename) for basename in basenames]
        return names

    def __getitem__(self, item) -> str:
        """Get a file from the group"""

        if not self.validated:
            self.validate()

        return self._paths[item]

    def _prevalidate(self) -> None:
        """Do some basic input validation upon initialization"""

        self._check_paths_not_none()
        self._check_paths_not_string()
        self._check_gzip()

    def _check_files_non_empty(self) -> None:
        if not self._possibly_empty:
            for path in self._paths:
                if os.path.getsize(path) < 3:
                    raise ValueError(path + " is empty")

    def _check_paths_not_none(self) -> None:
        """Check paths not an empty list"""
        if not self._paths:
            raise ValueError('Empty paths in BioFileGroup')

    def _check_paths_not_string(self) -> None:
        """Check that paths is not a single string"""

        if isinstance(self._paths, str):
            raise TypeError("""
            Input to BioFileGroup is a string (expects list)
            """)

    def _check_extensions(self) -> None:
        """Check that the file extensions match the accepted extensions

        The superclass should not have accepted extensions, but subclasses will
        use this method for extension validation.
        """
        if self._accepted_extensions:
            for extension in self._extensions:
                # Extension check is not caps sensitive
                if not extension.lower() in self._accepted_extensions:
                    raise ValueError("Invalid file extension")

    def _check_gzip(self) -> None:
        """Check that the gzipped flag parameter is correct"""

        if self.gzipped:
            if not fmcheck.any_endswith(self._extensions, 'gz'):
                raise TypeError("Files are not gzipped")
        else:
            if fmcheck.any_endswith(self._extensions, 'gz'):
                raise TypeError("Files are gzipped, but gzipped flag not set")

    def __len__(self) -> int:
        """Length of file list"""
        return len(self._paths)


class FastqGroup(BioFileGroup):
    """BioFileGroup for holding Fastq files.

    Parameters inherited from superclass
    """

    def __init__(self, *args, **kwargs) -> None:

        super().__init__(*args, **kwargs)

        self._accepted_extensions = ['fastq', 'fq']
        if self.gzipped:
            self._accepted_extensions += ['fastq.gz', 'fq.gz']


class FastaGroup(BioFileGroup):
    """BioFileGroup for holding .fasta files."""

    def __init__(self, *args, **kwargs) -> None:
        """Initialization"""
        super().__init__(*args, **kwargs)
        self._accepted_extensions = ['fa', 'fasta', 'mfa']


class SamGroupGroup(BioFileGroup):
    """BioFileGroup for holding .sam files"""

    def __init__(self, *args, **kwargs) -> None:
        """Initialization """
        super().__init__(*args, **kwargs)
        self._accepted_extensions = ['sam']


class BamGroup(BioFileGroup):
    """BioFileGroup for holding .bam files"""

    def __init__(self, *args, **kwargs) -> None:
        """Initialization """
        super().__init__(*args, **kwargs)
        self._accepted_extensions = ['bam']


class SamtoolsFAIndexGroup(BioFileGroup):
    """BioFileGroup for holding samtools .fai files"""
    def __init__(self, *args, **kwargs):
        """Initialization """
        super().__init__(*args, **kwargs)
        self._accepted_extensions = ['fai']


class Bowtie2IndexGroup(BioFileGroup):
    """BioFileGroup for holding Bowtie2 .bt2 fasta index files"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._accepted_extensions = ['bt2']
        self._index_prefixes = self._get_index_prefixes()

        # The names needs to be set differently for Bowtie2 indices because
        # each fasta file has multiple indices, but Bowtie2 refers to the group
        # of indices by a single name
        self.names = self._get_names()

    @property
    def index_files(self) -> List[str]:
        """Simple getter function to retrieve the paths of all index files."""

        return self._paths()

    def _get_index_prefixes(self) -> List[str]:
        """ """
        names = []

        # Get the prefixes of each .bt2 index
        for path in self._paths:
            # Reverse bowtie2 indices need three extensions removed, the rest
            # just need two
            if '.rev.' in os.path.basename(path):
                names.append(fmpaths.remove_suffix(path, 3))
            else:
                names.append(fmpaths.remove_suffix(path, 2))

        # Get the unique elements
        names = fmstring.get_unique(names)

        return names

    def __getitem__(self, item) -> str:
        """Get an index prefix from the list """
        if not self.validated:
            self.validate()

        return self._index_prefixes[item]


"""
TypeVar for typing module which groups different types of fasta index classes
into a single type
"""
FastaIndexGroup = Union[SamtoolsFAIndexGroup, Bowtie2IndexGroup]


class PairedFastqGroup(object):
    """Stores two groups of paired FastqGroup files

    Attributes
    ----------
    fwd - FastqGroup, rev - FastqGroup
        Fastq objects to be paired

    """
    def __init__(
            self,
            forward_fastq: FastqGroup,
            reverse_fastq: FastqGroup
            ) -> None:
        """Initialization

        Parameters
        ----------
        forward_fastq - FastqGroup, reverse_fastq - FastqGroup
            FastqGroup objects to be paired
        """
        self.fwd = forward_fastq
        self.rev = reverse_fastq
        self._prevalidate()

    def _prevalidate(self) -> None:
        """Check that the PairedFastqGroup is valid"""
        self.__check_lengths()

    def __check_lengths(self) -> None:
        """Raise exception if lengths of the two Fastq objects are not equal"""

        if len(self.fwd) != len(self.rev):
            raise ValueError("FastqGroup objects are not the same length")

    def __len__(self) -> int:
        """Length of the PairedFastqGroup """
        return len(self.fwd)

    def __getitem__(self, item) -> List[FastqGroup]:
        """Get the pair of Fastq files indexed by 'item' """
        return list([self.fwd[item], self.rev[item]])


class IndexedFastaGroup(object):
    """Represents a Fasta file grouped with its indices"""
    def __init__(
            self,
            fasta: FastaGroup,
            *indices: FastaIndexGroup
            ) -> None:
        self.fasta = fasta
        self.indices = indices

    def __getitem__(self, item) -> List[str]:
        index_items = [index[item] for index in self.indices]
        return [self.fasta[item]] + index_items
