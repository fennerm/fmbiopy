"""Classes for storing and validating various types of bioinformatics files.

Most common use case will be for pipeline scripts, which benefit from strict
validation of files between analysis steps.

Alternatively the classes can simply be used to validate groups of files for
function arguments etc. E.g if FastaGroup[list_of_fasta_files] runs without
error, you have some confidence that these files can be passed to analysis
modules for processing without error. However file validation is targeted
mostly at catching user error - e.g passing a Sam file rather than a Bam file,
it does almost no internal checks for file format correctness. Consequently,
errors introduced to files during analysis will not be caught.

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
import collections
from itertools import compress
import os
from typing import List
from typing import Sequence
from typing import Union

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths


class BioFileValidationError(Exception):
    """Basic exception for errors raised by BioFile validation checks

    Parameters
    ----------
    filenames
        The files which caused the error

    Attributes
    ----------
    msg
        The formatted error message
    """
    def __init__(self, filenames : Sequence[str] = []) -> None:
        self.filenames = filenames
        self._construct_msg()
        super().__init__(self.msg)

    def _formatted_filenames(self) -> str:
        """Construct the part of the error message which lists files"""
        return '\n'.join(['Files:', ', '.join(self.filenames)])

    def _err_description(self) -> str:
        """Construct the error description to output"""
        return ''

    def _construct_msg(self) -> None:
        """Construct the combined error message"""
        self.msg : str = '\n'.join([
                self._formatted_filenames(),
                self._err_description()])


class MatchedPrefixGroupValidationError(BioFileValidationError):
    """Basic exception for errors raised by `MatchedPrefixGroup` validation"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)

    def _formatted_filenames(self) -> None:
        """Format the file group names for printing"""
        formatted_filenames = 'Groups:\n'

        for fs in self.filenames:
            formatted_filenames += '[' + ', '.join(fs) + ']\n'
        return formatted_filenames


class EmptyFileError(BioFileValidationError):
    """Exception raised when one or more of the stored files are empty"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)

    def _err_description(self) -> str:
        return "File is empty but possibly_empty is False"


class GzipStatusError(BioFileValidationError):
    """Exception raised when one or more of the stored files are empty"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)

    def _err_description(self) -> str:
        return "Gzip status of files does not match the gzip argument"


class FileExtensionError(BioFileValidationError):
    """Exception raised when files do not have the expected file extension"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)

    def _err_description(self) -> str:
        return "Unexpected file extension"


class DuplicateFilegroupError(MatchedPrefixGroupValidationError):
    """Exception raised when a `BioFileGroup` is matched with itself"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)


class PrefixMatchError(MatchedPrefixGroupValidationError):
    """Exception raised the files do not share a prefix"""
    def __init__(self, filenames) -> None:
        super().__init__(filenames)


class GroupLengthError(MatchedPrefixGroupValidationError):
    """Exception raised `MatchedPrefixGroups`s don't have equal lengths"""

    def __init__(self, filenames) -> None:
        super().__init__(filenames)


class BioFileGroup(collections.abc.Sized):
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
    names : typing.List[str]
        The unique IDs of each file.
    paths : typing.List[str]
        The paths of the files
    """

    def __init__(
            self,
            files: Sequence[str],
            gzipped: bool = False,
            possibly_empty: bool = False) -> None:

        # Store paramaters
        self.paths = files
        self.gzipped = gzipped
        self.validated = False

        # List of acceptable extensions for files in BioFileGroup
        self._accepted_extensions : typing.List[str] = None

        # True if it is acceptable for the files to be empty
        self._possibly_empty = possibly_empty

        # Get the actual file extensions
        self._extensions = self._get_extensions()

        self._prevalidate()

        # These steps are done after prevalidation to prevent unexpected errors
        self.paths = [os.path.abspath(f) for f in files]
        self.paths = sorted(files)
        self.names = self._get_names()

    def __getitem__(self, item) -> str:
        """Get a file from the group"""

        if not self.validated:
            self.validate()

        return self.paths[item]

    def __eq__(self, other) -> bool:
        """Test for BioFileGroup is equal to another"""

        # Cannot be equal if lengths are different
        if len(self) != len(other):
            return False

        # Test that all paths are the same
        for me, you in zip(self, other):  # type: ignore
            if me != you:
                return False

        return True

    def __len__(self) -> int:
        """Length of file list"""
        return len(self.paths)

    def validate(self) -> None:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt
        """
        self._check_extensions()
        self._check_files_non_empty()
        self.validated = True

    def _check_paths_not_none(self) -> None:
        """Check paths not an empty list"""
        if not self.paths:
            raise ValueError('Empty paths in BioFileGroup')

    def _check_paths_not_string(self) -> None:
        """Check that paths is not a single string"""

        if isinstance(self.paths, str):
            raise TypeError("""
            Input to BioFileGroup is a string (expects list)""")

    def _get_extensions(self) -> List[str]:
        """Get the file extensions of the files.

        Returns
        -------
        If gzipped then it returns the two part extension E.g fq.gz. Otherwise
        just returns the final extension.
        """

        # If gzipped we return the two part extension E.g 'fq.gz'
        if self.gzipped:
            return [fmpaths.last_two_suffix(f) for f in self.paths]
        else:
            return [fmpaths.get_final_suffix(f) for f in self.paths]

    def _get_names(self) -> List[str]:
        """Get the names of the files in the group

        Returns
        -------
        The substring of the paths' basenames before the first dot. E.g
        'FA_SC.trim.map.sam' -> 'FA_SC'
        """
        basenames = [os.path.basename(path) for path in self.paths]
        names = [fmpaths.get_prefix(basename) for basename in basenames]
        return names

    def _prevalidate(self) -> None:
        """Do some basic input validation upon initialization"""

        self._check_paths_not_none()
        self._check_paths_not_string()
        self._check_gzip()

    def _check_files_non_empty(self) -> None:
        if not self._possibly_empty:
            empty = []
            for path in self.paths:
                if os.path.getsize(path) < 3:
                    empty.append(path)
            if len(empty) > 0:
                raise EmptyFileError(empty)

    def _check_extensions(self) -> None:
        """Check that the file extensions match the accepted extensions

        The superclass should not have accepted extensions, but subclasses will
        use this method for extension validation.
        """
        if self._accepted_extensions:
            incorrect_extensions = []
            for i, extension in enumerate(self._extensions):
                # Extension check is not caps sensitive
                if not extension.lower() in self._accepted_extensions:
                    incorrect_extensions.append(self.paths[i])

            if len(incorrect_extensions) > 1:
                raise FileExtensionError([self.paths[i]])

    def _check_gzip(self) -> None:
        """Check that the gzipped flag parameter is correct"""

        which_gzipped = [ext.endswith('gz') for ext in self._extensions]
        is_gzipped = list(compress(self.paths, which_gzipped))
        which_not_gzipped = [not x for x in is_gzipped]
        not_gzipped = list(compress(self.paths, which_not_gzipped))

        if self.gzipped:
            if any(which_not_gzipped):
                raise GzipStatusError(not_gzipped)
        else:
            if any(which_gzipped):
                raise GzipStatusError(is_gzipped)


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
    def index_files(self) -> Sequence[str]:
        """Simple getter function to retrieve the paths of all index files."""
        return self.paths

    def _get_index_prefixes(self) -> List[str]:
        """Get the prefixes of the bowtie2 indices"""
        names = []

        # Get the prefixes of each .bt2 index
        for path in self.paths:
            # Reverse bowtie2 indices need three extensions removed, the rest
            # just need two
            if '.rev.' in os.path.basename(path):
                names.append(fmpaths.remove_suffix(path, 3)[0])
            else:
                names.append(fmpaths.remove_suffix(path, 2)[0])

        # Get the unique elements
        names = fmlist.get_unique(names)

        return names

    def __getitem__(self, item) -> str:
        """Get an index prefix from the list """
        if not self.validated:
            self.validate()

        return self._index_prefixes[item]


"""
typing module TypeVar which groups different types of fasta index classes
into a single type
"""
FastaIndexGroup = Union[SamtoolsFAIndexGroup, Bowtie2IndexGroup]


class MatchedPrefixGroup(object):
    """Stores groups of matched BioFileGroups with same prefix

    Prefixes are checked for equality. Groups are checked for equal length.

    Parameters
    ----------
    groups
        typing.List of matched BioFileGroups

    Attributes
    ----------
    groups : typing.List[BioFileGroup]
        Same as Parameters """
    def __init__(self, groups: List[BioFileGroup]) -> None:
        self.groups = groups
        self._prevalidate()

    def _prevalidate(self) -> None:
        """Check that the MatchedPrefixGroup is valid"""
        self.__check_files_not_same()
        self.__check_lengths_match()
        self.__check_same_file_prefix()

    def __check_files_not_same(self) -> None:
        """Check that none of the filegroups are the exact same"""
        for i, group1 in enumerate(self.groups):
            for j, group2 in enumerate(self.groups):
                if i != j and group1 == group2:
                    raise DuplicateFilegroupError(self.groups)

    def __check_lengths_match(self) -> None:
        group_lengths = [len(g) for g in self.groups]
        if not fmcheck.all_equal(group_lengths):
            raise GroupLengthError(self.groups)

    def __check_same_file_prefix(self) -> None:
        """Check that the stored BioFileGroups all have the same prefixes"""
        names = [group.names for group in self.groups]
        if not fmcheck.all_equal(names):
            raise PrefixMatchError(self.groups)

    def __len__(self) -> int:
        """Length of the `MatchedPrefixGroup`"""
        return len(self.groups[0])

    def __getitem__(self, item) -> List[str]:
        """Index the `MatchedPrefixGroup`"""
        return [g[item] for g in self.groups]


class PairedFastqGroup(MatchedPrefixGroup):
    """Stores two groups of paired FastqGroup files

    Parameters
    ----------
    forward_fastq, reverse_fastq
        Fastq objects to be paired

    Attributes
    ----------
    forward_fastq - FastqGroup, reverse_fastq - FastqGroup
        Same as Parameters

    """
    def __init__(
            self,
            forward_fastq: FastqGroup,
            reverse_fastq: FastqGroup) -> None:
        self.fwd = forward_fastq
        self.rev = reverse_fastq
        self.groups = [self.fwd, self.rev]
        super()._prevalidate()


class IndexedFastaGroup(object):
    """Represents a Fasta file grouped with its indices"""
    def __init__(
            self,
            fasta: FastaGroup,
            *indices: FastaIndexGroup) -> None:
        self.fasta = fasta
        self.indices = indices

    def __getitem__(self, item) -> List[str]:
        index_items = [index[item] for index in self.indices]
        return [self.fasta[item]] + index_items
