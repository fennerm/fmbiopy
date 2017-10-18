"""Classes for storing and validating various types of bioinformatics files.

Most common use case will be for pipeline scripts, which benefit from strict
validation of files between analysis steps.

Alternatively the classes can simply be used to validate groups of files for
function arguments etc. E.g if `BiofileGroup`[list_of_fasta_files,
type='fasta'] runs without error, you have some confidence that these files can
be passed to analysis modules for processing without error. However file
validation is targeted mostly at catching user error - e.g passing a `Sam` file
rather than a `Bam` file, it does almost no internal checks for file format
correctness. Consequently, errors introduced to files during analysis will not
be caught.

Usage
-----
Fastq(<fastq_file>)
Sam(<sam_file>)
BiofileGroup([file, ...], type=<type>)

Notes
-----
Files need not necessarily exist upon initialization, but they must exist upon
the first access attempt. In this way, a file group class can act as a promise
that files will exist later.

of the file before the first dot:
    <unique_ID>.<processing_step>.<...>.<extension>
E.g FA_SC.trim.map.sam
This could be changed with a little work, but it would probably be easier to
just rename the files.

FileGroups are designed to be initialized and then used without editing. Once
a group has been initialized, it is not recommended to attempt to change the
files it maps to.
"""
from collections.abc import Sized
from os import path
from pathlib import (
        Path,
        PurePath,
        )
from typing import (
        List,
        Sequence,
        Type,
        )

from fmbiopy.fmcheck import all_equal
from fmbiopy.fmpaths import (
        as_strs,
        prefix,
        )

# -----------------------------------------------------------------------------
# Biofile Objects
# -----------------------------------------------------------------------------


class Biofile():
    """Superclass for storing and validating bioinformatics files.

    Classes of more specific filetypes inherit the majority of their attributes
    and validation methods from it.

    Parameters
    ----------
    path
        File to be stored in the group (need not exist)
    gzipped, Optional
        If True, file should be gzipped and have the .gz extension. Detecting
        that files are gzipped is obviously trivial, this flag is just used to
        ensure that the output function knows the gzip state. If this flag is
        unset but the files are gzipped, this indicates that an unintended zip
        step has occured or vice versa.
    possibly_empty, Optional
        If True, empty files will not cause validation errors.

    Attributes
    ----------
    input_type : str
        (class) Filetype which is stored
    accepted_extensions : List[str]
        (class) List of acceptable input file extensions
    gzipped : bool
        Same as parameter
    """
    input_type: str = 'ANY'
    accepted_extensions: List[str] = ['ANY']

    def __init__(
            self,
            path: PurePath,
            gzipped: bool = False,
            possibly_empty: bool = False) -> None:

        # Store paramaters
        self._path = path
        self.gzipped = gzipped
        self.validated = False

        # True if it is acceptable for the files to be empty
        self.possibly_empty = possibly_empty

        # Get the actual file extension
        self.extension = self._get_extension()

        self._prevalidate()


    @property
    def path(self) -> Path:
        """The path to the stored file"""
        if not self.validated:
            self.validate()
        return Path(self._path)

    @property
    def name(self) -> str:
        """The basename of the stored file"""
        if not self.validated:
            self.validate()
        return self._path.name

    def validate(self) -> bool:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt

        Returns
        -------
        True if successful
        """
        self._check_file_not_empty()
        self.validated = True
        return True

    def _check_not_dir(self) -> bool:
        """Check that input is not a directory

        Returns
        -------
        True if successful
        """
        if Path(self._path).is_dir():
            raise TypeError('File cannot be a directory')
        return True

    def _get_extension(self) -> str:
        """Get the file extension

        Returns
        -------
        If gzipped then it returns the two part extension E.g fq.gz. Otherwise
        just returns the final extension.
        """

        # If gzipped we return the two part extension E.g 'fq.gz'
        if self.gzipped:
            return ''.join(self._path.suffixes[-2:])
        return self._path.suffix

    def _prevalidate(self) -> None:
        """Do some basic input validation upon initialization"""

        self._check_not_dir()
        self._check_gzip()
        self._check_extension()

    def _check_file_not_empty(self) -> bool:
        """Check that the file has contents"""
        if not self.possibly_empty:
            if path.getsize(self._path) < 3:
                raise EmptyFileError(self._path)
        return True

    def _check_extension(self) -> bool:
        """Check that the file extension matches the accepted extensions

        The superclass should not have accepted extensions, but subclasses will
        use this method for extension validation.
        """
        if self.accepted_extensions != ['ANY']:
            # Extension check is not caps sensitive
            if self.extension.lower() not in self.accepted_extensions:
                raise FileExtensionError(self._path)
        return True

    def _check_gzip(self) -> bool:
        """Check that the gzipped flag parameter is correct"""

        if self.gzipped:
            if '.gz' not in self.extension:
                raise GzipStatusError(self._path)
        else:
            if 'gz' in self.extension:
                raise GzipStatusError(self.path)
        return True

    def __eq__(self, other) -> bool:
        """Test for BiofileGroup is equal to another"""
        return self.path == other.path


class Fastq(Biofile):
    """Biofile class for holding .fastq files."""

    input_type = 'fastq'
    accepted_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']


class Fasta(Biofile):
    """Biofile for class holding .fasta files."""

    input_type = 'fasta'
    accepted_extensions = ['.fasta', '.fa', '.fasta.gz', '.fa.gz']


# class Sam(Biofile):
#     """Biofile class for holding .sam files"""
#
#     input_type = 'sam'
#     accepted_extensions = ['.sam']
#
#
# class Bam(Biofile):
#     """Biofile class for holding .bam files"""
#
#     input_type = 'bam'
#     accepted_extensions = ['.bam']
#
#     def __init__(self, *args, **kwargs) -> None:
#         """Initialization """
#         super().__init__(*args, **kwargs)
#         self._accepted_extensions = ['bam']
#
#
# class SamtoolsFAIndex(Biofile):
#     """Biofile class for holding samtools .fai files"""
#     input_type = 'fai'
#     accepted_extension = ['.fai']


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def type_to_class(typ: str) -> Type[Biofile]:
    """Convert a type variable to its corresponding Biofile object"""
    if typ:
        try:
            cls = eval(typ.capitalize())
        except NameError:
            cls = Biofile
        return cls
    return None

# -----------------------------------------------------------------------------
# BiofileGroup Objects
# -----------------------------------------------------------------------------


class BiofileGroup(Sized):
    """Superclass for storing and validating groups of bioinformatics files.

    All parameters except files are passed to `Biofile` to initialize each
    individual file.

    Parameters
    ----------
    files
        The list of files to be stored in the group.
    filetype
        The stored filetype

    Attributes
    ----------
    gzipped : bool
        Same as parameter
    names : List[str]
        The unique IDs of each file.
    paths : List[str]
        The paths of the files
    cls : Type[Biofile]
        The stored `Biofile` class
    """
    input_type: List[str] = ['ANY*']
    accepted_extensions = ['ANY']

    def __init__(
            self,
            paths: Sequence[Path],
            filetype: str,
            *args,
            **kwargs,
            ) -> None:

        # Store paramaters
        self._paths = paths
        self.filetype = filetype
        self.validated = False

        self.cls = type_to_class(self.filetype)
        self._biofiles = self._initialize_biofiles(*args, **kwargs)
        self._prevalidate()

    def validate(self) -> None:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt
        """
        for biofile in self._biofiles:
            if not biofile.validated:
                biofile.validate()
        self.validated = True

    @property
    def paths(self) -> Sequence[Path]:
        """The stored filepaths"""
        if not self.validated:
            self.validate()
        return self._paths

    @property
    def biofiles(self) -> Sequence[Biofile]:
        """The group of `Biofile` objects"""
        if not self.validated:
            self.validate()
        return self._biofiles

    def __getitem__(self, item) -> str:
        """Get a file from the group"""

        if not self.validated:
            self.validate()

        return self._paths[item]

    def __eq__(self, other) -> bool:
        """Test for BiofileGroup is equal to another"""

        # Cannot be equal if lengths are different
        if len(self) != len(other):
            return False

        # Test that all paths are the same
        for me, you in zip(self, other):  # type: ignore
            if me != you:
                return False

        return True

    def _initialize_biofiles(self, *args, **kwargs) -> List[Biofile]:
        """Initalize a set of biofiles for the input file list"""
        return [self.cls(p, *args, **kwargs) for p in self._paths]

    def __len__(self) -> int:
        """Length of file list"""
        return len(self._paths)

    def _check_paths_not_none(self) -> bool:
        """Check paths not an empty list"""
        if not self._paths:
            raise ValueError('Empty paths in BiofileGroup')
        return True

    def _check_extensions_same(self) -> bool:
        """Check that the stored file extensions are all the same"""
        extensions = [f.extension for f in self._biofiles]
        if not all_equal(extensions):
            raise FileExtensionsNotSameError(self._paths)
        return True

    def _prevalidate(self) -> bool:
        """Do some basic input validation upon initialization"""

        self._check_paths_not_none()
        self._check_extensions_same()
        return True

#
#
#
# class Bowtie2Indices(BiofileGroup):
#     """Biofile class for holding Bowtie2 .bt2 fasta index files"""
#     input_type = ['bt2']
#     accepted_extension = ['bt2']
#
#     def __init__(self, *args, **kwargs):
#         """Initalize class"""
#         super().__init__()
#         self._names = None
#         self._name = self._get_name()
#
#     @property
#     def name(self) -> Sequence[str]:
#         """Simple getter function to retrieve the paths of all index
#            files."""
#
#         return self.paths
#
#     def _get_name(self) -> str:
#         """Get the prefixes of the bowtie2 indices"""
#
#         # Reverse bowtie2 indices need three extensions removed, the rest
#         # just need two
#         path = os.path.basename(self._paths[0])
#         if '.rev.' in path:
#             return fmpaths.remove_suffix(path, 3)[0]
#         return fmpaths.remove_suffix(path, 2)[0]
#
#     def _check_extensions(self):
#         assert False
#
#     def __getitem__(self, item) -> str:
#         """Get an index prefix from the list """
#         if not self.validated:
#             self.validate()
#
#         return self._index_prefixes[item]
#
#
# """
# module TypeVar which groups different types of fasta index classes
# into a single type
# """
# FastaIndexGroup = Union[SamtoolsFAIndex, Bowtie2Indices]
#
#


class MatchedPrefixGroup(object):
    """Stores groups of matched BiofileGroups with same prefix

    Prefixes are checked for equality. Groups are checked for equal length.

    Parameters
    ----------
    groups
        List of matched BiofileGroups

    Attributes
    ----------
    groups : List[BiofileGroup]
    """
    input_type = ['ANY*']
    accepted_extensions = ['ANY']

    def __init__(self, groups: List[BiofileGroup]) -> None:
        self.groups = groups
        self.validated = False
        self._prevalidate()

    def _prevalidate(self) -> None:
        """Check that the MatchedPrefixGroup is valid"""
        self.__check_files_not_same()
        self.__check_lengths_match()
        self.__check_same_file_prefix()

    def validate(self)-> bool:
        """Validate the `BiofileGroup`s"""
        if not self.validated:
            for group in self.groups:
                group.validate()
        self.validated = True
        return True

    def __check_files_not_same(self) -> None:
        """Check that none of the filegroups are the exact same"""
        for i, group1 in enumerate(self.groups):
            for j, group2 in enumerate(self.groups):
                if i != j and group1 == group2:
                    raise DuplicateFilegroupError(self.groups)

    def __check_lengths_match(self) -> None:
        group_lengths = [len(g) for g in self.groups]
        if not all_equal(group_lengths):
            raise GroupLengthError(self.groups)

    def __check_same_file_prefix(self) -> None:
        """Check that the stored BiofileGroups all have the same prefixes"""
        group_paths = [group._paths for group in self.groups]
        prefixes = []
        for paths in group_paths:
            prefixes.append([prefix(path) for path in paths])
        if not all_equal(prefixes):
            raise PrefixMatchError(self.groups)

    def __len__(self) -> int:
        """Length of the `MatchedPrefixGroup`"""
        return len(self.groups[0])

    def __getitem__(self, item) -> List[str]:
        """Index the `MatchedPrefixGroup`"""
        return [g[item] for g in self.groups]
#
#
# class PairedFastqGroup(MatchedPrefixGroup):
#     """Stores two groups of paired FastqGroup files
#
#     Parameters
#     ----------
#     forward_fastq, reverse_fastq
#         Fastq objects to be paired
#
#     Attributes
#     ----------
#     forward_fastq - FastqGroup, reverse_fastq - FastqGroup
#         Same as Parameters
#
#     """
#
#     input_type = ['fastq', 'fastq']
#
#
# class IndexedFastaGroup(object):
#     """Represents a Fasta file grouped with its indices"""
#     def __init__(
#             self,
#             fasta: FastaGroup,
#             *indices: FastaIndexGroup) -> None:
#         self.fasta = fasta
#         self.indices = indices
#
#     def __getitem__(self, item) -> List[str]:
#         index_items = [index[item] for index in self.indices]
#         return [self.fasta[item]] + index_items


# -----------------------------------------------------------------------------
# Biofile Exceptions
# -----------------------------------------------------------------------------


class BiofileValidationError(Exception):
    """Basic exception for errors raised by Biofile validation checks

    Parameters
    ----------
    name
        The filename which caused the error

    Attributes
    ----------
    msg
        The formatted error message
    """
    def __init__(self, name: PurePath = None) -> None:
        self.name = name
        self.msg = self._construct_msg()
        super().__init__(self.msg)

    def _formatted_filename(self) -> str:
        """Construct the part of the error message which lists the file"""
        return '\n'.join(['File:', '\t' + str(self.name)])

    def _err_description(self) -> str:
        """Construct the error description to output"""
        return ''

    def _construct_msg(self) -> str:
        """Construct the combined error message"""
        return '\n'.join([
            self._formatted_filename(), self._err_description()])


class EmptyFileError(BiofileValidationError):
    """Exception raised when one or more of the stored files are empty"""

    def _err_description(self) -> str:
        return "File is empty but `possibly_empty` is False"


class GzipStatusError(BiofileValidationError):
    """Exception raised when one or more of the stored files are empty"""

    def _err_description(self) -> str:
        return "Gzip status of files does not match the gzip argument"


class FileExtensionError(BiofileValidationError):
    """Exception raised when files do not have the expected file extension"""

    def _err_description(self) -> str:
        return "Unexpected file extension"

# -----------------------------------------------------------------------------
# BiofileGroup Exceptions
# -----------------------------------------------------------------------------


class BiofileGroupValidationError(Exception):
    """Basic exception for errors raised by `BiofileGroup` validation checks

    Parameters
    ----------
    name
        The filenames which caused the error

    Attributes
    ----------
    msg
        The formatted error message
    """
    def __init__(self, names: Sequence[PurePath] = None) -> None:
        self.names = names
        super().__init__()

    def _formatted_filenames(self) -> str:
        """Construct the part of the error message which lists the file"""
        return '\n'.join(['Files:', '\t' + ', '.join(str(self.names))])


class FileExtensionsNotSameError(BiofileGroupValidationError):
    """Exception raised when files do not have the expected file extension"""

    def _err_description(self) -> str:
        return "File extensions are not all equal"

# -----------------------------------------------------------------------------
# MatchedPrefixGroup Exceptions
# -----------------------------------------------------------------------------


class MatchedPrefixGroupValidationError(BiofileGroupValidationError):
    """Basic exception for errors raised by `MatchedPrefixGroup` validation

    Parameters
    ----------
    groups
        Nested list of filenames involved in the error.
    """
    def __init__(self, groups: Sequence[BiofileGroup] = None) -> None:
        self.groups = groups
        super().__init__()

    def _formatted_filenames(self) -> str:
        """Format the file group names for printing"""
        formatted_filenames = 'Groups:\n'

        for group in self.groups:
            filename_str = ', '.join(as_strs(group._paths))
            formatted_filenames += ''.join(['[', filename_str, ']\n'])
        return formatted_filenames


class DuplicateFilegroupError(MatchedPrefixGroupValidationError):
    """Exception raised when a `BiofileGroup` is matched with itself"""
    pass


class PrefixMatchError(MatchedPrefixGroupValidationError):
    """Exception raised the files do not share a prefix"""
    pass


class GroupLengthError(MatchedPrefixGroupValidationError):
    """Exception raised `MatchedPrefixGroups`s don't have equal lengths"""
    pass
