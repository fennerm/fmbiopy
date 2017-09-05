"""
Classes for storing and validating various types of bioinformatics files
"""

import os
from typing import Sequence, List
from fmbiopy.fmcheck import check_all_suffix, all_equal, any_endswith
from fmbiopy.fmpaths import final_suffix, last_two_suffix, remove_suffix

class BioFileGroup(object):
    """
    Basic class for storing and validating groups of bioinformatics files
    Classes for specific filetypes inherit from it.

    Attributes
    ----------
    paths
        Paths to files. The files need not exist upon initialization, but they
        are checked before the first access attempt.
    extensions
        The extensions of each file
    gzipped
        True if the files are gzipped
    possibly_empty
        Set this to True if it is acceptable for the files to be empty. This
        disables non-empty file checking.
    """

    def __init__(
            self,
            files: Sequence[str],
            gzipped: bool = False,
            possibly_empty: bool = False
            ) -> None:
        """ initialization """

        self._paths = files
        self._validated = False
        self.possibly_empty = possibly_empty
        self.gzipped = gzipped
        self._accepted_extensions = None
        self.extensions = self._get_extensions()

        self._prevalidate()

        # Sorting is done last because if it is done before prevalidation, the
        # sort will convert strings into lists of characters
        self._paths = sorted(files)

    def _get_extensions(self) -> Sequence[str]:
        """
        Get the file extensions of the files. If gzipped, return the two-part
        extension (E.g 'fq.gz')
        """
        if self.gzipped:
            return [last_two_suffix(f) for f in self._paths]
        else:
            return [final_suffix(f) for f in self._paths]

    def __getitem__(self, item) -> str:
        """ Get a file from the group using list syntax """
        if not self._validated:
            self._validate()

        return self._paths[item]

    def _prevalidate(self) -> None:
        """ Do some basic input validation upon initialization """

        self._check_paths_not_none()
        self._check_paths_not_string()
        self._check_gzip()

    def _validate(self) -> None:
        """ Validation function. Child classes implement their own """
        self._check_extensions()
        self._check_files_non_empty()
        self._validated = True

    def _check_files_non_empty(self) -> None:
        if not self.possibly_empty:
            for path in self._paths:
                if os.path.getsize(path) < 3:
                    raise ValueError(path + " is empty")

    def _check_paths_not_none(self) -> None:
        """ Check paths not an empty list """
        if not self._paths:
            raise ValueError('Empty paths in BioFileGroup')

    def _check_paths_not_string(self) -> None:
        """ Check that paths is not a single string """

        if isinstance(self._paths, str):
            raise TypeError('Input to BioFileGroup is a string (expects \
                    list)')

    def _check_extensions(self) -> None:
        """ Check that the file extensions match the accepted extensions """
        if self._accepted_extensions:
            for extension in self.extensions:
                if not extension in self._accepted_extensions:
                    raise ValueError("Invalid file extension")

    def _check_gzip(self) -> None:
        """ Check that the gzipped flag parameter matches the extensions of
        the files """

        if self.gzipped:
            if not any_endswith(self.extensions, 'gz'):
                raise TypeError("Files are not gzipped")
        else:
            if any_endswith(self.extensions, 'gz'):
                raise TypeError("Files are gzipped, but gzipped flag not set")

    def __len__(self) -> int:
        """ Length of file list """
        return len(self._paths)

class Fastq(BioFileGroup):
    """
    BioFileGroup for holding Fastq files. Reads can either be paired or
    unpaired.
    """
    def __init__(
            self,
            files: Sequence[str],
            gzipped: bool = False,
            possibly_empty: bool = False
            ) -> None:
        """ initialization """

        super().__init__(files, gzipped, possibly_empty)

        self._accepted_extensions = ['fastq', 'fq']
        if self.gzipped:
            self._accepted_extensions += ['fastq.gz', 'fq.gz']


class PairedFastq(object):
    """
    Stores two groups of paired Fastq files

    Attributes
    ----------
    fwd - Fastq, rev - Fastq
        Fastq objects to be paired

    """
    def __init__(
            self,
            forward_fastq: Fastq,
            reverse_fastq: Fastq
            ) -> None:
        """
        Initialization

        Parameters
        ----------
        forward_fastq - Fastq, reverse_fastq - Fastq
            Fastq objects to be paired
        """
        self.fwd = forward_fastq
        self.rev = reverse_fastq
        self._prevalidate()

    def _prevalidate(self) -> None:
        """
        Check that the PairedFastq is valid
        """
        self.__check_lengths()

    def __check_lengths(self) -> None:
        """
        Raise exception if lengths of the two Fastq objects are not equal
        """
        if len(self.fwd) != len(self.rev):
            raise ValueError("Fastq objects are not the same length")

    def __len__(self) -> int:
        """ Length of the PairedFastq """
        return len(self.fwd)

    def __getitem__(self, item) -> List[Fastq]:
        """ Get the pair of Fastq files indexed by 'item' """
        return list([self.fwd[item], self.rev[item]])

class Bam(BioFileGroup):
    def __init__(self, files, indices):
        pass

class Sam(BioFileGroup):
    def __init__():
        pass

class Fasta(BioFileGroup):
    def __init__():
        pass
