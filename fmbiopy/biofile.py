"""
Classes for storing and validating various types of bioinformatics files
"""

from fmbiopy.fmcheck import check_all_suffix, all_equal
from fmbiopy.fmpaths import final_suffix, remove_suffix

class BioFileGroup():
    """
    Basic class for storing and validating groups of bioinformatics files
    Classes for specific types of files inherit from it.
    """

    def __init__(self, files, gzipped = False):
        self.paths = files
        self.extensions = [final_suffix(x) for x in files]
        self._validated = False
        self.gzipped = gzipped
        self._prevalidate()

    def __getitem__(self, item):
        """ Get an file from the group using list syntax """

        return self.paths[item]

    def _check_paths_not_empty(self):
        """ Check paths not an empty list """
        if not self.paths:
            raise ValueError('Empty paths in BioFileGroup')

    def _check_paths_not_string(self):
        """ Check that paths is not a single string """

        if isinstance(self.paths, str):
            raise TypeError('Input to BioFileGroup is a string (expects \
                    list)')

    def _check_paths_same_extension(self):
        """ Check that paths all have the same extension """

        if self.gzipped:
            names = [remove_suffix(p, 1) for p in self.paths]
        else:
            names = self.paths

        extensions = [final_suffix(x) for x in names]

        if not all_equal(extensions):
            raise ValueError("Input paths to BioFileGroup do not have the \
                same extension")

    def _check_gzip(self):
        """ Check that the gzipped flag parameter matches the extensions of
        the files """

        self._check_paths_same_extension()
        if self.gzipped:
            if self.extensions[0] != 'gz':
                raise TypeError("Files are not gzipped")
        else:
            if self.extensions[0] == 'gz':
                raise TypeError("Files are gzipped, but gzip flag not set")

    def _prevalidate(self):
        """ Do some basic input validation upon initialization """

        self._check_paths_not_empty()
        self._check_paths_not_string()
        self._check_paths_same_extension()
        self._check_gzip()

    def __len__(self):
        """ Length of file list """
        return len(self.paths)

class FastqFileGroup(BioFileGroup):
    def __init__(self, files1, files2 = None):

        if files2 is None:
            super().__init__(paired)
            self.paired = False
        else:
            pass
    def _validate(self):
        pass

class BamFileGroup(BioFileGroup):
    def __init__(self, files, indices):
        pass

class SamFileGroup(BioFileGroup):
    def __init__():
        pass

class FastaFileGroup(BioFileGroup):
    def __init__():
        pass
