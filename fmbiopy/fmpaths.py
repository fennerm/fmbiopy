""" Path manipulation utilities """

from contextlib import contextmanager
from errno import ENOENT
from os import remove
from warnings import warn

from plumbum import local


from fmbiopy.fmlist import (
    is_non_string_sequence,
    not_empty,
)

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

@not_empty
def all_exist(paths):
    """Return True if all paths in list exist """
    return all(apply_exists(paths))


@not_empty
def any_dont_exist(paths):
    """Return True if any path in list doesnot exist """
    return not all(apply_exists(paths))


@not_empty
def any_exist(paths):
    """Return True if any path in list exists """
    return any(apply_exists(paths))


@not_empty
def all_empty(paths):
    """Return True if all paths in list exist, and are non-empty"""
    check_all_exist(paths)
    return all(apply_is_empty(paths))


def apply_is_empty(paths):
    """Check whether each file in list has a nonzero file size.

    Returns
    -------
    List[bool]
        True if empty for each path"""
    return [is_empty(p) for p in paths]


def apply_exists(paths):
    """Apply os.path.exists to a list of paths"""
    return [p.exists() for p in paths]


def as_dict(directory):
    """Represent the contents of a directory as a dictionary

    Parameters
    ----------
    directory
        A directory

    Returns
    -------
    Dict[str, List[LocalPath]]
        Output dictionary is of the form Dict[x, y] where x is a subdirectory
        of `direc` and and y is the contents of x as a list.
    """
    subdirs = listdirs(directory)
    output = {}
    for d in subdirs:
        output[d.name] = d // '*'
    return output


def as_paths(strs):
    """Convert a list of `str` to `plumbum.LocalPath`

    Returns
    -------
    List[Path]
        The list as paths
    """
    return [local.path(string) for string in strs]


def as_strs(paths):
    """Convert a list of `plumbum.LocalPath` to string"""
    return [str(path) for path in paths]



def check_all_exist(paths):
    """Raise OSError if any paths in list do not exist """
    if not all_exist(paths):
        raise FileNotFoundError(
                "Not all paths exist: \n" + ' '.join(as_strs(paths)))


def create_all(paths):
    """Given a list of nonexistant paths, create them"""
    for path in paths:
        path.touch()


@contextmanager
def delete(paths):
    """Context manager for deletion of temporary files.

    Context used for making sure that files are deleted even if an attempted
    action raises an exception. Useful for cleaning up temporary files.

    Usage
    -----
        with delete(paths):
            <code>
    """

    try:
        yield
    except Exception:
        remove_all(paths)
        raise
    finally:
        remove_all(paths)


def extension_without_gz(path):
    """Get the file extension ignoring .gz if present"""
    suffixes = path.suffixes
    last = len(suffixes)
    if suffixes[last] == '.gz':
        extension = suffixes[last-1]
    else:
        extension = suffixes[last]
    return extension


def find(directory, extensions=None, substring=None):
    """List all files with a given file extension and/or substring

    Parameters
    ----------
    directory
        The name of the directory to search in
    extensions, optional
        List of target file extensions (e.g py, txt)
    substring, optional
        Target substring

    Returns
    -------
    List[LocalPath]
        List of files with the given extensions and substrings"""
    # Filter by file extension
    hits = []
    if extensions is not None:
        for ext in extensions:
            hits += directory // ('*' + ext)
    else:
        hits = directory // '*'

    # Filter by substring
    if substring is not None:
        out = []
        for hit in hits:
            if substring in hit.name:
                out.append(hit)
    else:
        out = hits
    return sorted(out)


def get_bowtie2_indices(prefix):
    """Given the bowtie2 index prefix, return the bowtie2 indices"""
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                       '.rev.2.bt2']
    paths = [local.path(prefix + suffix) for suffix in bowtie_suffixes]
    return paths


def is_empty(path):
    """Return True if `path` is empty"""
    return size(path) < 1


def listdirs(directory):
    """List all the subdirectories of a directory"""
    contents = directory // '*'
    out = []
    for path in contents:
        if path.is_dir():
            out.append(path)
    return out

def move(src, dest):
    """Move a file

    Differs from `plumbum.local.path.move` in that when a symlink is
    encountered, the contents are copied to the new location and the link is
    destroyed"""
    if src.is_symlink():
        # copyfile(str(src), str(dest))
        src.copy(dest)
        src.unlink()
    else:
        src.move(dest)


def prefix(path):
    """Get the part of a string before the first dot"""
    return root(path).name


def remove_all(names, silent=False):
    """Remove all files given as either a string or list"""
    if silent:
        remove_func = silent_remove
    else:
        remove_func = remove

    if is_non_string_sequence(names):
        for name in names:
            if name:
                try:
                    remove_func(name)
                except IsADirectoryError:
                    if silent:
                        warn('Attempted to delete a directory. Skipping')
                    else:
                        raise
    else:
        remove_func(names)


def rm_gz_suffix(path):
    """Gzip aware suffix removal

    If .gz is the final path extension then two extensions are removed."""
    if path.suffix == '.gz':
        nremove = 2
    else:
        nremove = 1
    return path.with_suffix('', nremove)


def root(path):
    """Return the root name of a path (with directory included)"""
    return path.with_suffix('', None)


def silent_remove(filename):
    """Try to remove a file, ignore exception if doesn't exist """
    if filename is not None:
        try:
            filename.delete()
        except OSError as err:
            if err.errno != ENOENT:
                raise


def size(path):
    """Return the filesize of a path in bytes"""
    return path.stat().st_size
