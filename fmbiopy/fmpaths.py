""" Path manipulation utilities """

from contextlib import contextmanager
from errno import ENOENT
from os import chdir
from typing import (
        Dict,
        Generator,
        Iterable,
        List,
        )
from warnings import warn

from plumbum import (
        local,
        LocalPath,
        )

from fmbiopy.fmlist import not_empty


@not_empty
def all_exist(paths: Iterable[LocalPath])-> bool:
    """Return True if all paths in list exist """
    return all(apply_exists(paths))


@not_empty
def any_dont_exist(paths: Iterable[LocalPath]) -> bool:
    """Return True if any path in list does not exist """
    return not all(apply_exists(paths))


@not_empty
def any_exist(paths: Iterable[LocalPath]) -> bool:
    """Return True if any path in list exists """
    return any(apply_exists(paths))


@not_empty
def all_empty(paths: Iterable[LocalPath]) -> bool:
    """Return True if all paths in list exist, and are non-empty"""
    check_all_exist(paths)
    return all(apply_is_empty(paths))


def apply_is_empty(paths: Iterable[LocalPath]) -> List[bool]:
    """Check whether each file in list has a nonzero file size.

    Returns
    -------
    List[bool]
        True if empty for each path"""
    return [is_empty(p) for p in paths]


def apply_exists(paths: Iterable[LocalPath])-> List[bool]:
    """Apply os.path.exists to a list of paths"""
    return [p.exists() for p in paths]


def as_dict(directory: LocalPath) -> Dict[str, List[LocalPath]]:
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


def as_paths(strs: Iterable[str])-> List[LocalPath]:
    """Convert a list of `str` to `plumbum.LocalPath`

    Returns
    -------
    List[Path]
        The list as paths
    """
    return [local.path(string) for string in strs]


def as_strs(paths: Iterable[LocalPath])-> List[str]:
    """Convert a list of `plumbum.LocalPath` to string"""
    return [str(path) for path in paths]


def check_all_exist(paths: Iterable[LocalPath])-> None:
    """Raise OSError if any paths in list do not exist """
    if not all_exist(paths):
        raise FileNotFoundError(
                "Not all paths exist: \n" + ' '.join(as_strs(paths)))


def create_all(paths: Iterable[LocalPath])-> None:
    """Given a list of nonexistant paths, create them"""
    for path in paths:
        path.touch()


@contextmanager
def delete(paths: Iterable[LocalPath]) -> Generator[None, None, None]:
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


def extension_without_gz(path: LocalPath)-> str:
    """Get the file extension ignoring .gz if present"""
    suffixes = path.suffixes
    last = len(suffixes)
    if suffixes[last] == '.gz':
        extension = suffixes[last-1]
    else:
        extension = suffixes[last]
    return extension


def find(
        directory: LocalPath,
        extensions: Iterable[str] = None,
        substring: str = None)-> List[LocalPath]:
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
    hits: List[LocalPath] = []
    if extensions is not None:
        for ext in extensions:
            hits += directory // ('*' + ext)
    else:
        hits = directory // '*'

    # Filter by substring
    if substring is not None:
        out: List[LocalPath] = []
        for hit in hits:
            if substring in hit.name:
                out.append(hit)
    else:
        out = hits
    return sorted(out)


def get_bowtie2_indices(prefix: str)-> List[LocalPath]:
    """Given the bowtie2 index prefix, return the bowtie2 indices"""
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                       '.rev.2.bt2']
    paths = [LocalPath(prefix + suffix) for suffix in bowtie_suffixes]
    return paths


def is_empty(path: LocalPath)-> bool:
    """Return True if `path` is empty"""
    return size(path) < 1


def listdirs(directory: LocalPath)-> List[LocalPath]:
    """List all the subdirectories of a directory"""
    contents = directory // '*'
    out = []
    for path in contents:
        if path.is_dir():
            out.append(path)
    return out

def move(src: LocalPath, dest: LocalPath)-> None:
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


def prefix(path: LocalPath)-> str:
    """Get the part of a string before the first dot"""
    return root(path).name


def remove_all(names: Iterable[LocalPath], silent: bool = False)-> None:
    """Remove all files given as either a string or list"""
    if silent:
        remove_func = silent_remove
    else:
        remove_func = LocalPath.delete

    if isinstance(names, Iterable):
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


def rm_gz_suffix(path: LocalPath)-> LocalPath:
    """Gzip aware suffix removal

    If .gz is the final path extension then two extensions are removed."""
    if path.suffix == '.gz':
        nremove = 2
    else:
        nremove = 1
    return path.with_suffix('', nremove)


def root(path: LocalPath)-> LocalPath:
    """Return the root name of a path (with directory included)"""
    return path.with_suffix('', None)


def silent_remove(filename: LocalPath) -> None:
    """Try to remove a file, ignore exception if doesn't exist """
    if filename is not None:
        try:
            filename.delete()
        except OSError as err:
            if err.errno != ENOENT:
                raise


def size(path: LocalPath)-> int:
    """Return the filesize of a path in bytes"""
    return path.stat().st_size
