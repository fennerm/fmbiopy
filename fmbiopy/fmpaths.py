""" Path manipulation utilities """

from pathlib import (
        Path,
        PurePath,
        )
from typing import (
        Dict,
        Iterable,
        List,
        Sequence,
        )


def absglob(directory: Path, pattern: str)-> List[Path]:
    """Like normal glob except absolute paths are returned"""
    relative = sorted(list(directory.glob(pattern)))
    absolute = [path.resolve() for path in relative]
    return absolute


def add_suffix(path: PurePath, suffix: str)-> Path:
    """Append a suffix to a path"""
    return Path(str(path) + suffix)


def all_exist(paths: Iterable[Path]) -> bool:
    """Return True if all paths in list exist """
    return all(apply_exists(paths))


def any_dont_exist(paths: Iterable[Path]) -> bool:
    """Return True if any path in list does not exist """
    return not all(apply_exists(paths))


def any_exist(paths: Iterable[Path]) -> bool:
    """Return True if any path in list exists """
    return any(apply_exists(paths))


def all_empty(paths: Iterable[Path]) -> bool:
    """Return True if all paths in list exist, and are non-empty"""
    check_all_exist(paths)
    return all(apply_is_empty(paths))


def apply_exists(paths: Iterable[Path]) -> List[bool]:
    """Apply os.path.exists to a list of paths"""
    return [p.exists() for p in paths]


def apply_is_empty(paths: Iterable[Path]) -> List[bool]:
    """Check whether each file in list has a nonzero file size.

    Returns
    -------
    A list of bools"""
    return [is_empty(p) for p in paths]


def as_dict(directory: Path) -> Dict[str, List[Path]]:
    """Represent the contents of a directory as a dictionary

    Parameters
    ----------
    directory
        A directory

    Returns
    -------
    Output dictionary is of the form Dict[x, y] where x is a subdirectory of
    `direc` and and y is the contents of x as a list.
    """
    subdirs = listdirs(directory)
    output = {}
    for d in subdirs:
        output[d.name] = absglob(d, '*')
    return output


def as_paths(
        strs: Iterable[str],
        absolute: bool = True,
        )-> List[Path]:
    """Convert a list of `str` to `pathlib.Path`

    Parameters
    ----------
    absolute
        If True, an absolute path is returned
    strict
        Passed to `pathlib.Path.resolve()`. Does nothing if `absolute` is
        False.
    """
    paths = [Path(string) for string in strs]
    if absolute:
        paths = [path.absolute() for path in paths]
    return paths


def as_strs(paths: Iterable[PurePath])-> List[str]:
    """Convert a list of `pathlib.Path` to string"""
    return [str(path) for path in paths]


def check_all_exist(paths: Iterable[Path]) -> None:
    """Raise OSError if any paths in list do not exist """
    if not all_exist(paths):
        raise OSError("Not all paths exist: \n" + ' '.join(as_strs(paths)))


def find(
        directory: Path,
        extensions: Iterable[str] = None,
        substring: str = None)-> List[Path]:
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
    List of files with the given extensions and substrings"""
    # Filter by file extension
    hits: List[Path] = []
    if extensions is not None:
        for ext in extensions:
            hits += absglob(directory, '*' + ext)
    else:
        hits = absglob(directory, '*')

    # Filter by substring
    if substring is not None:
        out: List[Path] = []
        for hit in hits:
            if substring in hit.name:
                out.append(hit)
    else:
        out = hits
    return sorted(out)


def get_bowtie2_indices(prefix: str)-> List[Path]:
    """Given the bowtie2 index prefix, return the bowtie2 indices"""
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                       '.rev.2.bt2']
    paths = [Path(prefix + suffix) for suffix in bowtie_suffixes]
    return paths


def is_empty(path: Path)-> bool:
    """Return True if `path` is empty"""
    return path.stat().st_size < 1


def listdirs(directory: Path) -> List[Path]:
    """List all the subdirectories of a directory"""
    contents = absglob(directory, '*')
    out = []
    for path in contents:
        if path.is_dir():
            out.append(path)
    return out


def prefix(path: PurePath)-> str:
    """Get the part of a string before the first dot"""
    return path.name.split('.')[0]


def resolve(path: str, *args, **kwargs)-> str:
    """Return the absolute path of a string path

    Additional parameters are passed to `Path.resolve()`
    """
    return str(Path(path).resolve(*args, **kwargs))

def root(path: Path)-> Path:
    """Return the root name of a path (with directory included)"""
    while '.' in path.name:
        path = path.parent / path.stem
    return path
