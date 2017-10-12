"""Functions for checking that various conditions are met.

Functions which start with 'check_' raise exceptions, all others return bools
"""

import os
from typing import List
from typing import Sequence


def check_non_empty(items: Sequence) -> None:
    """Raises ValueError if list is empty"""

    if not items:
        raise ValueError("List is empty")


def all_equal(items: Sequence) -> bool:
    """Test whether all items in list are equal """

    return all(item == items[0] for item in items)


def exists(paths: Sequence[str]) -> List[bool]:
    """Apply os.path.exists to a list of paths"""
    if isinstance(paths, str):
        return [os.path.exists(paths)]
    return [os.path.exists(p) for p in paths]


def any_exist(paths: Sequence[str]) -> bool:
    """Return True if any path in list exists """
    return any(exists(paths))


def any_dont_exist(paths: Sequence[str]) -> bool:
    """Return True if any path in list does not exist """
    return not all(exists(paths))


def all_exist(paths: Sequence[str]) -> bool:
    """Return True if all paths in list exist """
    return all(exists(paths))


def filesize_nonzero(paths: Sequence[str]) -> Sequence[bool]:
    """Check whether each file in list has a nonzero file size.

    Returns
    -------
    A list of bools"""
    if isinstance(paths, str):
        return [os.path.getsize(paths) > 0]
    return [os.path.getsize(p) > 0 for p in paths]


def all_filesize_nonzero(paths: Sequence[str]) -> bool:
    """Return True if all paths in list exist, and are non-empty"""
    check_all_exist(paths)
    return all(filesize_nonzero(paths))


def any_endswith(items: Sequence[str], suffix) -> bool:
    """Return True if any item in list ends with the given suffix """
    return any([item.endswith(suffix) for item in items])


def check_all_exist(paths: Sequence[str]) -> None:
    """Raise OSError if any paths in list do not exist """
    if not all_exist(paths):
        raise OSError("Not all paths exist: \n" + ' '.join(paths))


def check_suffix(name: str, suffixes: Sequence[str]) -> None:
    """Raise ValueError if name does not end with any of a list of suffixes """

    has_suffix = any([name.endswith(suffix) for suffix in suffixes])

    if not has_suffix:
        raise ValueError(name + " does not have the correct suffix " +
                         ' '.join(suffixes))


def check_all_suffix(names: Sequence[str], suffixes: Sequence[str]) -> None:
    """Raise ValueError if any name does not have any of a list of suffixes"""

    for name in names:
        check_suffix(name, suffixes)
