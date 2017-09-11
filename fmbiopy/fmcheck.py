"""
Functions for checking that various conditions are met, or. Functions which
start with 'check_' raise exceptions, all others return bools
"""

from fmbiopy.fmtype import StringOrSequence
import os
import typing


def check_non_empty(items: typing.Sequence) -> None:
    """Raises ValueError if list is empty"""

    if not items:
        raise ValueError("List is empty")


def all_equal(items: typing.Sequence) -> bool:
    """Test whether all items in list are equal """

    return all(item == items[0] for item in items)


def exists(paths: typing.Sequence[str]) -> typing.List[bool]:
    """Apply os.path.exists to a list of paths"""
    return map(os.path.exists, paths)


def any_exist(paths: typing.Sequence[str]) -> bool:
    """Return True if any path in list exists """
    return any(exists(paths))


def any_dont_exist(paths: typing.Sequence[str]) -> bool:
    """Return True if any path in list does not exist """
    return not all(exists(paths))


def all_exist(paths: typing.Sequence[str]) -> bool:
    """Return True if all paths in list exist """
    return all(exists(paths))


def any_endswith(items: typing.Sequence[str], suffix) -> bool:
    """Return True if any item in list ends with the given suffix """
    return any([item.endswith(suffix) for item in items])


def check_all_exist(paths: typing.Sequence[str]) -> None:
    """Raise OSError if any paths in list do not exist """
    if not all_exist(paths):
        raise OSError("Not all paths exist: \n" + paths)


def check_suffix(name: str, suffixes: StringOrSequence) -> None:
    """Raise ValueError if name does not end with any of a list of suffixes """

    has_suffix = any([name.endswith(suffix) for suffix in suffixes])

    if not has_suffix:
        raise ValueError(name + " does not have the correct suffix " + \
                         ' '.join(suffixes))


def check_all_suffix(names: typing.Sequence[str],
                     suffixes: StringOrSequence) -> None:
    """Raise ValueError if any name does not have any of a list of suffixes"""

    for name in names:
        check_suffix(name, suffixes)
