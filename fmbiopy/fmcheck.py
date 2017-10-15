"""Functions for checking that various conditions are met.

Functions which start with 'check_' raise exceptions, all others return bools
"""

from typing import Sequence


def check_non_empty(items: Sequence) -> None:
    """Raises ValueError if list is empty"""

    if not items:
        raise ValueError("List is empty")


def all_equal(items: Sequence) -> bool:
    """Test whether all items in list are equal """

    return all(item == items[0] for item in items)


def any_endswith(items: Sequence[str], suffix) -> bool:
    """Return True if any item in list ends with the given suffix """
    return any([item.endswith(suffix) for item in items])


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
