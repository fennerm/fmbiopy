"""Utility functions for carrying out list manipulation"""

from collections import Iterable as Iterable_
from typing import (
        Any,
        Callable,
        Iterable,
        Iterator,
        List,
        Sequence,
        )

from fmbiopy.fmerr import EmptyListError
from fmbiopy.fmtype import T


def any_endswith(items: Sequence[str], suffix)-> bool:
    """Return True if any item in list ends with the given suffix """
    return any([item.endswith(suffix) for item in items])


def all_equal(items: Sequence)-> bool:
    """Test whether all items in list are equal """
    return all(item == items[0] for item in items)


def as_strs(x: Sequence[Any])-> List[str]:
    """Convert a list of items to their string representations"""
    return [str(xi) for xi in x]


def ensure_list(x: Any)-> List:
    """If not a sequence, convert to one"""
    if isinstance(x, str):
        return [x]

    try:
        x[0]
    except TypeError:
        return [x]

    return x


def exclude_blank(seq: Iterable[T])-> List[T]:
    """Remove empty and None items from an iterable"""
    return list(filter(bool, seq))


def flatten(sequence: Iterable)-> List[Any]:
    """Convert a sequence of sequences to a single flat sequence.

    Works on dictionaries, tuples, lists.
    """
    def _get_flat(sequence: Iterable) -> Iterator[Any]:
        for item in sequence:
            if (isinstance(item, Iterable_) and not
                    isinstance(item, (str, bytes))):
                yield from _get_flat(item)
            else:
                yield item

    return list(_get_flat(sequence))


def get_unique(items: Sequence[str]) -> List[str]:
    """ Get the unique elements of a string list """
    return sorted(list(set(items)))


def interleave(list1: Sequence[T], list2: Sequence[T]) -> List[T]:
    """Convert two lists to a single interleaved list"""
    if len(list1) != len(list2):
        raise ValueError('Lists are not the same length')
    return [val for pair in zip(list1, list2) for val in pair]


def none(x: Iterable[bool])-> bool:
    """Return True if all elements in `x` are False"""
    return all([not i for i in x])


def not_empty(func: Callable)-> Callable:
    """Function decorator for functions which require nonempty list input"""
    def _wrapper(*args, **kwargs):
        if not args[0]:
            raise EmptyListError
        return func(*args, **kwargs)
    return _wrapper


def split_list(x: Sequence[T], at: T) -> List[List[T]]:
    """Split a list into sublists by a specific item value

    E.g split_list(['a', 'b', '.', 'c', '.', 'd') = [['a', 'b'], ['c'], ['d']

    Parameters
    ----------
    x
        A list
    at
        The item value to split by

    Returns
    -------
    A list of sublists split by `at`
    """
    if at in x:
        # Stores the intermediate sublists for output
        output: List[List[T]] = []
        # Stores items in between split values
        buffer_: List[T] = []

        for item in x:
            if item == at:
                output.append(buffer_)
                buffer_ = []
            else:
                buffer_.append(item)

        output.append(buffer_)
        return output

    return [list(x)]
