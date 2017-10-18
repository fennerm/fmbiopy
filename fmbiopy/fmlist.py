"""Utility functions for carrying out list manipulation"""

from collections import Iterable as Iterable_
from typing import Any
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Sequence

from fmbiopy.fmtype import T


def exclude_blank(sequence: Iterable[T])-> List[T]:
    """Remove empty and None items from an iterable"""
    return list(filter(bool, sequence))


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
