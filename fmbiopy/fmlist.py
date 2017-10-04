"""Utility functions for carrying out list manipulation"""

import collections
import typing


def get_unique(items: typing.Sequence[str]) -> typing.List[str]:
    """ Get the unique elements of a list """
    return sorted(list(set(items)))


def interleave(list1: typing.Sequence, list2: typing.Sequence) -> typing.List:
    """Convert two lists to a single interleaved list"""
    if len(list1) != len(list2):
        raise ValueError('Lists are not the same length')
    return [val for pair in zip(list1, list2) for val in pair]


def flatten(sequence: typing.Sequence) -> typing.List:
    """Convert a sequence of sequences to a single flat sequence.

    Works on dictionaries, tuples, lists.
    """
    def _get_flat(sequence: typing.Sequence) -> typing.List:
        for item in sequence:
            if (isinstance(item, collections.Iterable) and not
                    isinstance(item, (str, bytes))):
                yield from _get_flat(item)
            else:
                yield item

    return list(_get_flat(sequence))


def split_list(x: typing.Sequence, at: str) -> typing.List[typing.List]:
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
        output = []
        # Stores items in between split values
        buffer = []

        for item in x:
            if item == at:
                output.append(buffer)
                buffer = []
            else:
                buffer.append(item)

        output.append(buffer)
        return output

    return x
