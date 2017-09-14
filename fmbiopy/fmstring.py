"""Utility functions for carrying out string/list manipulation"""

import typing


def get_unique(items: typing.Sequence[str]) -> typing.List[str]:
    """ Get the unique elements of a list """
    return sorted(list(set(items)))


def interleave(list1: typing.Sequence, list2: typing.Sequence) -> typing.List:
    """Convert two lists to a single interleaved list"""
    if len(list1) != len(list2):
        raise ValueError('Lists are not the same length')
    return [val for pair in zip(list1, list2) for val in pair]
