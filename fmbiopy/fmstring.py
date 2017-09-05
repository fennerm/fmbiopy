""" Utility functions for carrying out string manipulation """

from typing import Sequence, List

def get_unique(items: Sequence[str]) -> List[str]:
    """ Get the unique elements of a list """
    return sorted(list(set(items)))
