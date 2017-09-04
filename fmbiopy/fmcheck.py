"""
Functions for checking that various conditions are met, or. Functions which
start with 'check_' raise exceptions, all others return Bools
"""

import os

def check_non_empty(items):
    """ Raise exception if list is empty """

    if not items:
        raise ValueError("List is empty")

def all_equal(items):
    """ Test whether all items in list are equal """

    return all(item == items[0] for item in items)

def any_dont_exist(paths):
    """ Return True if any path in list does not exist """

    check_non_empty(paths)
    exists = map(os.path.exists, paths)
    return any(not x for x in exists)


def check_all_exist(paths):
    """ Raise an exception if any paths in list do not exist """

    for path in paths:
        if not os.path.exists(path):
            raise OSError("Path doesn't exist: " + path)

def check_suffix(name, suffixes):
    """ Check if a string x ends with any of a list of suffixes """

    # A single suffix string is also allowed.
    if isinstance(suffixes, str):
        suffixes = [suffixes]

    correct = any(name.endswith(suffix) for suffix in suffixes)

    if not correct:
        raise ValueError(name + " does not have the correct suffix " + \
                         ' '.join(suffixes))

def check_all_suffix(names, suffixes):
    """ Check if all strings in a list end with one of a list of suffixes.
        If not, raise an exception """

    for name in names:
        check_suffix(name, suffixes)
