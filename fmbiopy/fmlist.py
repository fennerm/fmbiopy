"""Utility functions for carrying out list manipulation"""

try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence


from fmbiopy.fmerr import EmptyListError


def any_endswith(items, suffix):
    """Return True if any item in list ends with the given suffix """
    return any([item.endswith(suffix) for item in items])


def all_equal(items):
    """Test whether all items in list are equal """
    return all(item == items[0] for item in items)


def as_strs(x):
    """Convert a list of items to their string representations"""
    return [str(xi) for xi in x]


def ensure_list(x):
    """If not a list, convert to one"""
    if isinstance(x, str):
        return [x]

    try:
        x[0]
    except TypeError:
        return [x]

    return x


def exclude_blank(seq):
    """Remove empty and None items from an iterable"""
    return [x for x in seq if x]


def flatten(items):
    """Convert a sequence of sequences to a single flat sequence.

    Works on dictionaries, tuples, lists.
    """
    result = []
    for item in items:
        if isinstance(item, list):
            result += flatten(item)
        else:
            result.append(item)
    return result


def get_unique(items):
    """ Get the unique elements of a string list """
    unique = list(set(items))
    try:
        unique = sorted(unique)
    except TypeError:
        pass
    return unique


def interleave(list1, list2):
    """Convert two lists to a single interleaved list"""
    if len(list1) != len(list2):
        raise ValueError('Lists are not the same length')
    return [val for pair in zip(list1, list2) for val in pair]


def is_non_string_sequence(x):
    """Return True if x is a non-string sequence"""
    return isinstance(x, Sequence) and not isinstance(x, (str, bytes))


def none(x):
    """Return True if all elements in `x` are False"""
    return all([not i for i in x])


def not_empty(func):
    """Function decorator for functions which require nonempty list input"""

    def _wrapper(*args, **kwargs):
        if not args[0]:
            raise EmptyListError
        return func(*args, **kwargs)

    return _wrapper


def split_list(x, at):
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
        buffer_ = []

        for item in x:
            if item == at:
                output.append(buffer_)
                buffer_ = []
            else:
                buffer_.append(item)

        output.append(buffer_)
        return exclude_blank(output)

    return [x]
