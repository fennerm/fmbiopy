"""Functions for checking that various conditions are met.

Functions which start with 'check_' raise exceptions, all others return bools
"""

def all_equal(items):
    """Test whether all items in list are equal"""
    return all(item == items[0] for item in items)


def any_endswith(items, suffix):
    """Return True if any item in list ends with the given suffix

    Parameters
    ----------
    items: Sequence[str]
        List of strings
    suffix: str
        A suffix string

    Returns
    -------
    bool
        True if any element in `items` ends with `suffix`
    """
    return any([item.endswith(suffix) for item in items])


def check_suffix(x, suffixes):
    """Check that a string has one of a list of suffixes

    Parameters
    ----------
    x: str
        Possibly suffixed string
    suffixes: Sequence[str]
        List of suffixes to check

    Raises
    ------
    ValueError
        If x does not have one of `suffixes`
    """

    has_suffix = any([x.endswith(suffix) for suffix in suffixes])

    if not has_suffix:
        raise ValueError(
                ' '.join([
                    x, "does not have the correct suffix", ' '.join(suffixes)]))


def check_all_suffix(names, suffixes):
    """Raise ValueError if any name does not have any of a list of suffixes"""
    for name in names:
        check_suffix(name, suffixes)
