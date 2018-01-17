"""Parsing utilities"""
from collections import defaultdict
import sys

from docopt import (
        docopt,
        DocoptExit,
        )

def helpful_docopt(doc, *args, **kwargs):
    """Wrapper around docopt which prints the full usage message upon error

    Parameters
    ----------
    doc: str
        Docstring to parse

    Returns
    -------
    Dict[str, Union[str, List[str]]]
        Same output as normal docopt

    """
    try:
        opt = docopt(doc, help=False, *args, **kwargs)
        return opt
    except DocoptExit:
        print(doc)
        sys.exit()


def optget(args, short, long=None):
    """Parse command line arguments to fetch the value of an option

    Parameters
    ----------
    args: List[str]
        Command line arguments
    short: str
        The short name of the option (including dash)
    long: str, optional
        The long name of the option (including dashes)

    Returns
    -------
    str
        The option value.

    Raises
    ------
    ValueError
        If option is not present in args
    """
    try:
        # Try parse as a short option
        opt = args[args.index(short) + 1]
    except ValueError:
        if long:
            # Try parse as a long option
            for arg in args:
                if arg.startswith(long):
                    if arg == long:
                        opt = args[args.index(arg) + 1]
                    else:
                        opt = arg.split("=")[1]
            try:
                opt
            except NameError:
                raise ValueError("Opt not present in args")
        raise
    return opt
