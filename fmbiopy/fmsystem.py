"""Utilities for system manipulation

Moving/creating files, running commands etc.
"""

from collections import Iterable as Iterable_
from contextlib import contextmanager
from logging import getLogger
from os import chdir
from subprocess import (
        Popen,
        PIPE,
        STDOUT,
        )
from typing import (
        Any,
        Callable,
        Dict,
        Generator,
        Iterable,
        Sequence,
        Tuple,
        )

from plumbum import (
        local,
        LocalPath,
        )
from plumbum.cmd import cat

from fmbiopy.fmlist import (
        as_strs,
        exclude_blank,
        )


def concat(filenames: Iterable[LocalPath], outpath: LocalPath):
    """Concatenate a list of files """
    (cat.__getitem__(filenames) > outpath)()


def parse_param_dict(param: Dict[str, str]) -> str:
    """Convert a parameter dictionary to a string BASH commands

    Parameters
    ----------
    param
        A dictionary with argument flags (-x, --long etc.) as the keys and
        BASH parameter values as the values

    Returns
    -------
    str
        A Bash command substring containing the parameters
    """
    if param:
        bash_string = ''
        for key, value in param.items():
            bash_string += (' '.join([key, value]))

        return bash_string
    return ''
