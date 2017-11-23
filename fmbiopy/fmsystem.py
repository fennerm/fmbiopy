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


# def bash(
#         command: Sequence,
#         log: Tuple[bool, bool] = [True, True],
#         logfile: Path = None,
#         shell: bool = True,
#         console: bool = True)-> Tuple[int, str, str]:
#     """Run a bash command with logging support
#
#     Command can either be a list or a string
#
#     Parameters
#     ----------
#     command
#         Bash command to be run
#     logger_id
#         Name to use for logging handler. By default, no logger is used. Use ''
#         for the root logger.
#     logfile
#         Logfile to write output to. Note that this may be a different file than
#         that of the logger_id, depending on the logger setup.
#     log
#         If log[0] is true, stdout will be logged. If log[1] is true, stderr will
#         be logged. (Both can be true).
#     shell
#         If True, the command is run directly in the shell rather than the
#         python interpreter. Useful for Bash commands with piping.
#
#     Returns
#     -------
#     Tuple[int, str, str]
#         A tuple of the form (return code, standard out, standard error)
#     """
#     if not isinstance(command, str):
#         command = exclude_blank(command)
#         command = as_strs(command)
#
#     if shell:
#         # If run in shell, command needs to be a string, not a list
#         if not isinstance(command, str):
#             command = ' '.join(command)
#     else:
#         if isinstance(command, str):
#             command = command.split(' ')
#
#     process = Popen(
#             command,
#             stdout=PIPE,
#             stderr=PIPE,
#             universal_newlines=True,  # UTF-8 encoding specification
#             shell=shell)
#
#     output_stdout = []
#     output_stderr = []
#     # While receiving piped out and err
#     while process.poll() is not None:
#         out = process.stdout.readline()
#         err = process.stderr.readline()
#         for i, pipe, storage, stream in [
#                 (1, out, output_stdout, sys.stdout),
#                 (2, err, output_stderr, sys.stderr)]:
#             if log[0]:
#                 # Write to console
#                 if console:
#                     stream.write(pipe)
#                 # Write to logfile
#                 if logfile is not None:
#                     print(pipe, logfile.open('w'))
#                 # Store in variable for output
#                 storage += pipe
#
#     return (int(process.returncode), output_stdout, output_stderr)



# def run_command(*args, **kwargs):
#     """Deprecated. Use bash"""
#     return bash(*args, **kwargs)
#
# def run_silently(command: Sequence[str])-> Tuple[int, str, str]:
#     """Run a command without logging results """
#     return bash(command, log=(False, False))


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
