"""Utilities for system manipulation

Moving/creating files, running commands etc.
"""

import collections
from contextlib import contextmanager
import errno
import logging
import os
import subprocess
from typing import Any
from typing import Dict
from typing import Generator
from typing import List
from typing import Sequence
from typing import Tuple

import fmbiopy.fmtype as fmtype


class IncorrectCommandFormatError(Exception):
    """Raised when a command argument cannot be parsed"""


def run_command(
        command: Sequence,
        logger_id: str = '',
        log_stdout: bool = True,
        log_stderr: bool = True,
        # We don't specify type explicitely here due to circular import
        mutex_log=None,
        shell: bool = False) -> Tuple[int, str, str]:
    """Run a bash command with logging support

    Parameters
    ----------
    command
        Bash command to be run
    logger_id
        Name to use for logging handler (ignored if mutex_log given).
        By default, root logger is used.
    log_stdout, log_stderr
        Should standard out and standard error be logged?
    mutex_log
        If running in parallel, a RuffusLog instance can be passed to log
        output with a mutex lock.
    shell
        If True, the command is run directly in the shell rather than the
        python interpreter. Useful for Bash commands with piping.

    Returns
    -------
    A triple of the form (return code, standard out, standard error)
    """

    if shell:
        # If run in shell, command needs to be a string, not a list
        if not isinstance(command, str):
            command = list(filter(None, command))
            command = ' '.join(command)

        process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,  # UTF-8 encoding specification
                shell=shell)
    else:
        # If not run in shell, command needs to be a list
        if isinstance(command, str):
            command = command.split()

        # Remove empty list items
        command = list(filter(None, command))

        # Run the command
        process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,  # UTF-8 encoding specification
                shell=shell)

    stdout, stderr = process.communicate()

    # Log results
    if stdout or stderr:
        if mutex_log:
            if log_stdout and stdout:
                mutex_log.write(stdout)
            if log_stderr and stderr:
                mutex_log.write(stderr)
        else:
            logger = logging.getLogger(logger_id)
            if log_stdout and stdout:
                logger.info(stdout)
            if log_stderr and stderr:
                logger.info(stderr)

    return (int(process.returncode), stdout, stderr)


@contextmanager
def working_directory(directory: str) -> Generator:
    """Change working directory context safely.

    Usage
    -----
        with working_directory(directory):
            <code>
    """

    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


def remove_all(names: fmtype.StringOrSequence, silent: bool = False)-> None:
    """Remove all files given as either a string or list"""
    if silent:
        remove_func = silent_remove
    else:
        remove_func = os.remove  # type: ignore

    if isinstance(names, collections.Iterable) and not isinstance(names, str):
        for name in names:
            remove_func(name)
    else:
        remove_func(names)


@contextmanager
def delete(paths: Sequence[str]) -> Generator:
    """Context manager for deletion of temporary files.

    Context used for making sure that files are deleted even if an attempted
    action raises an exception. Useful for cleaning up temporary files.

    Usage
    -----
        with delete(paths):
            <code>
    """

    try:
        yield
    except Exception:
        remove_all(paths)
        raise
    finally:
        remove_all(paths)


def run_silently(command: Sequence[str]) -> Tuple[int, str, str]:
    """Run a command without logging results """
    return run_command(command, log_stdout=False, log_stderr=False)


def concat(filenames: Sequence[str], outpath: str) -> None:
    """Concatenate a list of files """
    filenames = ' '.join(filenames)
    command = 'cat ' + filenames + ' > ' + outpath
    err_code = run_silently(command)[0]
    if err_code != 0:
        raise OSError("File concatenation failed")


def mkdir(path: str) -> str:
    """Create a directory if it doesn't exist

    Returns
    -------
    Absolute path of the created directory
    """

    path = str(path)
    if not os.path.exists(path):
        os.makedirs(path)

    return os.path.abspath(path)


def mkdirs(dirnames: Sequence[str], output_directory: str) -> Sequence[str]:
    """Create a list directories

    Parameters
    ----------
    dirnames
        Names of directories to create
    output_directory
        Name of directory in which to create the directories

    Returns
    -------
        The paths of the created directories
    """

    with working_directory(str(output_directory)):
        abspaths = [mkdir(dirname) for dirname in dirnames]

    return abspaths


def silent_remove(filename : Any) -> None:
    """Try to remove a file, ignore exception if doesn't exist """
    try:
        os.remove(filename)
    except OSError as err:
        if err.errno != errno.ENOENT:
            raise


def dict_to_list(dictionary: Dict[str, str]) -> List[str]:
    """Convert a dictionary to a flat list """
    return [e for l in dictionary for e in l]


def parse_param_dict(param: Dict[str, str]) -> str:
    """Convert a parameter dictionary to a string BASH commands

    Parameters
    ----------
    param:
        A dictionary with argument flags (-x, --long etc.) as the keys and
        BASH parameter values as the values

    Returns
    -------
        A Bash command substring containing the parameters
    """
    if param:
        bash_string = ''
        for key, value in param.items():
            bash_string += (' '.join([key, value]))

        return bash_string
    return ''
