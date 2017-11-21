"""Utilities for system manipulation

Moving/creating files, running commands etc.
"""

from collections import Iterable as Iterable_
from contextlib import contextmanager
from errno import ENOENT
from logging import getLogger
from os import chdir
from pathlib import Path
from subprocess import (
        Popen,
        PIPE,
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
from warnings import warn

from fmbiopy.fmlist import (
        as_strs,
        exclude_blank,
        )


def bash(
        command: Sequence,
        logger_id: str = '',
        log: Tuple[bool, bool] = (False, True),
        shell: bool = True)-> Tuple[int, str, str]:
    """Run a bash command with logging support

    Command can either be a list or a string

    Parameters
    ----------
    command
        Bash command to be run
    logger_id
        Name to use for logging handler. By default, root logger is used.
    log
        If log[0] is True, stdout will be logged. If log[1] is True, stderr
        will be logged. By default only stderr is logged.
    shell
        If True, the command is run directly in the shell rather than the
        python interpreter. Useful for Bash commands with piping.

    Returns
    -------
    Tuple[int, str, str]
        A tuple of the form (return code, standard out, standard error)
    """
    if not isinstance(command, str):
        command = exclude_blank(command)
        command = as_strs(command)

    if shell:
        # If run in shell, command needs to be a string, not a list
        if not isinstance(command, str):
            command = ' '.join(command)
    else:
        if isinstance(command, str):
            command = command.split(' ')

    process = Popen(
            command,
            stdout=PIPE,
            stderr=PIPE,
            universal_newlines=True,  # UTF-8 encoding specification
            shell=shell)

    output = process.communicate()

    # Write log
    for msg, stream in zip(output, log):
        if msg and stream:
            getLogger(logger_id).info(msg)

    return (int(process.returncode), output[0], output[1])


@contextmanager
def working_directory(directory: Path) -> Generator[Path, None, None]:
    """Change working directory context safely.

    Usage
    -----
        with working_directory(directory):
            <code>
    """

    owd = Path.cwd()
    try:
        chdir(str(directory))
        yield directory
    finally:
        chdir(str(owd))


def remove_all(names: Iterable[Path], silent: bool = False)-> None:
    """Remove all files given as either a string or list"""
    if silent:
        remove_func: Callable[[Path], Any] = silent_remove
    else:
        remove_func = Path.unlink

    if isinstance(names, Iterable_):
        for name in names:
            if name:
                try:
                    remove_func(name)  # type: ignore
                except IsADirectoryError:
                    if silent:
                        warn('Attempted to delete a directory. Skipping')
                    else:
                        raise
    else:
        remove_func(names)  # type: ignore


@contextmanager
def delete(paths: Iterable[Path]) -> Generator[None, None, None]:
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

def run_command(*args, **kwargs):
    """Deprecated. Use bash"""
    return bash(*args, **kwargs)

def run_silently(command: Sequence[str])-> Tuple[int, str, str]:
    """Run a command without logging results """
    return bash(command, log=(False, False))


def concat(filenames: Iterable[Path], outpath: Path)-> bool:
    """Concatenate a list of files """
    command = ' '.join(['cat', *as_strs(filenames), '>', outpath.name])
    err_code = run_silently(command)[0]
    if err_code != 0:
        raise OSError("File concatenation failed")
    return True


def silent_remove(filename: Path) -> None:
    """Try to remove a file, ignore exception if doesn't exist """
    if filename is not None:
        try:
            filename.unlink()
        except OSError as err:
            if err.errno != ENOENT:
                raise


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
