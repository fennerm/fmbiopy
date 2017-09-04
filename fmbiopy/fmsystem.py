"""
Utilities for system manipulation (moving/creating files, running commands etc.
"""

import os
import logging
import errno
from contextlib import contextmanager
from subprocess import Popen, PIPE

def run_command(command, logger_id=None, log_stdout=True, log_stderr=True):
    """
    Run a bash command with logging support

    Parameters
    ----------
    command (List)
        Bash command to be run
    logger_id (String)
        Name to use for logging handler
    log_stdout, log_stderr (Bool)
        Should standard out and standard error be logged?

    Returns
    -------
    A triple of the form (return code, standard out, standard error)

    """

    # If command is passed as a string, convert to list
    if isinstance(command, str):
        command = command.split()

    # Remove empty list items
    command = filter(None, command)

    # Run the command
    process = Popen(command, stdout=PIPE, stderr=PIPE,
                    universal_newlines=True)

    # UTF-8 encoding specification reqd for python 3
    stdout, stderr = process.communicate()

    # Log results
    logger = logging.getLogger(logger_id)
    if log_stdout and stdout:
        logger.info(stdout)
    if log_stderr and stderr:
        logger.info(stderr)

    return (process.returncode, stdout, stderr)

@contextmanager
def working_directory(directory):
    """
    Change working directory context safely.

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

def run_silently(command):
    """ Run a command without logging results """
    return run_command(command, log_stdout=False, log_stderr=False)

def concat(filenames, outpath):
    """ Concatenate a list of files """
    command = ' '.join(['cat'] + filenames + ['>', outpath])
    err_code = run_silently(command)[0]
    if err_code != 0:
        raise OSError("File concatenation failed")

# Create directory if it doesn't already exist
def mkdir(path):
    """
    Create a directory if it doesn't exist

    Returns
    -------
    Absolute path of the created directory
    """

    if not os.path.exists(path):
        os.makedirs(path)

    return os.path.abspath(path)




def mkdirs(dirnames, output_directory):
    """
    Create a list directories

    Parameters
    ----------
    dirnames - List
        Names of directories to create
    output_directory - String
        Name of directory in which to create the directories

    Returns
    -------
        The paths of the created directories
    """

    with working_directory(output_directory):
        abspaths = [mkdir(dirname) for dirname in dirnames]

    return abspaths


def silent_remove(filename):
    """ Try to remove a file, ignore exception if doesn't exist """
    try:
        os.remove(filename)
    except OSError as err:
        if err.errno != errno.ENOENT:
            raise
