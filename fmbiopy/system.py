"""Utilities for system manipulation

Moving/creating files, running commands etc.
"""
from plumbum.cmd import cat

from fmbiopy.iter import exclude_blank


def concat(filenames, outpath):
    """Concatenate a list of files"""
    (cat.__getitem__(filenames) > outpath)()


def capture_stdout(command):
    '''Run a plumbum command and return the stdout as a list of strings'''
    stdout = command.run()[1]
    stdout = stdout.split('\n')
    stdout = exclude_blank(stdout)
    return stdout
