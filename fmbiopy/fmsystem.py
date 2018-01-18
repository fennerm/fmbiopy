"""Utilities for system manipulation

Moving/creating files, running commands etc.
"""
from plumbum.cmd import cat


def concat(filenames, outpath):
    """Concatenate a list of files """
    (cat.__getitem__(filenames) > outpath)()
