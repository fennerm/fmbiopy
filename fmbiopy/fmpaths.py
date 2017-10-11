""" Path manipulation utilities """

from glob import glob
import os
from re import sub
from typing import List
from typing import Sequence

import fmbiopy.fmcheck as fmcheck


def get_prefix(path: str) -> str:
    """Get the part of a string before the first dot"""
    return path.split(".")[0]


def get_suffix(path: str) -> str:
    """Get the part of a string after the first dot"""
    return ''.join(path.split(".")[1:])


def get_final_suffix(path: str) -> str:
    """Get the part of a string after the last dot"""
    dot_split = path.split(".")
    final_index = len(dot_split) - 1
    return dot_split[final_index]


def last_two_suffix(path: str) -> str:
    """Get the last two suffixes of a filename as a string"""
    dot_split = path.split(".")
    suffix_list = [dot_split[len(dot_split) - 2],
                   dot_split[len(dot_split) - 1]]
    return '.'.join(suffix_list)


def replace_suffix(path: str, old_suffix: str, new_suffix: str) -> str:
    """Replace the suffix of a string"""
    # Check that the given path has the given suffix
    fmcheck.check_suffix(path, old_suffix)
    # Replace the suffix
    return sub(old_suffix, new_suffix, path)


def add_suffix(name: str, suffix: str) -> str:
    """Append a suffix to a string"""
    return name + suffix


def remove_suffix(name: str, nremove=1) -> str:
    """Remove n suffixes from a string"""
    if nremove == 0:
        return name

    unsuffixed = os.path.splitext(name)[0]
    nremove -= 1
    return remove_suffix(unsuffixed, nremove)


def listdirs(directory: str) -> Sequence[str]:
    """List all the subdirectories of a directory"""
    # Get all paths in directory including regular files
    contents = glob(directory + '/*')

    # Convert to absolute paths
    abs_paths = [os.path.abspath(item) for item in contents]

    # Select only the directories
    dirs = []
    for path in abs_paths:
        if os.path.isdir(path):
            dirs.append(path)
    return dirs


def get_bowtie2_indices(prefix: str) -> List[str]:
    """Given the bowtie2 index prefix, return the bowtie2 indices"""
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                       '.rev.2.bt2']
    return sorted([prefix + suf for suf in bowtie_suffixes])


def match_files(
        directory : str,
        types : Sequence[str] = None,
        substring : str = None) -> List[str]:
    """List all files with a given file extension and/or substring

    Parameters
    ----------
    directory
        The name of the directory to search in
    types, optional
        List of target file extensions (e.g py, txt)
    substring, optional
        Target substring

    Returns
    -------
    List[str]
        List of files with the given extensions and substrings"""

    hits : List[str] = []

    # Filter by file extension
    if types:
        for typ in types:
            hits += glob(directory + '/*.' + typ)
    else:
        hits = glob(directory + '/*')

    # Filter by substring

    if substring:
        out = []
        for hit in hits:
            if substring in hit:
                out.append(hit)
    else:
        out = hits
    return sorted(out)
