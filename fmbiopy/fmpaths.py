""" Path manipulation utilities """

from glob import glob
import os
from re import sub
import typing

import fmbiopy.fmcheck as fmcheck
from fmbiopy.fmtype import StringOrSequence


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


def add_suffix(names: str, suffix: str) -> str:
    """Append a suffix to a string or list of strings"""
    return names + suffix


def remove_suffix(names: StringOrSequence, nremove=1) -> typing.List[str]:
    """Remove n suffixes from a list of strings """
    if nremove == 0:
        # To ensure that return type is consistent, we return a list if given a
        # string
        if isinstance(names, str):
            return [names]
        return list(names)

    if isinstance(names, str):
        # If x is a string its simple
        unsuffixed = [os.path.splitext(names)[0]]
    else:
        # If x is a list, recursively apply
        unsuffixed = [remove_suffix(name, 1)[0] for name in names]

    nremove -= 1
    return remove_suffix(unsuffixed, nremove)


def listdirs(directory: str) -> typing.Sequence[str]:
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


def get_bowtie2_indices(prefix: str) -> typing.List[str]:
    """Given the bowtie2 index prefix, return the bowtie2 indices"""
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
                       '.rev.2.bt2']
    return sorted([prefix + suf for suf in bowtie_suffixes])
