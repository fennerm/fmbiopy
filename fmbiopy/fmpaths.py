""" Path manipulation utilities """

import os
from re import sub
from typing import Sequence
from fmbiopy.fmcheck import check_suffix
from fmbiopy.fmtype import PathsOrStrings, StringOrSequence

def paths_to_string(paths: PathsOrStrings) -> str:
    """ Convert a list of Paths to a space separated string """
    return ' '.join(str(p) for p in paths)

def get_prefix(path: str) -> str:
    """ Get the part of a string before the first dot """
    return path.split(".")[0]

def get_suffix(path: str) -> str:
    """ Get the part of a string after the first dot """
    return ''.join(path.split(".")[1:])

def final_suffix(path: str) -> str:
    """ Get the part of a string after the last dot """
    dot_split = path.split(".")
    final_index = len(dot_split) - 1
    return dot_split[final_index]

def replace_suffix(path: str, old_suffix: str, new_suffix: str) -> str:
    """ Replace the suffix of a string """

    check_suffix(path, old_suffix)

    return sub(path, old_suffix, new_suffix)

def add_suffix(names: StringOrSequence, suffix: str) -> Sequence[str]:
    """ Append a suffix to a string or list of strings """

    if isinstance(names, str):
        return names + suffix

    return [add_suffix(name, suffix) for name in names]

def remove_suffix(names: StringOrSequence, nremove=1) -> StringOrSequence:
    """ Remove n suffixes from a list of strings """

    if nremove == 0:
        return names

    if isinstance(names, str):

        # If x is a string its simple
        unsuffixed = os.path.splitext(names)[0]
    else:
        # If x is a list, recursively apply
        unsuffixed = [remove_suffix(name, 1) for name in names]

    nremove -= 1

    return remove_suffix(unsuffixed, nremove)

def abs_paths(paths: Sequence[str]) -> Sequence[str]:
    """ Convert list of relative paths to absolute """
    return [os.path.abspath(path) for path in paths]
