import os, sys
from pathlib import Path
import re

def paths_to_string(paths):
    """ Convert a list of Paths to a space separated string """
    s = ' '.join(str(p) for p in paths)
    return s

def prefix(x):
    """ Get the part of a string before the first dot """
    return x.split(".")[0]

def suffix(x):
    """ Get the part of a tring after the first dot """
    return ''.join(x.split(".")[1:])

def final_suffix(x):
    """ Get the part of a string after the last dot """
    sp = x.split(".")
    return sp[len(sp)-1]

def replace_suffix(x, old_suffix, new_suffix):
    """ Replace the suffix of a string """
    if not x.endswith(old_suffix):
        raise ValueError('Given suffix does not match the actual suffix (' + 
                x + ', ' + old_suffix)
    else:
        x.sub(old_suffix, new_suffix)

def add_suffix(x, suffix):
    """ Append a suffix to a string """

    if isinstance(x, basestring):
        return x + suffix
    else:
        return [add_suffix(item, suffix) for item in x]

def remove_suffix(x, n=1):
    """ Remove n suffixes from a list of strings """
    if isinstance(x, basestring):
        # If x is a string
        unsuffixed = os.path.splitext(x)[0]
    else:
        unsuffixed = [remove_suffix(item, 1) for item in x]

    n = n-1
    if n == 0:
        return unsuffixed
    else:
        return remove_suffix(unsuffixed, n)
