import os, sys
from pathlib import Path

# Convert a list of paths to a space separated string
def paths_to_string(paths):
    s = ' '.join(str(p) for p in paths)
    return s

# Get the part of string before the first dot
def prefix(x):
    return x.split(".")[0]

# Get the part of string after the first dot
def suffix(x):
    return ''.join(x.split(".")[1:])

# Get the part of string after the last dot
def final_suffix(x):
    sp = x.split(".")
    return sp[len(sp)-1]

## Replace the suffix of a string
## Param-
##   x  String
##   old_suffix  Suffix to be replaced
##   new_suffix  Suffix to replace with
## Return-
##   The string with updated suffix
def replace_suffix(x, old_suffix, new_suffix):
    if not x.endswith(old_suffix):
        raise ValueError('Given suffix does not match the actual suffix (' + 
                x + ', ' + old_suffix)
    else:
        return x[:-len(old_suffix)] + new_suffix

## Add a suffix to a list of strings
def add_suffix(lst, suffix):
    suffixed = []
    for item in lst:
        suffixed.append(item + suffix)
    return suffixed

## Remove n suffixes from a list of string
def remove_suffix(x, n=1):
    if isinstance(x, basestring):
        # If x is a string
        unsuffixed = os.path.splitext(x)[0]
    else:
        # If x is a list
        unsuffixed = []
        for item in x:
            unsuffixed.append(os.path.splitext(item)[0])

    n = n-1
    if n == 0:
        return unsuffixed
    else:
        return remove_suffix(unsuffixed, n)
