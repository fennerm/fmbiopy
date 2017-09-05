import os
from glob import glob
from fmbiopy.fmpaths import abs_paths, listdirs

def get_dat():
    testdirs = abs_paths(listdirs('testdat'))
    dat = {}
    for d in testdirs:
        base = os.path.basename(d)
        dat[base] = sorted(glob(d + '/*'))
    return dat
