from fmbiopy.fmpaths import listdirs
from glob import glob
import os


def get_dat():
    testdirs = [os.path.abspath(d) for d in listdirs('testdat')]
    dat = {}
    for d in testdirs:
        base = os.path.basename(d)
        dat[base] = sorted(glob(d + '/*'))
    return dat
