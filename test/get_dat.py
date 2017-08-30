#!/usr/bin/env python

import os
from glob import glob
def get_dat():
    dat = {
        'fwd_reads' : sorted(map(os.path.abspath, glob('testdat/reads/fwd/*'))),
        'rev_reads' : sorted(map(os.path.abspath, glob('testdat/reads/rev/*'))),
        'assemblies' : sorted(map(os.path.abspath, 
                                  glob('testdat/assembly/*fa*')))
    }
    return dat
