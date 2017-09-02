#!/usr/bin/env python

import os
from glob import glob
def get_dat():
    dat = {
        'fwd_reads' : sorted([os.path.abspath(x) for x in
                              glob('testdat/reads/fwd/*')]),
        'rev_reads' : sorted([os.path.abspath(x) for x in
                              glob('testdat/reads/rev/*')]),
        'assemblies' : sorted([os.path.abspath(x) for x in
                               glob('testdat/assembly/*fa*')])
    }
    return dat
