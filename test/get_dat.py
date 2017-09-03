import os
from glob import glob

def get_dat():
    dat = {
        'fwd_reads' : sorted([os.path.abspath(x) for x in
                              glob('testdat/reads/fwd/*.fq.gz')]),
        'rev_reads' : sorted([os.path.abspath(x) for x in
                              glob('testdat/reads/rev/*.fq.gz')]),
        'assemblies' : sorted([os.path.abspath(x) for x in
                               glob('testdat/assembly/*.fasta')])
    }
    return dat
