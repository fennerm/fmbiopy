import os
from glob import glob
from fmbiopy.fmpaths import abs_paths

def get_dat():
    dat = {
        'fwd_reads' : sorted(abs_paths(glob('testdat/reads/fwd/*.fq.gz'))),
        'rev_reads' : sorted(abs_paths(glob('testdat/reads/rev/*.fq.gz'))),
        'assemblies' : sorted(abs_paths(glob('testdat/assembly/*.fasta'))),
        'empty' : sorted(abs_paths(glob('testdat/reads/empty/*.fq')))
    }
    return dat
