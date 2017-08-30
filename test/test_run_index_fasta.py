## Pytests for run_bowtie2.py

import os
import pytest
from get_dat import get_dat
from fmbiopy.run_index_fasta import run_index_fasta

def test_run_index_fasta():
    testdat = get_dat()
    assemblies = testdat['assemblies']
    # Just test if it runs without error.
    indices = run_index_fasta(assemblies)
    map(os.remove, indices[0]+indices[1])
