## Pytests for run_bowtie2.py

import os
import pytest
from get_dat import get_dat
from fmbiopy.run_index_fasta import run_index_fasta
import fmbiopy.fen_util as fen_util

def test_run_index_fasta():
    testdat = get_dat()
    assemblies = testdat['assemblies']
    # Test if it runs without error.
    indices = run_index_fasta(assemblies)
    indices = sorted(indices[0] + indices[1])

    root_indices = fen_util.remove_suffix(assemblies)
    expected_indices = sorted(fen_util.add_suffix(assemblies, '.fai') +
            fen_util.add_suffix(root_indices, '.1.bt2') +
            fen_util.add_suffix(root_indices, '.2.bt2') +
            fen_util.add_suffix(root_indices, '.3.bt2') +
            fen_util.add_suffix(root_indices, '.4.bt2') +
            fen_util.add_suffix(root_indices, '.rev.1.bt2') +
            fen_util.add_suffix(root_indices, '.rev.2.bt2'))

    # Test if indices actually produced
    assert indices == expected_indices

    # Clean up
    map(os.remove, indices)
