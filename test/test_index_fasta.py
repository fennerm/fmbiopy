## Pytests for run_bowtie2.py

import os
import pytest
from get_dat import get_dat
from fmbiopy.index_fasta import index_fasta
from fmbiopy.fmsystem import delete
from fmbiopy.fmpaths import remove_suffix, add_suffix

def test_index_fasta():
    testdat = get_dat()
    assemblies = testdat['assemblies']
    root_indices = remove_suffix(assemblies)
    expected_indices = sorted(
            add_suffix(assemblies, '.fai') +
            add_suffix(root_indices, '.1.bt2') +
            add_suffix(root_indices, '.2.bt2') +
            add_suffix(root_indices, '.3.bt2') +
            add_suffix(root_indices, '.4.bt2') +
            add_suffix(root_indices, '.rev.1.bt2') +
            add_suffix(root_indices, '.rev.2.bt2'))

    # Test if it runs without error.
    with delete(expected_indices):
        indices = index_fasta(assemblies)
        indices = sorted(indices[0] + indices[1])

        # Test if indices actually produced
        assert indices == expected_indices
