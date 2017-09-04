## Pytests for run_bowtie2.py

import os
import pytest
from get_dat import get_dat
from fmbiopy.run_index_fasta import run_index_fasta
from fmbiopy import fmpaths

def test_run_index_fasta():
    testdat = get_dat()
    assemblies = testdat['assemblies']

    # Test if it runs without error.
    try:
        indices = run_index_fasta(assemblies)
        indices = sorted(indices[0] + indices[1])

        root_indices = fmpaths.remove_suffix(assemblies)
        expected_indices = sorted(
                fmpaths.add_suffix(assemblies, '.fai') +
                fmpaths.add_suffix(root_indices, '.1.bt2') +
                fmpaths.add_suffix(root_indices, '.2.bt2') +
                fmpaths.add_suffix(root_indices, '.3.bt2') +
                fmpaths.add_suffix(root_indices, '.4.bt2') +
                fmpaths.add_suffix(root_indices, '.rev.1.bt2') +
                fmpaths.add_suffix(root_indices, '.rev.2.bt2'))

        # Test if indices actually produced
        assert indices == expected_indices
    finally:
        # Clean up
        map(os.remove, indices)

