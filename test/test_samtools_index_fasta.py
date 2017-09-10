import os
import pytest
from get_dat import get_dat
from fmbiopy.samtools_index_fasta import samtools_index_fasta
from fmbiopy.fmsystem import delete
from fmbiopy.fmpaths import remove_suffix, add_suffix
from glob import glob

def test_samtools_index_fasta():
    assemblies = get_dat()['assemblies']
    root_indices = remove_suffix(assemblies)
    expected_indices = sorted(add_suffix(assemblies, '.fai'))

    # Test if it runs without error.
    with delete(expected_indices):
        samtools_index_fasta(assemblies, expected_indices)
        assembly_dir = os.path.dirname(assemblies[0])
        indices = sorted(glob(assembly_dir + "/*.fai"))

        # Test if indices actually produced
        assert indices == expected_indices
