from fmbiopy.bowtie_index_fasta import bowtie_index_fasta
from fmbiopy.fmpaths import remove_suffix, add_suffix
from fmbiopy.fmsystem import delete
from get_dat import get_dat
from glob import glob
import os
import pytest

def test_bowtie_index_fasta():
    assemblies = get_dat()['assemblies']
    root_indices = remove_suffix(assemblies)
    expected_indices = sorted(
            add_suffix(root_indices, '.1.bt2') +
            add_suffix(root_indices, '.2.bt2') +
            add_suffix(root_indices, '.3.bt2') +
            add_suffix(root_indices, '.4.bt2') +
            add_suffix(root_indices, '.rev.1.bt2') +
            add_suffix(root_indices, '.rev.2.bt2'))

    # Test if it runs without error.
    with delete(expected_indices):
        bowtie_index_fasta(assemblies)
        assembly_dir = os.path.dirname(assemblies[0])
        indices = sorted(glob(assembly_dir + "/*.bt2"))

        # Test if indices actually produced
        assert indices == expected_indices
