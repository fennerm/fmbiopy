import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
from fmbiopy.fmsystem import delete
from get_dat import get_dat
from glob import glob
import os
import pytest


def test_samtools_index_fasta():
    assemblies = get_dat()['assemblies']
    root_indices = fmpaths.remove_suffix(assemblies)
    expected_indices = sorted(fmpaths.add_suffix(assemblies, '.fai'))

    # Test if it runs without error.
    with delete(expected_indices):
        samtools_index_fasta(assemblies, expected_indices)
        assembly_dir = os.path.dirname(assemblies[0])
        indices = sorted(glob(assembly_dir + "/*.fai"))

        # Test if indices actually produced
        assert indices == expected_indices


def test_bowtie_index_fasta():
    assemblies = get_dat()['assemblies']
    root_indices = fmpaths.remove_suffix(assemblies)
    expected_indices = sorted(
            fmpaths.add_suffix(root_indices, '.1.bt2') +
            fmpaths.add_suffix(root_indices, '.2.bt2') +
            fmpaths.add_suffix(root_indices, '.3.bt2') +
            fmpaths.add_suffix(root_indices, '.4.bt2') +
            fmpaths.add_suffix(root_indices, '.rev.1.bt2') +
            fmpaths.add_suffix(root_indices, '.rev.2.bt2'))

    # Test if it runs without error.
    with delete(expected_indices):
        bowtie_index_fasta(assemblies)
        assembly_dir = os.path.dirname(assemblies[0])
        indices = sorted(glob(assembly_dir + "/*.bt2"))

        # Test if indices actually produced
        assert indices == expected_indices

"""
## Pytests for bowtie2_align.py

import os
import pytest
from get_dat import get_dat
from fmbiopy.index_fasta import index_fasta
from fmbiopy.bowtie2_align import bowtie2_align
from fmbiopy.fmsystem import delete
from fmbiopy.fmpaths import remove_suffix

def test_bowtie2_align():
    testdat = get_dat()
    fwd_reads = testdat['fwd_reads']
    rev_reads = testdat['rev_reads']
    assemblies = testdat['assemblies']
    indices = index_fasta(assemblies)

    bowtie2_indices = remove_suffix(indices[0], 2)

    indices = indices[0] + indices[1]

    out_paths = []
    for a in assemblies:
        base = os.path.basename(a)
        out_path = os.path.join(os.getcwd(), base+'.sam')
        out_paths.append(out_path)

    with delete(out_paths + indices):
        exit_codes = bowtie2_align(fwd_reads, rev_reads, bowtie2_indices,
                out_paths, threads=2)
        if not all(code == 0 for code in exit_codes):
            raise RuntimeError("Not all bowtie2 exit codes non-zero")
"""
