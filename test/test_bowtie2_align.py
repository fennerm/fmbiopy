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
