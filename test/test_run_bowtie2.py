## Pytests for run_bowtie2.py

import os, pytest
from get_dat import get_dat
import fmbiopy.fen_util as fen_util
from fmbiopy.run_index_fasta import run_index_fasta
from fmbiopy.run_bowtie2 import run_bowtie2

def test_run_bowtie2():
    testdat = get_dat()
    fwd_reads = testdat['fwd_reads']
    rev_reads = testdat['rev_reads']
    assemblies = testdat['assemblies']
    indices = run_index_fasta(assemblies)
    
    bowtie2_indices = fen_util.remove_suffix(indices[0], 2)

    indices = indices[0] + indices[1]

    out_paths = []
    for a in assemblies:
        base = os.path.basename(a)
        out_path = os.path.join(os.getcwd(), base+'.sam')
        out_paths.append(out_path)

    try:
        exit_codes = run_bowtie2(fwd_reads, rev_reads, bowtie2_indices, 
                out_paths, threads=2)
        if not all(code == 0 for code in exit_codes):
            raise RuntimeError("Not all bowtie2 exit codes non-zero")

    finally:

        # Cleanup
        flag = 0
        for path in out_paths:
            if os.path.isfile(path):
                os.remove(path)
            else:
                flag = 1

        for idx in indices:
            if os.path.isfile(idx):
                os.remove(idx)

        assert flag == 0

