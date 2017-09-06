#!/usr/bin/env/ python
import os
from glob import glob
from fmbiopy.fmcheck import check_all_exist, check_all_suffix
from fmbiopy.fmsystem  import run_command

## Index a list of fasta files using samtools faidx and bowtie2-build
## Param
##   references  List; Fasta files
## Return
##   A tuple (list of created .fai indices, list of created bowtie2 indices)
def index_fasta(references):

    # Check arguments
    check_all_exist(references)
    check_all_suffix(references, [".fasta", ".fa", ".fna", ".mfa"])

    # Convert to absolute paths
    references = [os.path.abspath(x) for x in references]

    command = ['bowtie2_index'] + references
    run_command(command, "build_index", False)

    faidx_indices = []
    bt2_indices = []
    for ref in references:
        faidx_indices.append(ref+'.fai')

        ref_no_suffix = os.path.splitext(ref)[0]

        bt2_indices = bt2_indices + glob(ref_no_suffix+'*.bt2')

    return (sorted(faidx_indices), sorted(bt2_indices))
