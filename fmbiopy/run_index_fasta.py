#!/usr/bin/env/ python
import fen_util
import os
from glob import glob

## Index a list of fasta files using samtools faidx and bowtie2-build 
## Param
##   references  List; Fasta files
## Return 
##   A tuple (list of created .fai indices, list of created bowtie2 indices)
def run_index_fasta(references):

    # Check arguments
    fen_util.check_all_exist(references)
    for ref in references:
        fen_util.check_file_extension(references, [".fasta", ".fa", ".fna", ".mfa"])

    # Convert to absolute paths
    references = map(os.path.abspath, references)

    command = ['bowtie2_index'] + references
    fen_util.run_command(command, "build_index", False) 

    faidx_indices = []
    bt2_indices = []

    for ref in references:
        faidx_indices.append(ref+'.fai')
        bt2_indices = bt2_indices + glob(ref+'*.bt2')

    return (sorted(faidx_indices), sorted(bt2_indices))
