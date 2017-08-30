#!/usr/bin/env/ python
import fen_util
import os

def run_index_fasta(references):
    command = ['bowtie_index'] + references
    fen_util.run_command(command, "build_index") 
    indices = map(fen_util.final_suffix, references) 
    return indices
