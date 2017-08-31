#!/usr/bin/env/ python
import os, shutil
from multiprocessing import cpu_count
from itertools import izip
import fen_util

## Python wrapper around bowtie2
## Param:
##   fwd_read, rev_read = List of paths to fastq files (.fastq/.fastq.gz)
##   indices = List of paths to the reference fasta indices
##   out_paths = List of output file paths (Default: input[-.fastq].align.sam)
##   All other arguments refer directly to bowtie2 arguments
## Return:
##   The exit codes for each bowtie2 run as a list
def run_bowtie2(fwd_read, rev_read, indices, out_paths, threads=None, 
        maxins=None, no_discordant=None, no_mixed=None, preset=None):

    # Load default arguments if not specified

    if threads is None:
        threads = str(cpu_count() - 2)

    if maxins is None:
        maxins = '1000'

    if no_discordant is None:
        no_discordant = ''

    if no_mixed is None:
        no_mixed = ''

    if preset is None:
        preset = '--sensitive'

    # Construct the optional portion of the bash command
    opt_command = ['-p', str(threads), '-X', str(maxins), no_discordant, 
            no_mixed, preset]

    # Paths for unfinished output
    temp_out = fen_util.add_suffix(out_paths, '.part')

    # Run bowtie
    exit_codes = []
    for r1, r2, idx, tmp, out in izip(fwd_read, rev_read, indices, temp_out, 
            out_paths):

        # Bowtie2 sometimes segfaults with gzipped fastq files, so we gunzip
        # them
        if fen_util.final_suffix(r1) == '.gz':
            r1 = '<(gunzip -c ' + r1 + ')'
        if fen_util.final_suffix(r2) == '.gz':
            r2 = '<(gunzip -c ' + r2 + ')'
            
        # Construct command
        command = (['bowtie2', '-1', r1, '-2', r2, '-x', idx, '-S', tmp] + 
                opt_command)

        # Run command
        exit_code = fen_util.run_command(command, ("bowtie2." +
            fen_util.remove_suffix(os.path.basename(idx))))[0]
        exit_codes.append(exit_code)

        # Remove the temp file
        shutil.move(tmp, out)

    return exit_codes
