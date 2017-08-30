#!/usr/bin/env/ python
from multiprocessing import cpu_count
from itertools import izip
from subprocess import run
import fen_util
import os

## Python wrapper around bowtie2
## Param:
##   fwd_read, rev_read = List of paths to fastq files (.fastq/.fastq.gz)
##   indices = List of paths to the reference fasta indices
##   out_paths = List of output file paths (Default: input[-.fastq].align.sam)
##   All other arguments refer directly to bowtie2 arguments
def run_bowtie2(fwd_read, rev_read, indices, out_paths, threads=None, 
                maxins=None, no_discordant=None, no_mixed=None, preset=None):

    # Load default arguments if not specified
            
    if threads is None:
        threads = str(multiprocessing.cpu_count() - 2)
        
    if maxins is None:
        maxins = '1000'

    if no_discordant is None:
        no_discordant = ''

    if no_mixed is None:
        no_mixed = ''

    if preset is None:
        preset = 'sensitive'

    # Delete previous logfile
    os.remove(logfile)

    # Construct the optional portion of the bash command
    opt_command = ['-p', threads, '-X', maxins, no_discordant, no_mixed, 
                   preset]

    # Run bowtie
    for r1, r2, idx, out in itertools.izip(fwd_read, rev_read, indices,
                                           out_paths):

        command = (['bowtie2', '-1', r1, '-2', r2, '-x', idx, '-S', out] + 
                   opt_command)

        fen_util.run_command(command, "bowtie2" + basename(index))
