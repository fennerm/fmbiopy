""" Various broadly useful Ruffus task definitions

Ruffus: http://www.ruffus.org.uk/
"""

import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
import logging
import os
from shutil import move
import typing


def bowtie_index_fasta(input_fasta: str,
        output_index: str,
        param: typing.Dict[str, str] = None,
        logger: fmruffus.RuffusLog = None,
        ) -> None:
    """Index a list of .fasta files using samtools faidx"""

    parsed_param = fmsystem.dict_to_list(param)

    # Construct the samtools command
    command = ['bowtie2-build'] + parsed_param + [input_fasta, output_index]

    # Run the command
    logger.info("Indexing fasta sequence with Bowetie2")
    fmsystem.run_command(command, mutex_log=logger)


def samtools_index_fasta(input_fasta: str,
        output_index: str,
        logger: fmruffus.RuffusLog = None,
        ) -> None:
    """Index a list of .fasta files using samtools faidx"""

    # Construct the samtools command
    command = ['samtools', 'faidx', input_fasta]

    # Run the command
    logger.info("Indexing fasta sequence with samtools")
    fmsystem.run_command(command, mutex_log=logger)


def bowtie2_align(fwd_read, rev_read, indices, out_paths, threads=None,
        maxins=None, no_discordant=None, no_mixed=None, preset=None):
    """
    Python wrapper around bowtie2

    Parameters
    ----------
        fwd_read, rev_read = List of paths to fastq files (.fastq/.fastq.gz)
        indices = List of paths to the reference fasta indices
        out_paths = List of output file paths (Default: input[-.fastq].align.sam)
        All other arguments refer directly to bowtie2 arguments
    Returns
    -------
       The exit codes for each bowtie2 run as a list
    """

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
    temp_out = add_suffix(out_paths, '.part')

    # Run bowtie
    exit_codes = []
    for r1, r2, idx, tmp, out in zip(fwd_read, rev_read, indices, temp_out,
            out_paths):

        # Bowtie2 sometimes segfaults with gzipped fastq files, so we gunzip
        # them
        if get_final_suffix(r1) == '.gz':
            r1 = '<(gunzip -c ' + r1 + ')'
        if get_final_suffix(r2) == '.gz':
            r2 = '<(gunzip -c ' + r2 + ')'

        # Construct command
        command = (['bowtie2', '-1', r1, '-2', r2, '-x', idx, '-S', tmp] +
                opt_command)

        # Run command
        exit_code = run_command(command, ("bowtie2." +
            remove_suffix(os.path.basename(idx))))[0]
        exit_codes.append(exit_code)

        # Remove the temp file
        move(tmp, out)

    return exit_codes
