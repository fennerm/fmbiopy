""" Various broadly useful Ruffus task definitions

Ruffus: http://www.ruffus.org.uk/
"""

import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
import logging
import os
from shutil import move
import typing

def _run_ruffus_command(
        command: typing.List[str],
        logger: fmruffus.RuffusLog = None,
        log_results: bool = True) -> int:
    """Helper function for running a command with or without logging"""
    if log_results:
        exit = fmsystem.run_command(
                command, mutex_log=logger, log_stdout=False,
                log_stderr=False)[0]
    else:
        exit = fmsystem.run_command(command, mutex_log=logger)[0]
    return exit


def bowtie_index_fasta(
        input_fasta: str,
        output_prefix: str,
        param: typing.List[str] = [],
        logger: fmruffus.RuffusLog = None,
        log_results: bool = False
        ) -> int:
    """Index a .fasta file using samtools faidx

    Parameters
    ----------
    input_fasta
        Path to fasta file to be indexed
    output_prefix
        Prefix for the output Bowtie2 indices
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    logger
        The Ruffus logging instance
    log_results
        If False, results will not be logged to file

    Returns
    -------
    The error code of the process
    """

    # Construct the samtools command
    command = ['bowtie2-build'] + param + [input_fasta, output_prefix]

    if logger:
        logger.write_header("Indexing " + input_fasta + " with Bowtie2")

    # Run the command
    return _run_ruffus_command(command, logger, log_results)


def samtools_index_fasta(input_fasta: str,
        output_index: str,
        logger: fmruffus.RuffusLog = None,
        log_results: bool = False
        ) -> int:
    """Index a list of .fasta files using samtools faidx

    Parameters
    ----------
    input_fasta
        Path to fasta file to be indexed
    output_index
        Path to the output .fai file
    logger
        The Ruffus logging instance
    log_results
        If False, results will not be logged to file
    Returns
    -------
    The error code of the process

    """

    # Construct the samtools command
    command = ['samtools', 'faidx', input_fasta]

    if logger:
        logger.write_header("Indexing " + input_fasta + " with samtools")

    # Run the command
    exit = _run_ruffus_command(command, logger, log_results)
    return exit

def gunzip(
        input_file: str,
        output_file: str,
        param: typing.List[str] = []) -> int:
    """Gunzip a file and keep the original

    Parameters
    ----------
    input_file
        Path to input file
    output_file
        Path to the gunzipped output
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    """
    command = ['gunzip', '-c'] + param + [input_file]
    exit, stdout, _ = fmsystem.run_command(
            command, log_stdout=False, log_stderr=False)
    with open(output_file, "w") as f:
        f.write(stdout)
    return exit

def paired_bowtie2_align(
        input_files: typing.Tuple[str, str, str],
        output_sam: str,
        param: typing.List[str]) -> int:
    """Align a pair of fastq files to a fasta file using bowtie2

    Parameters
    ----------
    input_files
        A Tuple of the form (Forward reads, Reverse reads, Bowtie2 index)
    output_index
        Path to the output .sam file
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    logger
        The Ruffus logging instance
    log_results
        If False, results will not be logged to file
    Returns
    -------
    The exit code of the process
    """

    fwd_reads = input_files[0]
    rev_reads = input_files[1]
    fasta = input_files[2]
    # Bowtie2 sometimes segfaults with gzipped fastq files, so we gunzip
    # them
    if fmpaths.get_final_suffix(fwd_reads) == '.gz':
        gunzipped_fwd_reads = fmpaths.remove_suffix(fwd_reads)
        raise OSError

        fwd_reads = '<(gunzip -c ' + r1 + ')'
    if fmpaths.get_final_suffix(rev_reads) == '.gz':
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
