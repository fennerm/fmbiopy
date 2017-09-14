""" Various broadly useful Ruffus task definitions

Ruffus: http://www.ruffus.org.uk/
"""

import fmbiopy.fmsam as fmsam
import fmbiopy.fmpaths as fmpaths
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
        logger.write_header(['Running:'] + command)
        exit_code = fmsystem.run_command(
                command, mutex_log=logger, log_stdout=False,
                log_stderr=False)[0]
    else:
        exit_code = fmsystem.run_command(command, mutex_log=logger)[0]
    return exit_code


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

    # Run the command
    exit_code = _run_ruffus_command(command, logger, log_results)
    return exit_code


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

    # Run the command
    exit_code = _run_ruffus_command(command, logger, log_results)
    return exit_code


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
    exit_code, stdout, _ = fmsystem.run_command(
            command, log_stdout=False, log_stderr=False)
    with open(output_file, "w") as f:
        f.write(stdout)
    return exit_code


def gzip(
        input_file: str,
        output_file: str,
        param: typing.List[str] = []) -> int:
    """Gzip a file

    Parameters
    ----------
    input_file
        Path to input file
    output_file
        Path to the gzipped output
    param
        A list of bash parameters. E.g ['-x', 'foo', '--long', 'bar']
    """
    command = ['gzip'] + param + [input_file]
    exit_code = fmsystem.run_command(
            command, log_stdout=False, log_stderr=False)[0]
    return exit_code

def paired_bowtie2_align(
        input_files: typing.Tuple[str, str, str],
        output_bam: str,
        param: typing.List[str] = [],
        logger: fmruffus.RuffusLog = None,
        log_results: bool = False) -> int:
    """Align a pair of fastq files to a fasta file using bowtie2

    Parameters
    ----------
    input_files
        A Tuple of the form (Forward reads, Reverse reads, Bowtie2 index)
    output_bam
        Path to the output .bam file
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
    bowtie2_index = input_files[2]

    # Construct command
    output_sam = fmpaths.replace_suffix(output_bam, '.bam', '.sam')
    command = ['bowtie2', '-1', fwd_reads, '-2', rev_reads, '-x',
            bowtie2_index, '-S', output_sam]

    # Run command
    exit_code = _run_ruffus_command(command, logger, log_results)

    # Convert to a sorted Bam
    fmsam.sam_to_bam(output_sam, output_bam)

    return exit_code
