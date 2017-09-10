#!/usr/bin/env/ python

import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
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
