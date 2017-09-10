#!/usr/bin/env/ python

import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem

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
