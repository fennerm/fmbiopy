#!/usr/bin/env/ python

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmsystem as fmsystem
from glob import glob
import os
from ruffus.proxy_logger import AcquirerProxy
from ruffus.proxy_logger import LoggerProxy
from typing import Dict
from typing import Sequence
from typing import Tuple

def samtools_index_fasta(input_fasta: str,
        output_indices: str,
        param: Dict[str, str],
        logger: LoggerProxy,
        logging_mutex: AcquirerProxy
        ) -> None:
    """Index a list of .fasta files using samtools faidx"""

    # Convert the parameter dictionary to a flat list
    parsed_param = fmsystem.dict_to_list(param)

    # Construct the samtools command
    command = ['samtools', 'faidx'] + parsed_param + [input_fasta]

    # Run the command
    fmsystem.run_command(command, "build_index", False)
