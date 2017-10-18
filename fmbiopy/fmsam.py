"""Functions for manipulating Sam and Bam files"""

from pathlib import Path

from fmbiopy.fmpaths import is_empty
from fmbiopy.fmsystem import run_command


def sam_to_bam(sam: Path, bam: Path)-> None:
    """Convert a .sam file to a sorted, indexed .bam file

    Parameters
    ----------
    sam
        Path to a .sam file
    bam
        Location to output the .bam file
    """

    # Convert to bam and sort
    sam_to_bam_command = [
            'samtools', 'view', '-bS', str(sam), '|', 'samtools', 'sort', '-o',
            str(bam), '-']

    run_command(sam_to_bam_command, shell=True)

    # Samtools index
    index_sam_command = ['samtools', 'index', str(bam)]
    run_command(index_sam_command, shell=True)

    # Delete the sam file
    if bam.exists() and not is_empty(bam):
        sam.unlink()
    else:
        raise OSError('sam_to_bam system command failed')
