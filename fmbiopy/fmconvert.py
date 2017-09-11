"""Functions for converting file types"""

import fmbiopy.fmsystem as fmsystem
import os


def sam_to_bam(sam: str, bam: str) -> None:
    """Convert a sam file to a sorted, indexed bam file"""

    # Convert to bam and sort
    sam_to_bam_command = ('samtools view -bS ' + sam +
                          ' | samtools sort -o ' + bam + ' -')
    os.system(sam_to_bam_command)

    # Samtools index
    index_sam_command = ('samtools index ' + bam)
    os.system(index_sam_command)

    # Delete the sam file
    if os.path.exists(bam):
        fmsystem.silent_remove(sam)
    else:
        raise OSError('sam_to_bam system command failed')
