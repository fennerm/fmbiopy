from plumbum import local
from fmbiopy.fmtest import validate_bam_file


def test_sam_to_bam(sam, sandbox):
    bam = sandbox / sam.with_suffix('.bam')
    sam_to_bam = local["bin/sam_to_bam.py"]
    sam_to_bam(sam, bam)
    assert bam.exists()
    validate_bam_file(bam)
