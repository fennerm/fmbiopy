from plumbum import local

from test.helpers import validate_bam_file


def test_sam_to_bam(trimmed_sam, sandbox):
    bam = sandbox / trimmed_sam.with_suffix('.bam')
    sam_to_bam = local["bin/sam_to_bam.py"]
    sam_to_bam(trimmed_sam, bam)
    assert bam.exists()
    validate_bam_file(bam)
