from plumbum import local

def test_sam_to_bam(dat, sandbox):
    sam = dat['tiny']['sam'][0]
    bam = sandbox / sam.with_suffix('.bam')
    sam_to_bam = local["bin/sam_to_bam.py"]
    sam_to_bam(sam, bam)
    assert bam.exists()
