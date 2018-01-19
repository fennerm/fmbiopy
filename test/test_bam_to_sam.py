from fmbiopy.fmtest import (
    assert_script_produces_files,
    validate_bam_file,
    )

def test_bam_to_sam_usage(tiny_indexed_bam, tmpdir):
    outsam = tmpdir / "out.sam"
    assert_script_produces_files(script="bin/bam_to_sam",
                                 args=[tiny_indexed_bam], output=[outsam],
                                 redirect=outsam)
    validate_bam_file(outsam)
