from plumbum import local
from plumbum.cmd import samtools
from pytest import fixture

from fmbiopy.iter import exclude_blank
from fmbiopy.system import capture_stdout


@fixture(name="subsampled_bam")
def gen_subsampled_bam(dat, tmpdir):
    inbam = dat["small"]["bam"][0]
    outbam = tmpdir / "a.bam"
    n = "100"
    subsample_bam = local["bin/subsample_bam_by_csome.py"]
    subsample_bam("-o", outbam, "-n", n, inbam)
    samtools("index", outbam)
    return (inbam, outbam)


@fixture
def list_csomes():
    return local['bin/list_csomes']


def test_subsample_bam_by_csome(list_csomes, subsampled_bam):
    input_csomes = capture_stdout(list_csomes[subsampled_bam[0]])
    output_csomes = capture_stdout(list_csomes[subsampled_bam[1]])
    assert len(input_csomes) > len(output_csomes)
    # Test if subset
    assert set(output_csomes) < set(input_csomes)
