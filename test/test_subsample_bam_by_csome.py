from plumbum import local
from plumbum.cmd import samtools
from pytest import fixture

from fmbiopy.fmlist import exclude_blank


@fixture(name="subsampled_bam")
def gen_subsampled_bam(dat, tmpdir):
    inbam = dat["small"]["bam"][0]
    outbam = tmpdir / "a.bam"
    n = "100"
    subsample_bam = local["bin/subsample_bam_by_csome.py"]
    subsample_bam("-o", outbam, "-n", n, inbam)
    samtools("index", outbam)
    return (inbam, outbam)


def get_nonempty_csomes(bam):
    list_csomes = local["bin/list_csomes"]
    nonempty_csomes = list_csomes(bam)
    nonempty_csomes = exclude_blank(nonempty_csomes.split('\n'))
    return nonempty_csomes


def test_subsample_bam_by_csome(subsampled_bam):
    input_csomes = get_nonempty_csomes(subsampled_bam[0])
    output_csomes = get_nonempty_csomes(subsampled_bam[1])
    assert len(input_csomes) > len(output_csomes)
    # Test if subset
    assert set(output_csomes) < set(input_csomes)
