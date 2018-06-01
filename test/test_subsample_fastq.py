from fmbiopy.paths import size
from plumbum import local

def test_subsample_fastq(dat):
    fwd = dat['tiny']['fwd_reads'][0]
    rev = dat['tiny']['rev_reads'][0]
    subfwd = fwd.with_suffix('.sub.fq.gz')
    subrev = rev.with_suffix('.rev.fq.gz')
    subsample = local['bin/subsample_fastq.py']
    subsample['-n', 5, fwd, rev, subfwd, subrev]()
    for inp, out in [(fwd, subfwd), (rev, subrev)]:
        assert out.exists()
        assert size(out) < size(inp)
        out.delete()
