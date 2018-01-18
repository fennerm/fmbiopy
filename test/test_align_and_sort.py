from plumbum import local
from pytest import mark

from fmbiopy.fmpaths import is_empty


@mark.parametrize("threads", [None, 2])
def test_align_and_sort(threads, dat, sandbox):
    fwd_reads = dat['tiny']['fwd_reads'][0]
    rev_reads = dat['tiny']['rev_reads'][0]
    ref = dat['tiny']['assemblies'][0]
    idx = ref.with_suffix('')
    unpaired = dat['tiny']['fwd_reads'][1]
    out_bam = sandbox / 'out.bam'
    align_and_sort = local['bin/align_and_sort.py']

    if threads:
        (align_and_sort[
            '-1', fwd_reads, '-2', rev_reads, '-p', threads, '-x', idx, '-U',
            unpaired] > out_bam)()
    else:
        (align_and_sort[
            '-1', fwd_reads, '-2', rev_reads, '-x', idx, '-U', unpaired] >
         out_bam)()

    assert out_bam.exists()
    assert not is_empty(out_bam)
