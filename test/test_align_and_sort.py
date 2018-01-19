from plumbum import local
from pytest import mark

from fmbiopy.fmpaths import is_empty
from fmbiopy.fmtest import assert_script_produces_files


@mark.parametrize("threads", [None, 2])
def test_align_and_sort(threads, dat, sandbox):
    fwd_reads = dat['tiny']['fwd_reads'][0]
    rev_reads = dat['tiny']['rev_reads'][0]
    ref = dat['tiny']['assemblies'][0]
    idx = ref.with_suffix('')
    unpaired = dat['tiny']['fwd_reads'][1]
    out_bam = sandbox / 'out.bam'

    args = ['-1', fwd_reads, '-2', rev_reads, '-x', idx, '-U', unpaired]
    if threads:
        args = args + ['-p', threads]

    assert_script_produces_files(script='bin/align_and_sort.py', args=args,
                                 output=[out_bam],
                                 redirect=out_bam)
