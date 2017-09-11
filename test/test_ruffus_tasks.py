import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
import fmbiopy.ruffus_tasks as ruffus_tasks
from get_dat import get_dat
from glob import glob
import os
import pytest
import tempfile

def check_ruffus_task(task, input_files, output_files, to_delete):
    """Check that a Ruffus function runs as expected

    Specifically checks that the requested output_files are created, and that
    the exit code is zero"""

    with fmsystem.delete(to_delete):
        # Some ruffus tasks might not have a logger argument.
        exit_code = task(input_files, output_files)

        assert exit_code == 0
        if isinstance(output_files, str):
            assert os.path.exists(output_files)
        else:
            for f in output_files:
                assert os.path.exists(f)

def test_samtools_index_fasta():
    assembly = get_dat()['assemblies'][0]
    root_index = fmpaths.remove_suffix(assembly)
    expected_index = assembly + '.fai'
    fmsystem.silent_remove(expected_index)

    check_ruffus_task(
            ruffus_tasks.samtools_index_fasta, assembly, expected_index,
            expected_index)


def test_bowtie_index_fasta():
    assembly = get_dat()['assemblies'][0]
    root_index = fmpaths.remove_suffix(assembly)
    bowtie_indices = fmpaths.get_bowtie2_indices(root_index)

    # Cant use check_ruffus_task here because the outputs do not correspond to
    # files
    with fmsystem.delete(bowtie_indices):
        exit_code = ruffus_tasks.bowtie_index_fasta(assembly, root_index)
        assert exit_code == 0
        for f in bowtie_indices:
            assert os.path.exists(f)

def test_gunzip():
    reads = get_dat()['fwd_reads'][0]
    gunzipped_output = fmpaths.remove_suffix(reads)

    check_ruffus_task(ruffus_tasks.gunzip, reads, gunzipped_output,
            gunzipped_output)

def test_gzip():
    tmp = tempfile.NamedTemporaryFile()
    gzipped_output = fmpaths.add_suffix(tmp.name, '.gz')

    with fmsystem.delete(gzipped_output):
        check_ruffus_task(ruffus_tasks.gzip, tmp.name, gzipped_output,
                gzipped_output)

def test_paired_bowtie2_align():
    fwd_reads = get_dat()['fwd_reads'][0]
    rev_reads = get_dat()['rev_reads'][0]
    gunzip_fwd = fmpaths.remove_suffix(fwd_reads)
    gunzip_rev = fmpaths.remove_suffix(rev_reads)
    ruffus_tasks.gunzip(fwd_reads, gunzip_fwd)
    ruffus_tasks.gunzip(rev_reads, gunzip_rev)

    assembly = get_dat()['assemblies'][0]
    root_name = fmpaths.remove_suffix(fwd_reads, 2)
    output_index = fmpaths.remove_suffix(assembly)
    bowtie2_indices = fmpaths.get_bowtie2_indices(output_index)
    ruffus_tasks.bowtie_index_fasta(assembly, output_index)
    output_bam = fmpaths.add_suffix(root_name, '.bam')
    output_bai = fmpaths.add_suffix(output_bam, '.bai')
    input_files = (gunzip_fwd, gunzip_rev, output_index)
    output_files = bowtie2_indices + [output_bam, gunzip_fwd,
                                      gunzip_rev, output_bai]
    pytest.set_trace()
    check_ruffus_task(
            ruffus_tasks.paired_bowtie2_align, input_files, output_bam,
            output_files)
