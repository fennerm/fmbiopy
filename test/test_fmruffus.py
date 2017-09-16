"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""

import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
from get_dat import get_dat
from glob import glob
import os
import pytest
from tempfile import NamedTemporaryFile as tmp

class TestRuffusLog(object):
    def test_empty_name_returns_value_error(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("", tmp().name)

    def test_non_existant_path(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("foo", "bar/bar.log")

    def test_normal_usage(self):
        temp = tmp()
        ruflog = fmruffus.RuffusLog("foo", temp.name)
        ruflog.write("Test")
        assert os.path.getsize(temp.name) > 0

    def test_header(self):
        temp = tmp()
        ruflog = fmruffus.RuffusLog("foo", temp.name)
        head = "HEADER"
        ruflog.write_header(head)
        with open(temp.name, 'r') as f:
            lines = [line for line in f][:2]
        assert '-' * 30 in lines[0]
        assert head in lines[1]

@pytest.fixture
def ruffus_task(fasta_paths, samtools_index_paths):
    return fmruffus.RuffusTask(fasta_paths, samtools_index_paths,
                               log_results=True)


class TestRuffusTask(object):
    def test_initialization(self, ruffus_task):
        assert ruffus_task.logger.log is not None


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
            for f in to_delete:
                assert os.path.exists(f)

def test_samtools_index_fasta():
    assembly = get_dat()['assemblies'][0]
    root_index = fmpaths.remove_suffix(assembly)
    expected_index = assembly + '.fai'
    fmsystem.silent_remove(expected_index)

    check_ruffus_task(
            fmruffus.samtools_index_fasta, assembly, expected_index,
            expected_index)


def test_bowtie_index_fasta():
    assembly = get_dat()['assemblies'][0]
    root_index = fmpaths.remove_suffix(assembly)
    bowtie_indices = fmpaths.get_bowtie2_indices(root_index)

    # Cant use check_ruffus_task here because the outputs do not correspond to
    # files
    with fmsystem.delete(bowtie_indices):
        exit_code = fmruffus.bowtie_index_fasta(assembly, root_index)
        assert exit_code == 0
        for f in bowtie_indices:
            assert os.path.exists(f)


def test_gunzip():
    reads = get_dat()['fwd_reads'][0]
    gunzipped_output = fmpaths.remove_suffix(reads)

    check_ruffus_task(fmruffus.gunzip, reads, gunzipped_output,
            gunzipped_output)


def test_gzip():
    tmp = tempfile.NamedTemporaryFile()
    gzipped_output = fmpaths.add_suffix(tmp.name, '.gz')

    with fmsystem.delete(gzipped_output):
        check_ruffus_task(fmruffus.gzip, tmp.name, gzipped_output,
                gzipped_output)


def test_paired_bowtie2_align():
    fwd_reads = get_dat()['fwd_reads'][0]
    rev_reads = get_dat()['rev_reads'][0]
    gunzip_fwd = fmpaths.remove_suffix(fwd_reads)
    gunzip_rev = fmpaths.remove_suffix(rev_reads)
    fmruffus.gunzip(fwd_reads, gunzip_fwd)
    fmruffus.gunzip(rev_reads, gunzip_rev)

    assembly = get_dat()['assemblies'][0]
    root_name = fmpaths.remove_suffix(fwd_reads, 2)
    output_index = fmpaths.remove_suffix(assembly)
    bowtie2_indices = fmpaths.get_bowtie2_indices(output_index)
    fmruffus.bowtie_index_fasta(assembly, output_index)
    output_bam = fmpaths.add_suffix(root_name, '.bam')
    output_bai = fmpaths.add_suffix(output_bam, '.bai')
    input_files = (gunzip_fwd, gunzip_rev, output_index)
    output_files = bowtie2_indices + [output_bam, gunzip_fwd,
                                      gunzip_rev, output_bai]
    check_ruffus_task(
            fmruffus.paired_bowtie2_align, input_files, output_bam,
            output_files)
