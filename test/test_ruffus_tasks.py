import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmsystem as fmsystem
import fmbiopy.ruffus_tasks as ruffus_tasks
from get_dat import get_dat
from glob import glob
import inspect
import os
import pytest
import tempfile

def check_ruffus_task(task, input_files, output_files):
    """Check that a Ruffus function runs as expected

    Specifically checks that the requested output_files are created, and that
    the exit code is zero"""

    log = fmruffus.RuffusLog("Temp", "tmp.log")

    # Handle deletion differently depending on whether input is a string
    if isinstance(output_files, str):
        to_delete = [output_files, log.config['file_name']]
    else:
        to_delete = output_files + [log.config['file_name']]

    with fmsystem.delete(to_delete):
        # Some ruffus tasks might not have a logger argument.
        task_param = inspect.getfullargspec(task)
        if 'logger' in task_param:
            exit_code = task(input_files, output_files, logger=log)
        else:
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
            ruffus_tasks.samtools_index_fasta, assembly, expected_index)


def test_bowtie_index_fasta():
    assembly = get_dat()['assemblies'][0]
    root_index = fmpaths.remove_suffix(assembly)
    bowtie_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',
            '.rev.2.bt2']
    expected_indices = sorted([root_index + suf for suf in bowtie_suffixes])

    # Cant use check_ruffus_task here because the outputs do not correspond to
    # files
    with fmsystem.delete(expected_indices):
        exit_code = ruffus_tasks.bowtie_index_fasta(assembly, root_index)
        assert exit_code == 0
        for f in expected_indices:
            assert os.path.exists(f)

def test_gunzip():
    reads = get_dat()['fwd_reads'][0]
    gunzipped_output = fmpaths.remove_suffix(reads)

    check_ruffus_task(ruffus_tasks.gunzip, reads, gunzipped_output)

def test_gzip():
    tmp = tempfile.NamedTemporaryFile()
    gzipped_output = fmpaths.add_suffix(tmp.name, '.gz')

    with fmsystem.delete(gzipped_output):
        check_ruffus_task(ruffus_tasks.gzip, tmp.name, gzipped_output)

"""
def test_bowtie2_align():
    fwd_reads = get_dat()['fwd_reads'][0]
    rev_reads = get_dat()['rev_reads'][0]
    assembly = testdat['assemblies'][0]
    root_name = fmpaths.remove_suffix(fwd_reads, 2)
    output_index = fmpaths.remove_suffix(assembly)
    bowtie_index_fasta(assembly, output_index)
    output_sam = fmpaths.add_suffix(root_name, '.sam')

    with fmsystem.delete([output_index, output_sam]):
        input_files = (fwd_reads, rev_reads, output_index)
        exit_code = ruffus_tasks.paired_bowtie2_align(input_files, output_sam)
        assert os.path.exists(output_sam)
        assert exit_code == 0
"""
