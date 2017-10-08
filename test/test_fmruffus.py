"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""

import os
import pytest

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmclass as fmclass
import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmtest as fmtest


class TestRuffusLog(object):
    def test_non_existant_path(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("foo", "bar/bar.log")

    def test_normal_unbuffered_usage(self):
        tmp = fmtest.gen_tmp()
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=False)
        ruflog.write("Test")
        assert os.path.getsize(tmp) > 0

    def test_header(self):
        tmp = fmtest.gen_tmp()
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=False)
        head = "HEADER"
        ruflog.write_header(head)
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert '=' * 80 in lines[0]
        assert head in lines[1]

    def test_normal_buffered_usage(self):
        tmp = fmtest.gen_tmp()
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=True)
        ruflog.write("bar")
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert not lines
        ruflog.flush()
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert lines



def get_example_file(filetype):
    if filetype == 'fasta':
        return fmtest.get_dat()['assemblies'][0]
    elif filetype == ('fastq', 'fastq'):
        return (fmtest.get_dat()['fwd_reads'][0],
                fmtest.get_dat()['rev_reads'][0])
    elif filetype == 'fastq':
        return fmtest.get_dat()['fwd_reads'][0]
    elif filetype == 'fai':
        return fmtest.get_dat()['faindices'][0]
    elif filetype == 'sam':
        return fmtest.get_dat()['sam'][0]
    elif filetype == 'bam':
        return fmtest.get_dat()['bam'][0]
    elif filetype == 'gz':
        return fmtest.get_dat()['fwd_reads'][0]
    return fmtest.gen_tmp(empty=False)


def get_test_instance(task):
    """Given a RuffusTask name, return an instance of the task"""
    input_example = [get_example_file(t) for t in task.input_type]
    input_example = fmlist.flatten(input_example)

    if task.output_type == ['']:
        output_example = fmpaths.remove_suffix(input_example[0])
    else:
        output_suffix = '.' + task.output_type[0]
        output_example = [
                fmpaths.add_suffix(path, output_suffix) for path
                in input_example]
        # Match the length of the example output to the actual number of
        # outputs
        output_example = output_example[0:len(task.output_type)]

    return task(input_example, output_example)


@pytest.fixture(
        scope="module",
        params=fmclass.list_classes(
            'fmbiopy.fmruffus',
            exclude=[fmruffus.RuffusLog, fmruffus.RuffusTask]))
def task(request):
    ruffus_task = get_test_instance(request.param)
    yield ruffus_task
    ruffus_task._cleanup()


class TestAllTasks(object):
    def test_initialization(self, task):
        assert task._logger.log is not None
        assert task._construct_command is not None
        assert task._run_command is not None

    def test_construct_command_not_none(self, task):
        assert task._command is not None

    def test_construct_command_is_listlike(self, task):
        assert not isinstance(task._command, str)
        assert task._command[0]

    def test_construct_command_has_no_spaces(self, task):
        for word in task._command:
            assert ' ' not in word

    def test_run_command_produces_expected_output(self, task):
        if isinstance(task._output_files, str):
            os.path.exists(task._output_files)
        else:
            for path in task._output_files:
                assert os.path.exists(path)

    def test_run_command_produces_zero_exit_code(self, task):
        assert task.exit_code == 0

    def test_inplace_tasks_delete_their_inputs(self, task):
        if task._inplace:
            assert not fmcheck.any_exist(task._input_files)
