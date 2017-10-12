"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""

import os

import pytest

import fmbiopy.fmcheck as fmcheck
import fmbiopy.fmclass as fmclass
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmtest as fmtest
from fmbiopy.fmtest import instance_of


class TestRuffusLog(object):
    def test_non_existant_path(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("foo", "bar/bar.log")

    def test_normal_unbuffered_usage(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=False)
        ruflog.write("Test")
        assert os.path.getsize(tmp) > 0

    def test_header(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=False)
        head = "HEADER"
        ruflog.write_header(head)
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert '=' * 80 in lines[0]
        assert head in lines[1]

    def test_normal_buffered_usage(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog("foo", tmp, buffered=True)
        ruflog.write("bar")
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert not lines
        ruflog.flush()
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert lines


@pytest.fixture(
        params=fmclass.list_classes(
            'fmbiopy.fmruffus',
            exclude=['RuffusLog', 'RuffusTask']))
def task(request, instance_of):
    task_class = request.param
    ruffus_task = instance_of(task_class, 'tiny')
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
        for code in task.exit_code:
            assert code == 0

    def test_inplace_tasks_delete_their_inputs(self, task):
        if task._inplace:
            assert not fmcheck.any_exist(task._input_files)
