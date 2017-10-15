"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""

from pathlib import Path
import pytest

import fmbiopy.fmclass as fmclass
import fmbiopy.fmlist as fmlist
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus
import fmbiopy.fmtest as fmtest


@pytest.fixture
def instance_of(example_file):
    """Return example RuffusTask instances"""
    def _make_test_instance(cls, size):

        input_example = [example_file(t, size) for t in cls.input_type]
        input_example = fmlist.flatten(input_example)

        # If we have multiple inputs, the output suffix is added to
        # the first input as in ruffus
        output_example = []
        for typ in cls.output_type:
            if cls.output_type == ['']:
                output_example.append(
                        input_example[0].with_suffix(typ))
            else:
                output_example.append(
                        input_example[0].with_suffix('.' + typ))

        return cls(
                fmpaths.as_str(input_example),
                fmpaths.as_str(output_example))
    return _make_test_instance


class TestRuffusLog(object):
    def test_non_existant_path(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("foo", Path("bar/bar.log"))

    def test_normal_unbuffered_usage(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog('foo', tmp, buffered=False)
        ruflog.write("Test")
        # Check filesize non zero
        assert Path(tmp).read_text()

    def test_header(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog('foo', tmp, buffered=False)
        head = "HEADER"
        ruflog.write_header(head)
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert '=' * 80 in lines[0]
        assert head in lines[1]

    def test_normal_buffered_usage(self, tmpdir):
        tmp = fmtest.gen_tmp(directory=tmpdir)
        ruflog = fmruffus.RuffusLog('foo', tmp, buffered=True)
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
            assert Path(task._output_files).exists()
        else:
            for path in task._output_files:
                assert Path(path).exists()

    def test_run_command_produces_zero_exit_code(self, task):
        for code in task.exit_code:
            assert code == 0

    def test_inplace_tasks_delete_their_inputs(self, task):
        if task._inplace:
            assert not fmpaths.any_exist(fmpaths.as_path(task._input_files))
