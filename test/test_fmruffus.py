"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""

from pathlib import Path
import pytest

from fmbiopy.fmclass import list_classes
from fmbiopy.fmlist import flatten
from fmbiopy.fmpaths import (
        any_exist,
        as_paths,
        as_strs,
        )
from fmbiopy.fmruffus import *


@pytest.fixture(scope='module')
def instance_of(example_file):
    """Return example RuffusTask instances"""
    def _make_test_instance(cls, size):

        input_example = [example_file(t, size) for t in cls.input_type]
        input_example = flatten(input_example)

        # If we have multiple inputs, the output suffix is added to
        # the first input as in ruffus
        output_example = []
        for typ in cls.output_type:
            if typ == '':
                output_example.append(
                        input_example[0].with_suffix(typ))
            elif typ == 'SAME':
                # If output is same type create a new directory to place the
                # output
                parent_dir = input_example[0].parent
                outdir = parent_dir / 'tmp'
                outdir.mkdir(exist_ok=True)
                output_example.append(outdir / input_example[0].name)
            else:
                output_example.append(
                        input_example[0].with_suffix('.' + typ))

        return cls(
                as_strs(input_example),
                as_strs(output_example))
    return _make_test_instance


class TestRuffusLog(object):
    def test_non_existant_path(self):
        with pytest.raises(FileNotFoundError):
            RuffusLog("foo", Path("bar/bar.log"))

    def test_normal_unbuffered_usage(self, gen_tmp, tmpdir):
        tmp = gen_tmp(directory=tmpdir)
        ruflog = RuffusLog('foo', tmp, buffered=False)
        ruflog.write("Test")
        # Check filesize non zero
        assert Path(tmp).read_text()

    def test_header(self, gen_tmp, tmpdir):
        tmp = gen_tmp(directory=tmpdir)
        ruflog = RuffusLog('foo', tmp, buffered=False)
        head = "HEADER"
        ruflog.write_header(head)
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert '=' * 80 in lines[0]
        assert head in lines[1]

    def test_normal_buffered_usage(self, gen_tmp, tmpdir):
        tmp = gen_tmp(directory=tmpdir)
        ruflog = RuffusLog('foo', tmp, buffered=True)
        ruflog.write("bar")
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert not lines
        ruflog.flush()
        with open(tmp, 'r') as f:
            lines = [line for line in f][:2]
        assert lines


@pytest.fixture(
        params=list_classes(
            'fmruffus',
            'fmbiopy',
            of_type=['RuffusTask'],
            exclude=['RuffusLog', 'RuffusTask']))
def task(request, instance_of):
    task_class = request.param
    ruffus_task = instance_of(task_class, 'tiny')
    yield ruffus_task
    ruffus_task._cleanup()


class TestAllTasks(object):

    def test_command_not_none(self, task):
        assert task._command is not None

    def test_command_is_listlike(self, task):
        assert not isinstance(task._command, str)
        assert task._command[0]

    def test_command_has_no_spaces(self, task):
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
            assert not any_exist(as_paths(task._input_files))


def test_apply_symlink_produces_expected_output(full_dir):
    symlink_all = apply(SymlinkInputs)
    inputs = full_dir.glob('*')
    output_ids = ['tar.x', 'sar.x', 'lar.x']
    outputs = [full_dir / out for out in output_ids]
    symlink_all(as_strs(inputs), as_strs(outputs))
    for path in outputs:
        assert path.exists()
