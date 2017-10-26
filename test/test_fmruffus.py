"""Classes and functions for use in Ruffus pipeline infrastructure

Ruffus: http://www.ruffus.org.uk/
"""
from collections import namedtuple
from pathlib import Path
from pytest import (
        fixture,
        mark,
        raises,
        )
from uuid import uuid4

from fmbiopy.fmclass import list_classes
from fmbiopy.fmlist import flatten
from fmbiopy.fmpaths import (
        add_suffix,
        any_exist,
        as_paths,
        as_strs,
        )
from fmbiopy.fmruffus import *


class TestRuffusLog(object):
    @fixture
    def unbuffered_ruflog(self, nonexistant_path):
        return RuffusLog(
                name='foo',
                location=nonexistant_path,
                buffered=False,
                overwrite=True)

    def filename(self, ruflog):
        return Path(ruflog.config['file_name'])

    def test_possible_path_inputs(self, poss_paths):
        if poss_paths[0] == 'nonexistant_parent':
            with raises(FileNotFoundError):
                RuffusLog('foo', poss_paths[1])
        elif poss_paths[0] == 'tmpdir':
            with raises(IsADirectoryError):
                RuffusLog('foo', poss_paths[1])
        else:
            RuffusLog('foo', poss_paths[1])


    @mark.parametrize('buffered', [True, False])
    @mark.parametrize('overwrite', [True, False])
    def test_optional_parameters(self, buffered, overwrite, nonexistant_path):
        ruflog = RuffusLog(
                name='foo',
                location=nonexistant_path,
                buffered=buffered,
                overwrite=overwrite)
        ruflog.write("Test")

        if buffered:
            assert not nonexistant_path.exists()
            ruflog.flush()
            content = nonexistant_path.read_text()
            assert content
        else:
            content = nonexistant_path.read_text()
            assert content

        ruflog2 = RuffusLog(
                name='foo',
                location=nonexistant_path,
                buffered=buffered,
                overwrite=overwrite)
        if overwrite:
            assert not nonexistant_path.exists()
        else:
            content = nonexistant_path.read_text()
            assert content

    @mark.parametrize('subheader', [True, False])
    def test_header(self, subheader, unbuffered_ruflog):
        head = "HEADER"
        unbuffered_ruflog.write_header(head, subheader)
        location = self.filename(unbuffered_ruflog)
        with location.open('r') as f:
            lines = [line for line in f][:3]
        assert head in lines[1]
        if subheader:
            assert '-' * 80 in lines[0]
            assert '-' * 80 in lines[2]
        else:
            assert '=' * 80 in lines[0]
            assert '=' * 80 in lines[2]


class TestAllTasks(object):
    @fixture(scope='module')
    def get_example_files(self, example_file, sandbox):
        """Return example RuffusTask instances"""
        def _make_test_instance(cls, size):

            input_example = [example_file(t, size) for t in cls.input_type]
            input_example = flatten(input_example)

            # If we have multiple inputs, the output suffix is added to
            # the first input as in ruffus
            output_example = []
            for typ in cls.output_type:
                input_stem = sandbox / input_example[0].stem
                if typ == '':
                    output_example.append(input_stem)
                elif typ == 'SAME':
                    tempdir = sandbox / uuid4().hex
                    output_example.append(tempdir / input_example[0].name)
                else:
                    suffix = '.' + typ
                    output_example.append(add_suffix(input_stem, suffix))

            return (as_strs(input_example), as_strs(output_example))
        return _make_test_instance

    @fixture(
            params=list_classes(
                'fmruffus',
                'fmbiopy',
                of_type=['RuffusTransform'],
                exclude=[
                    'RuffusLog', 'RuffusTask', 'RuffusTransform',
                    'BuildCentrifugeDB']))
    def task(self, request, cd, get_example_files):
        cls = request.param
        inp, out = get_example_files(cls, 'tiny')
        instance = cls()
        instance.run(inp, out)
        task_info = namedtuple('task_info', 'cls instance inp out')
        yield task_info(cls, instance, inp, out)
        instance.cleanup()

    def test_param_are_added_to_command_or_exception_raised(self, task):
        try:
            param=['--nonexistant', 'foo']
            param_instance = task.cls(param=param)
            param_instance.run(task.inp, task.out)
            for subcommand in param_instance._command:
                if set(param).issubset(subcommand):
                    return
            assert False

        except ParameterError:
            pass

    def test_run_command_produces_expected_output(self, cd, task):
        for path in task.instance._output_paths:
            assert path.exists()

    def test_run_command_produces_zero_exit_code(self, cd, task):
        for code in task.instance.exit_code:
            assert code == 0

    def test_inplace_tasks_delete_their_inputs(self, cd, task):
        if task.instance._inplace:
            assert not any_exist(as_paths(task.instance.input_files))


def test_apply_symlink_produces_expected_output(cd, full_dir):
    symlink_all = apply_(SymlinkInputs)
    inputs = full_dir.glob('*')
    output_ids = ['tar.x', 'sar.x', 'lar.x']
    outputs = [full_dir / out for out in output_ids]
    symlink_all(as_strs(inputs), as_strs(outputs))
    for path in outputs:
        assert path.exists()


def test_format_():
    suffixes = [['fastq', 'fq', 'fastq.gz', 'fq.gz'], ['fa', 'fasta']]
    formatter = format_(suffixes)
    base_regex = ".*/(?P<PREFIX>[^\.]*)\.(?P<MID>.*)\.?(?P<SUFFIX>"

    actual_regex = formatter.args
    expected_regex = (
            ''.join([base_regex, 'fq.gz|fq|fastq.gz|fastq)']),
            ''.join([base_regex, 'fasta|fa)']))
    assert actual_regex == expected_regex

