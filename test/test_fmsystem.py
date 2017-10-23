"""Test suite for fmbiopy.fmsystem"""
from pathlib import Path
import tempfile

import pytest


from fmbiopy.fmpaths import (
        any_exist,
        all_exist,
        )
from fmbiopy.fmsystem import *

class TestRunCommand():

    def test_trivial_success(self, tmpdir, gen_tmp):
        pytest.set_trace()
        tmpfile = gen_tmp(directory=tmpdir)
        with tmpfile.open(mode='wt') as temp:
            command = ['wc', '-m', str(tmpfile)]
            temp.write('xyz')
            temp.flush()
            succ = run_command(command)
            assert succ[0] == 0
            assert succ[1] == '3 '+ str(tmpfile) + '\n'
            assert succ[2] == ''

    def test_trivial_mixed_success(self):
        with tempfile.NamedTemporaryFile(mode='wt') as temp:
            command = ['wc', '-m', '/', temp.name]
            temp.write('xyz')
            temp.flush()
            succ_and_fail = run_command(command)
            assert succ_and_fail[0] == 1
            assert succ_and_fail[1] == (
                    '      0 /\n      3 ' + temp.name + '\n      3 total\n')
            assert succ_and_fail[2] == 'wc: /: Is a directory\n'

    def test_trivial_failure(self):
        command = ['wc', '-m', '/']
        fail = run_command(command)
        assert fail[0] == 1
        assert fail[1] == '0 /\n'
        assert fail[2] == 'wc: /: Is a directory\n'

    def test_shell_true_doesnt_fail(self):
        command = ['echo', 'foo', '|', 'grep', 'foo']
        process = run_command(command, shell=True)
        assert process[0] == 0
        assert process[1] == "foo\n"
        assert not process[2]


class TestDelete():
    def test_normal_usage(self, gen_tmp):
        tmpfiles = [gen_tmp() for i in range(3)]
        assert all_exist(tmpfiles)
        with delete(tmpfiles):
            pass

        assert not any_exist(tmpfiles)


class TestParseParamDict():
    @pytest.fixture
    def example_dict(self):
        example = {
                '-a': 'a_val',
                '-b': 'b_val',
                '--long': 'long_val'}
        return example

    def test_parse_correct_output(self, example_dict):
        parsed = parse_param_dict(example_dict)
        assert '-a a_val' in parsed
        assert '-b b_val' in parsed
        assert '--long long_val' in parsed

    def test_empty_input(self):
        assert parse_param_dict({}) == ''

class TestWorkingDirectory():
    def test_changes_directory(self, tmpdir):
        with working_directory(tmpdir):
            assert Path.cwd() == tmpdir


