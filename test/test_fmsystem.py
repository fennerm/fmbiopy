## Pytests for run_bowtie2.py

import os
import tempfile
import pytest
from fmbiopy import fmsystem
from fmbiopy.fmcheck import any_exist, all_exist

class TestRunCommand():

    def test_trivial_success(self):
        with tempfile.NamedTemporaryFile(mode='wt') as temp:
            command = "wc -m " + temp.name
            temp.write('xyz')
            temp.flush()
            succ = fmsystem.run_command(command)
            assert succ[0] == 0
            assert succ[1] == '3 '+temp.name + '\n'
            assert succ[2] == ''

    def test_trivial_mixed_success(self):
        with tempfile.NamedTemporaryFile(mode='wt') as temp:
            command = "wc -m / " + temp.name
            temp.write('xyz')
            temp.flush()
            succ_and_fail = fmsystem.run_command(command)
            assert succ_and_fail[0] == 1
            assert succ_and_fail[1] == (
                    '      0 /\n      3 '+temp.name+'\n      3 total\n')
            assert succ_and_fail[2] == 'wc: /: Is a directory\n'

    def test_trivial_failure(self):
        command = "wc -m /"
        fail = fmsystem.run_command(command)
        assert fail[0] == 1
        assert fail[1] == '0 /\n'
        assert fail[2] == 'wc: /: Is a directory\n'

class TestDelete():
    def test_normal_usage(self):
        tmpfiles = [tempfile.NamedTemporaryFile() for i in range(3)]
        tmpfile_names = [tmpfile.name for tmpfile in tmpfiles]
        assert all_exist(tmpfile_names)
        with fmsystem.delete(tmpfile_names):
            raise OSError

        assert not any_exist(tmpfile_names)
