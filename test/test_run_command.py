## Pytests for run_bowtie2.py

import os
import tempfile
import pytest
from fmbiopy import fmsystem

def triv_success(f):
    x = "wc -m " + f
    return x

def triv_fail():
    x = "wc -m /"
    return x
def triv_succ_and_fail(f):
    x = "wc -m / " + f
    return x

def test_trivial():
    with tempfile.NamedTemporaryFile(mode='wt') as temp:
        temp.write('xyz')
        temp.flush()
        succ = fmsystem.run_command(triv_success(temp.name))
        assert succ[0] == 0
        assert succ[1] == '3 '+temp.name + '\n'
        assert succ[2] == ''

    with tempfile.NamedTemporaryFile(mode='wt') as temp:
        temp.write('xyz')
        temp.flush()
        succ_and_fail = fmsystem.run_command(triv_succ_and_fail(temp.name))
        assert succ_and_fail[0] == 1
        assert succ_and_fail[1] == '      0 /\n      3 '+temp.name+'\n      3 total\n'
        assert succ_and_fail[2] == 'wc: /: Is a directory\n'

    fail = fmsystem.run_command(triv_fail())
    assert fail[0] == 1
    assert fail[1] == '0 /\n'
    assert fail[2] == 'wc: /: Is a directory\n'

