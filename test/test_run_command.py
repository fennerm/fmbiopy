## Pytests for run_bowtie2.py

import os, tempfile, pytest
import fmbiopy.fen_util as fen_util

def triv_success(f):
    return "wc -m " + f
def triv_fail():
    return "wc -m /"
def triv_succ_and_fail(f):
    return "wc -m / " + f

def test_trivial():
    with tempfile.NamedTemporaryFile() as temp:
        temp.write('xyz')
        temp.flush()
        succ = fen_util.run_command(triv_success(temp.name))
        assert succ[0] == 0
        assert succ[1] == '3 '+temp.name + '\n'
        assert succ[2] == ''

    with tempfile.NamedTemporaryFile() as temp:
        temp.write('xyz')
        temp.flush()
        succ_and_fail = fen_util.run_command(triv_succ_and_fail(temp.name))
        assert succ_and_fail[0] == 1
        assert succ_and_fail[1] == '      0 /\n      3 '+temp.name+'\n      3 total\n'
        assert succ_and_fail[2] == 'wc: /: Is a directory\n'
            
    fail = fen_util.run_command(triv_fail())
    assert fail[0] == 1
    assert fail[1] == '0 /\n'
    assert fail[2] == 'wc: /: Is a directory\n'

