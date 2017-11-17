from pytest import fixture

from fmbiopy.fmparse import *

@fixture(name='argv')
def gen_argv():
    return ['prog', '--flag1', '-1', '-x', '--a1', 'b1', '--a2=b2', '-a3',
            'b3', 'b4', 'b5']

def test_simple_parse(argv):
    expected = {
            '--flag1': None,
            '-x': None,
            '--a1': 'b1',
            '--a2': 'b2',
            '-a3': 'b3',
            'args' : ['b4', 'b5']
            }
    assert simple_parse(argv, aliases = [('--alias', '-1')]) == expected

