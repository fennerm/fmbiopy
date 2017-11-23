"""Test suite for fmbiopy.fmsystem"""
from pathlib import Path
import tempfile

from pytest import (
        fixture,
        mark,
        )
from factory import (
        Factory,
        lazy_attribute,
        )
from pytest_factoryboy import (
        LazyFixture,
        register,
        )

from fmbiopy.fmlog import setup_log
from fmbiopy.fmpaths import (
        any_exist,
        all_exist,
        )
from fmbiopy.fmsystem import *

# class Command():
#     def __init__(self, **kwargs):
#         for key, value in kwargs.items():
#             setattr(self, key, value)
#
#         if self.logfile is not None:
#             setup_log(self.logfile)
#
#
# class CommandFactory(Factory):
#     success = True
#     as_list = True
#     shell = True
#     logfile = None
#
#     class Meta:
#         model = Command
#
#     @lazy_attribute
#     def command(self):
#         if self.success:
#             command = ['echo', 'foo']
#         else:
#             command = ['false']
#
#         if self.as_list:
#             command = ' '.join(command)
#         return command
#
# register(CommandFactory)

# class TestBash():
#     @mark.parametrize('command__as_list', [True, False])
#     def test_successful_command(self, command):
#         assert False # Rewrite logging in bash function
#         results = bash(command.command, shell=command.shell)
#         assert results == (0, 'foo\n', '')
#
#     @mark.parametrize('command__as_list', [True, False])
#     @mark.parametrize('command__success', [False])
#     @mark.parametrize('command__shell', [True, False])
#     def test_failed_command(self, command):
#         results = bash(command.command, shell=command.shell)
#         assert results == (1, '', '')
#
#     @mark.parametrize('command__logfile', [LazyFixture(
#         lambda randpath: randpath())])
#     @mark.parametrize('command__command', ['echo foo 1>&2'])
#     @mark.parametrize('command__shell', [True])
#     def test_logging(self, command, randpath):
#         bash(command.command, shell = command.shell)
#         assert 'foo' in command.logfile.read_text()



class TestParseParamDict():
    @fixture
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
