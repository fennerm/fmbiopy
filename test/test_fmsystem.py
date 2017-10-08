import tempfile
import pytest
import fmbiopy.fmsystem as fmsystem
from fmbiopy.fmcheck import any_exist, all_exist


class TestRunCommand():

    def test_trivial_success(self):
        with tempfile.NamedTemporaryFile(mode='wt') as temp:
            command = "wc -m " + temp.name
            temp.write('xyz')
            temp.flush()
            succ = fmsystem.run_command(command)
            assert succ[0] == 0
            assert succ[1] == '3 '+ temp.name + '\n'
            assert succ[2] == ''

    def test_trivial_mixed_success(self):
        with tempfile.NamedTemporaryFile(mode='wt') as temp:
            command = "wc -m / " + temp.name
            temp.write('xyz')
            temp.flush()
            succ_and_fail = fmsystem.run_command(command)
            assert succ_and_fail[0] == 1
            assert succ_and_fail[1] == (
                    '      0 /\n      3 ' + temp.name + '\n      3 total\n')
            assert succ_and_fail[2] == 'wc: /: Is a directory\n'

    def test_trivial_failure(self):
        command = "wc -m /"
        fail = fmsystem.run_command(command)
        assert fail[0] == 1
        assert fail[1] == '0 /\n'
        assert fail[2] == 'wc: /: Is a directory\n'

    def test_shell_true_doesnt_fail(self):
        command = "echo foo | grep foo"
        process = fmsystem.run_command(command, shell=True)
        assert process[0] == 0
        assert process[1] == "foo\n"
        assert not process[2]


class TestDelete():
    def test_normal_usage(self):
        tmpfiles = [tempfile.NamedTemporaryFile() for i in range(3)]
        tmpfile_names = [tmpfile.name for tmpfile in tmpfiles]
        assert all_exist(tmpfile_names)
        with fmsystem.delete(tmpfile_names):
            pass

        assert not any_exist(tmpfile_names)

class TestParseParamDict():
    @pytest.fixture
    def example_dict(self):
        example = {
                '-a':'a_val',
                '-b':'b_val',
                '--long':'long_val'
                }
        return example

    def test_parse_correct_output(self, example_dict):
        parsed = fmsystem.parse_param_dict(example_dict)
        assert '-a a_val' in parsed
        assert '-b b_val' in parsed
        assert '--long long_val' in parsed

    def test_empty_input(self):
        assert fmsystem.parse_param_dict({}) == ''
