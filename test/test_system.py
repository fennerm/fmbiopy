from plumbum import local
from plumbum.cmd import echo

from fmbiopy.system import *


def test_capture_stdout_single_item():
    assert capture_stdout(echo["foo"]) == ["foo"]


def test_capture_stdout_no_output():
    assert capture_stdout(echo) == []


def test_capture_stdout_only_stderr():
    gen_stderr = local["test/bin/gen_stderr"]
    assert capture_stdout(gen_stderr["foo"]) == []


def test_capture_stdout_multiple_items():
    assert capture_stdout(echo["foo\nbar\n"]) == ["foo", "bar"]


def test_capture_stdout_strips_whitespace():
    assert capture_stdout(echo["\nfoo\n\n"]) == ["foo"]
