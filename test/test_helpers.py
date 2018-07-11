from pytest import fixture, mark
from pytest_lazyfixture import lazy_fixture

from test.helpers import *


@fixture
def meta():
    # This is a hack which is required for file_generator to work.
    # file_generator always overwrites the meta variable in wrapped functions.
    # However, if fixture named meta is not defined, pytest will throw an error
    # before file_generator gets a chance to run.
    pass


@fixture
def dummy_fixture():
    """A blank/meaningless fixture.

    Used to test file_generator.
    """
    return "dummy"


@fixture
@file_generator
def dummy_func_with_fixture_dependencies(meta, dummy_fixture):
    """A function which is wrapped by both `file_generator` and `fixture`."""
    # Check that meta was passed intact
    meta["file"]


@fixture
@file_generator
def dummy_func_without_fixture_dependencies(meta):
    pass


@mark.parametrize(
    "dummy_func",
    [
        lazy_fixture("dummy_func_with_fixture_dependencies"),
        lazy_fixture("dummy_func_without_fixture_dependencies"),
    ],
)
def test_that_file_generator_injects_return(dummy_func):
    assert "file" in dummy_func
