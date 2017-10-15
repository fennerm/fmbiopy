"""Test suite for fmbiopy.fmclass.py"""
import os
from pathlib import Path
import pytest

import fmbiopy.fmclass as fmclass
import fmbiopy.fmpaths as fmpaths
from fmbiopy.fmtest import gen_tmp
from fmbiopy.fmtest import load_sandbox


@pytest.fixture(scope='module')
def minimal_classfile(load_sandbox):
    tmpfile = gen_tmp(empty=True, suffix='.py')

    eg = "class Foo(object):\n   def __init__(self):\n \
        super().__init__(self)\n\n"
    eg += "class Bar(Foo):\n   def __init__(self):\n       pass\n\n"
    eg += "class Car(Exception):\n def __init__(self):\
            \n      super().__init__()\n\n"

    with tmpfile.open('w') as f:
        f.write(eg)
    class_file = Path('.'.join(['test', 'sandbox', tmpfile.name]))
    class_name = class_file.with_suffix('').name
    return class_name


class TestListClasses(object):
    def test_simple_case(self, minimal_classfile):
        actual = fmclass.list_classes(minimal_classfile)
        actual = [fmclass.classname(cls) for cls in actual]
        expected = ['Bar', 'Car', 'Foo']
        assert actual == expected

    def test_exclude(self, minimal_classfile):
        actual = fmclass.list_classes(
                minimal_classfile,
                exclude='Foo')
        actual = [fmclass.classname(cls) for cls in actual]
        expected = ['Bar', 'Car']
        assert actual == expected

    def test_choose_type(self, minimal_classfile):
        actual = fmclass.list_classes(
                minimal_classfile,
                of_type='Foo')
        actual = [fmclass.classname(cls) for cls in actual]
        expected = ['Bar', 'Foo']
        assert actual == expected
