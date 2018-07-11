"""Test suite for fmbiopy.py"""
import sys

from plumbum import local
from pytest import fixture, mark
from pytest_lazyfixture import lazy_fixture

from fmbiopy.obj import *


@fixture(scope="class", name="class_text")
def gen_class_text():
    txt = "class Foo(object):\n    def __init__(self):\n\
            super().__init__(self)\n\n"
    txt += "class Bar(Foo):\n    def __init__(self):\n        pass\n\n"
    txt += "class Car(Exception):\n    def __init__(self):\n\
            super().__init__()\n\n"
    return txt


@fixture(scope="class")
def egclassfile(gen_tmp, testdir, sandbox, class_text):
    # We have to give each classfile a unique directory to prevent interference
    # from previous imports in the namespace
    tmpfile = gen_tmp(empty=True, directory=sandbox, suffix=".py")
    with tmpfile.open("w") as f:
        f.write(class_text)
    module_name = tmpfile.with_suffix("").name
    package_name = ".".join([testdir.name, sandbox.name])
    return (module_name, package_name)


@fixture(name="local_classfile")
def gen_local_classfile(gen_tmp, class_text):
    tmpfile = gen_tmp(empty=True, directory=local.cwd, suffix=".py")
    with tmpfile.open("w") as f:
        f.write(class_text)
    module_name = tmpfile.with_suffix("").name
    yield module_name
    tmpfile.delete()
    tmpfile_pyc = local.path(tmpfile + "c")
    if tmpfile_pyc.exists():
        tmpfile_pyc.delete()


@fixture
def egclass():
    class A(object):
        pass

    class B(object):
        pass

    class C(A, B):
        pass

    return C


def classnames(classes):
    return [classname(cls) for cls in classes]


@fixture(scope="class")
def all_classes(egclassfile):
    classes = list_classes(module=egclassfile[0], package=egclassfile[1])
    return classnames(classes)


@fixture(scope="class")
def no_foo(egclassfile):
    classes = list_classes(
        module=egclassfile[0], package=egclassfile[1], exclude=["Foo"]
    )
    return classnames(classes)


@fixture(scope="class")
def foo_type(egclassfile):
    classes = list_classes(
        module=egclassfile[0], package=egclassfile[1], of_type=["Foo"]
    )
    return classnames(classes)


class TestListClasses(object):
    def test_simple_case(self, all_classes):
        expected = ["Bar", "Car", "Foo"]
        assert all_classes == expected

    def test_exclude(self, no_foo):
        expected = ["Bar", "Car"]
        assert no_foo == expected

    def test_choose_type(self, foo_type):
        expected = ["Bar", "Foo"]
        assert foo_type == expected

    def test_module_only(self, local_classfile, all_classes):
        expected = ["Bar", "Car", "Foo"]
        sys.path.append(str(local.cwd))
        assert classnames(list_classes(local_classfile)) == expected


def test_parent_names(egclass):
    assert parent_names(egclass) == ["A", "B"]


@mark.parametrize(
    "func,params", [(lambda a, b: True, ["a", "b"]), (lambda: True, [])]
)
def test_get_param_names(func, params):
    assert get_param_names(func) == params


@mark.parametrize("func", [lambda a, b: True, lambda: True])
def test_replace_param_sig(func):
    def model_func(a, b):
        pass

    modified_func = replace_param_sig(func, ["a", "b"])
    assert signature(modified_func) == signature(model_func)
