import os, pytest
from fmbiopy.fmcheck import any_dont_exist

def test_all_exist():
    x = os.listdir(".")
    assert not any_dont_exist(x)

def test_some_exist():
    x = os.listdir(".") + ["doesnt_exist.txt"]
    assert any_dont_exist(x)

def test_none_exist():
    x = ["x.xy", "y.yz"]
    assert any_dont_exist(x)
