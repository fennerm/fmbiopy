import os, pytest
import fmbiopy.fmcheck as fmcheck

class TestAnyDontExist():
    def test_all_exist(self):
        x = os.listdir(".")
        assert not fmcheck.any_dont_exist(x)

    def test_some_exist(self):
        x = os.listdir(".") + ["doesnt_exist.txt"]
        assert fmcheck.any_dont_exist(x)

    def test_none_exist(self):
        x = ["x.xy", "y.yz"]
        assert fmcheck.any_dont_exist(x)

class TestCheckAllSuffix():
    def test_all_correct(self):
        fmcheck.check_all_suffix(['a.x', 'b.y.x'], ['x'])

    def test_some_correct(self):
        with pytest.raises(Exception):
            fmcheck.check_all_suffix(['a.y', 'b.y.x'], ['y'])

    def test_none_correct(self):
        with pytest.raises(Exception):
            fmcheck.check_all_suffix(['a.x', 'b.y.x'], ['y'])
