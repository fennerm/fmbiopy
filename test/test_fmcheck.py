"""Test module for fmbiopy.fmcheck"""

import os
import pytest

import fmbiopy.fmcheck as fmcheck


class TestCheckAllSuffix():
    def test_all_correct(self):
        fmcheck.check_all_suffix(['a.x', 'b.y.x'], ['x'])

    def test_some_correct(self):
        with pytest.raises(Exception):
            fmcheck.check_all_suffix(['a.y', 'b.y.x'], ['y'])

    def test_none_correct(self):
        with pytest.raises(Exception):
            fmcheck.check_all_suffix(['a.x', 'b.y.x'], ['y'])
