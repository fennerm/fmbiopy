"""Test suite for fmbiopy.fmlist"""
import pytest

import fmbiopy.fmlist as fmlist

class TestSplitList(object):
    @pytest.fixture
    def normal_case(self):
        return [1, 2, 0, 3, 0, 4]

    @pytest.fixture
    def split_at_end_case(self):
        return [1, 2, 0]

    def test_normal_usecase(self, normal_case):
        assert fmlist.split_list(normal_case, 0) == [[1, 2], [3], [4]]

    def test_empty_input(self):
        assert fmlist.split_list([], 0) == [[]]

    def test_no_split_found(self, normal_case):
        assert fmlist.split_list(normal_case, 20) == [normal_case]

    def test_nonlist_raises_type_error(self):
        with pytest.raises(TypeError):
            fmlist.split_list(1, 1)
