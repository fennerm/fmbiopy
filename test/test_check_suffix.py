import pytest
from fmbiopy.fmcheck import check_all_suffix

def test_all_correct():
    check_all_suffix(['a.x', 'b.y.x'], ['x'])

def test_some_correct():
    with pytest.raises(Exception):
        check_all_suffix(['a.y', 'b.y.x'], ['y'])

def test_none_correct():
    with pytest.raises(Exception):
        check_all_suffix(['a.x', 'b.y.x'], ['y'])
