"""Test fmbiopy.fmpaths"""

from pytest import (
    fixture,
    mark,
    raises,
)

from plumbum import (
    local,
    LocalPath,
)

from fmbiopy.fmerr import EmptyListError
from fmbiopy.fmpaths import *

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def test_all_exist(poss_path_lists):
    (name, value) = poss_path_lists
    if name == 'absolute_exist_paths':
        assert all_exist(value)
    elif name == 'relative_exist_paths':
        assert all_exist(value)
    elif name == 'absolute_some_exist_paths':
        assert not all_exist(value)
    elif name == 'absolute_nonexist_paths':
        assert not all_exist(value)
    elif name == 'empty_list':
        with raises(EmptyListError):
            all_exist(value)


def test_any_dont_exist(poss_path_lists):
    (name, value) = poss_path_lists
    if name == 'absolute_exist_paths':
        assert not any_dont_exist(value)
    elif name == 'relative_exist_paths':
        assert not any_dont_exist(value)
    elif name == 'absolute_some_exist_paths':
        assert any_dont_exist(value)
    elif name == 'absolute_nonexist_paths':
        assert any_dont_exist(value)
    elif name == 'relative_nonexist_paths':
        assert any_dont_exist(value)
    elif name == 'empty_list':
        with raises(EmptyListError):
            any_dont_exist(value)


def test_any_exist(poss_path_lists):
    (name, value) = poss_path_lists
    if name == 'absolute_exist_paths':
        assert any_exist(value)
    elif name == 'relative_exist_paths':
        assert any_exist(value)
    elif name == 'absolute_some_exist_paths':
        assert any_exist(value)
    elif name == 'absolute_nonexist_paths':
        assert not any_exist(value)
    elif name == 'relative_nonexist_paths':
        assert not any_exist(value)
    elif name == 'empty_list':
        with raises(EmptyListError):
            any_exist(value)


def test_all_empty(poss_path_lists):
    (name, value) = poss_path_lists
    if name == 'absolute_exist_paths':
        assert all_empty(value)
    elif name == 'absolute_nonexist_paths':
        with raises(FileNotFoundError):
            all_empty(value)
    elif name == 'empty_list':
        with raises(EmptyListError):
            all_empty(value)
    elif name == 'nonempty_paths':
        assert not all_empty(value)


def test_apply_is_empty():
    # Already covered
    pass


def test_apply_exists():
    # Already covered
    pass


def test_as_dict(nested_dir):
    dic = as_dict(nested_dir)
    assert isinstance(dic, dict)


def test_as_paths(randstrs):
    paths = as_paths(randstrs(3))
    for path in paths:
        assert isinstance(path, LocalPath)


def test_as_strs(empty_list):
    # Trival
    pass


def test_check_all_exist():
    # Already covered
    pass


def test_create_all(poss_path_lists):
    name, value = poss_path_lists
    if name == 'absolute_nonexist_paths':
        create_all(value)
        assert all_exist(value)
    elif name == 'empty_list':
        create_all(value)


class TestDelete():
    def test_normal_usage(self, gen_tmp):
        tmpfiles = [gen_tmp() for i in range(3)]
        assert all_exist(tmpfiles)
        with delete(tmpfiles):
            pass

        assert not any_exist(tmpfiles)


@mark.parametrize("ext,substr,expect",
                  [(None, None, ['a.x', 'b.y', 'b.x']), ('y', None, ['b.y']),
                   (None, 'a', ['a.x']), (None, 'k', []), ('x', 'b', ['b.x'])])
def test_find(ext, substr, expect, full_dir):
    if expect:
        expect = sorted([full_dir / exp for exp in expect])
    assert find(full_dir, extensions=ext, substring=substr) == expect


def test_get_bowtie2_indices(bowtie2_suffixes):
    prefix = 'foo'
    expected_paths = [prefix + suffix for suffix in bowtie2_suffixes]
    dot_paths = get_bowtie2_indices(prefix)
    actual_paths = [str(path) for path in dot_paths]
    assert actual_paths == expected_paths


def test_get_bowtie2_indices():
    # Trivial
    pass


def test_is_empty(empty_path):
    assert is_empty(empty_path)
    with empty_path.open('w') as f:
        f.write('0')
    assert not is_empty(empty_path)


def test_listdirs(gen_tmp, nested_dir):
    reg_file = gen_tmp(empty=False, directory=nested_dir)
    assert str(reg_file) not in as_strs(listdirs(nested_dir))


def test_move(gen_tmp, tmpdir, randpath):
    src = gen_tmp(empty=False, directory=tmpdir)
    link1, link2, dest = [randpath() for i in range(3)]
    src.symlink(link1)
    link1.symlink(link2)
    move(link2, dest)
    assert dest.exists()
    assert not link2.exists()
    assert src.exists()


def test_prefix(double_suffixed_path):
    pre = prefix(double_suffixed_path)
    assert '/' not in pre
    assert '.' not in pre


def test_recursive_list_contents(nested_dir):
    expected = [
        "foo", "bar", "car", "foo/a.x", "foo/b.y", "bar/a.x", "bar/b.y",
        "car/a.x", "car/b.y"
    ]
    expected = [nested_dir / x for x in expected]
    assert set(recursive_list_contents(nested_dir)) == set(expected)


def test_rm_gz_suffix(double_suffixed_path, gzipped_path):
    for path in [double_suffixed_path, gzipped_path]:
        nsuffixes = len(path.suffixes)
        nsuffixes_after = len(rm_gz_suffix(path).suffixes)
        diff = nsuffixes - nsuffixes_after
        if path.suffix == '.gz':
            expected_diff = 2
        else:
            expected_diff = 1
        assert diff == expected_diff


def test_root(double_suffixed_path):
    actual = str(root(double_suffixed_path))
    par = double_suffixed_path.up()
