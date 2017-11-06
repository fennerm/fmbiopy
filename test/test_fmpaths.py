"""Test fmbiopy.fmpaths"""

from pathlib import Path
from pytest import (
        fixture,
        mark,
        raises,
        )


from fmbiopy.fmerr import EmptyListError
from fmbiopy.fmpaths import *
from fmbiopy.fmsystem import working_directory


def test_absglob(full_dir):
    glb = absglob(full_dir, '*')
    assert all_absolute(glb)


def test_add_suffix(randsuffix, randpath):
    suff = randsuffix()
    suffixed = add_suffix(randpath(), suff)
    assert suffixed.suffix == suff
    assert suffixed.is_absolute()


def test_all_absolute(poss_path_lists):
    (name, value) = poss_path_lists
    if name == 'absolute_exist_paths':
        assert all_absolute(value)
    elif name == 'relative_nonexist_paths':
        assert not all_absolute(value)
    elif name == 'mixed_absolute_relative_paths':
        assert not all_absolute(value)
    elif name == 'empty_list':
        with raises(EmptyListError):
            all_absolute(value)


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
        assert isinstance(path, Path)


def test_as_strs(empty_list):
    # Trival
    pass


def test_check_all_exist():
    # Already covered
    pass


def test_create_all(poss_path_lists):
    name, value = poss_path_lists
    if name == 'absolute_exist_paths':
        with raises(FileExistsError):
            create_all(value)
    elif name == 'absolute_nonexist_paths':
        create_all(value)
        assert all_exist(value)
    elif name == 'empty_list':
        create_all(value)
    elif name == 'nonempty_paths':
        with raises(FileExistsError):
            create_all(value)


def test_extension_without_gz(gzipped_path, double_suffixed_path):
    pass

@mark.parametrize("ext,substr,expect", [
        (None, None, [Path('a.x'), Path('b.y'), Path('b.x')]),
        ('y', None, [Path('b.y')]),
        (None, 'a', [Path('a.x')]),
        (None, 'k', []),
        ('x', 'b', [Path('b.x')])])
def test_find(ext, substr, expect, full_dir):
    if expect:
        expect = sorted([full_dir.joinpath(exp) for exp in expect])
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
    link1.symlink_to(src)
    link2.symlink_to(link1)
    move(link2, dest)
    assert dest.exists()
    assert not link2.exists()
    assert src.exists()


def test_prefix(double_suffixed_path):
    pre = prefix(double_suffixed_path)
    assert '/' not in pre
    assert '.' not in pre


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
    par = double_suffixed_path.parent
    expected = str(par / Path(double_suffixed_path.stem).stem)
    assert actual == expected
