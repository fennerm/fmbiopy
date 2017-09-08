import fmbiopy.fmruffus as fmruffus
import os
import pytest
from tempfile import NamedTemporaryFile as tmp


class TestRuffusLog(object):
    def test_empty_name_returns_value_error(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("", tmp().name)

    def test_non_existant_path(self):
        with pytest.raises(ValueError):
            fmruffus.RuffusLog("foo", "bar/bar.log")

    def test_normal_usage(self):
        temp = tmp()
        ruflog = fmruffus.RuffusLog("foo", temp.name)
        with ruflog.mutex:
            ruflog.log.info("Test")
        assert os.path.getsize(temp.name) > 0
