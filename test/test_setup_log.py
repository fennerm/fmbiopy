import logging, pytest, tempfile
from tempfile import NamedTemporaryFile
import fmbiopy
from fmbiopy.fmlog import setup_log

def gen_tmp_logfile():
    tmp = NamedTemporaryFile()
    setup_log(tmp.name)
    return tmp

def test_root_logging():
    """ Test that logging.info writes a line to the logfile""" 

    logfile = gen_tmp_logfile()
    logging.info("foo")
    line_count = sum(1 for line in open(logfile.name))
    assert line_count == 1
