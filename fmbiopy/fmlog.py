#'''Logging utilities'''
import logging
import os


def setup_log():
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
