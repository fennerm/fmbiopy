""" Classes and functions for use in Ruffus pipelines

Ruffus: http://www.ruffus.org.uk/
"""

from fmbiopy.fmtype import StringOrSequence
import logging
import os
from ruffus.proxy_logger import make_shared_logger_and_proxy
from ruffus.proxy_logger import setup_std_shared_logger

"""Default logging format"""
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"

class RuffusLog(object):
    """A logger/mutex pair used for logging in Ruffus pipelines

    Usage
    -----
        ruflog = RuffusLog("foo", "bar.log")
        with ruflog.mutex:
            ruflog.log.info("Logged message")

    Parameters
    ----------
    name
        The ID used to refer to the logging instance
    location
        The path to where the logfile should be stored

    Other Parameters
    ----------------
    formatter
        Specifies logging format
    delay
        If True, file opening deferred until first log action
    level
        Logging level, see logging module docs
    Attributes
    ----------
    config: Dict[str, str, bool, int]
        Arguments passed to logger setup
    log: proxy_logger.LoggerProxy
        The logging proxy instance
    mutex: multiprocessing.managers.AcquirerProxy
        The mutex lock used for synchronous access to logfile
    """

    def __init__(
            self,
            name: str,
            location: str,
            formatter: str = DEFAULT_LOG_FORMAT,
            delay: bool = True,
            level: int = logging.DEBUG) -> None:

        if not name:
            raise ValueError("Name cannot be empty")
        logdir = os.path.dirname(location)
        if logdir:
            if not os.path.exists(logdir):
                raise ValueError("Logfile directory does not exist")

        self.config = {
                'file_name' : location,
                'formatter' : formatter,
                'delay' : delay,
                'level' : level
                }
        self.log, self.mutex = make_shared_logger_and_proxy(
                setup_std_shared_logger, name, self.config)

    def write(self, message: StringOrSequence) -> None:
        """Log a message with mutex lock

        Parameters
        ----------
        message
            Message to be logged to the Ruffus logfile. List is converted to
            string and spaces are added. If already a string just log directly
        header, optional
            Header string. Will precede message and be wrapped in dashes to
            stand out from other text.

        """
        with self.mutex:
            if not isinstance(message, str):
                message = ' '.join(message)
            self.log.info(message)

    def _divider(self):
        """Write a text divider to logfile"""
        self.log.info('-' * 30)

    def write_header(self, message: StringOrSequence) -> None:
        self._divider()
        self.write(message)
        self._divider()

