""" Logging utilities """

from logging import (
        basicConfig,
        DEBUG,
        Formatter,
        getLogger,
        INFO,
        StreamHandler,
        )
from multiprocessing.managers import AcquirerProxy  #type: ignore
from pathlib import Path
from typing import Tuple

from ruffus.proxy_logger import (
        LoggerProxy,
        make_shared_logger_and_proxy,
        setup_std_shared_logger,
        )

"""Default logging format"""
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"

class MutexLogger(object):
    """A logger/mutex pair used for logging in Duffus pipelines

    Examples
    --------
    >>>>ruflog = MutexLogger("foo", "bar.log")
    >>>>with ruflog.mutex:
            ruflog.log.info("logged message")

    Attributes
    ----------
    config: Dict[str, str, bool, int]
        Arguments passed to logger setup
    log: LoggerProxy
        The logging proxy instance
    mutex: AcquirerProxy
        The mutex lock used for synchronous access to logfile
    """
    def __init__(
            self,
            name: str,
            location: Path,
            formatter: str = DEFAULT_LOG_FORMAT,
            delay: bool = True,
            level: int = DEBUG,
            buffered: bool = True,
            overwrite: bool = True) -> None:
        """Initialization

        Parameters
        ----------
        name
            The ID of the logger
        location
            The path to where the logfile should be stored
        formatter
            Specifies logging format
        delay
            If True, file opening deferred until first log action
        level
            Logging level, see logging module docs
        buffered
            If True, the `write` statements are appended to the buffer, but
            nothing is logged until the buffer is `flush`ed.
        overwrite
            If True, logfile will be overwritten if it already exists

        Raises
        ------
        FileNotFoundError
            Parent directory of logfile location doesn't exist
        """


        self.config = {
                'file_name': str(location),
                'formatter': formatter,
                'delay': delay,
                'level': level,
                'buffered': buffered}

        if self.config['buffered']:
            self._buffer: str = ''

        if not location.parent.exists():
            raise FileNotFoundError('Parent directory of logfile doesnt exist')

        # Overwrite logfile if it already exists
        if overwrite:
            if location.exists():
                location.unlink()

        logprox: Tuple[LoggerProxy, AcquirerProxy] = \
                make_shared_logger_and_proxy(
                        setup_std_shared_logger, name, self.config)
        self.log, self.mutex = logprox

    def flush(self)-> None:
        """Write the buffer to the logfile and clear it"""
        self.log.info(self._buffer)
        self._clear_buffer()

    def write(self, msg: str)-> None:
        """Log a message with mutex lock

        Parameters
        ----------
        msg
            Message to be logged to the Duffus logfile. List is converted to
            string and spaces are added. If already a string just log directly
        """
        with self.mutex:
            if self.config['buffered']:
                self._new_msg(msg)
            else:
                self.log.info(msg)

    def write_header(self, msg: str, subheader: bool = False)-> None:
        """Write a header to the logfile

        Parameters
        ----------
        msg
            A string or list of words to be written as the header text
        sub, optional
            If True, write a subheader, instead of a header
        """
        self._write_divider(subheader)
        self.write(msg)
        self._write_divider(subheader)

    def _clear_buffer(self)-> None:
        self._buffer = ''

    def _new_msg(self, msg: str)-> None:
        """Add a new message to the buffer"""
        self._buffer += ''.join(['\n', msg])

    def _write_divider(self, subheader: bool)-> None:
        """Write a text divider to logfile

        Parameters
        ----------
        subheader
            If True, the divider is made up of dashes. Else, the divider is
            made up of equal signs
        """
        if subheader:
            self.write('-' * 80)
        else:
            self.write('=' * 80)


def setup_log(logfile: Path) -> None:
    """Set up a logging instance for a given logfile path. See:
    https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing\
    -to-a-file-at-the-same-time """
    basicConfig(level=DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s \
                                %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=str(logfile),
                        filemode='w')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = StreamHandler()
    console.setLevel(INFO)

    # set a format which is simpler for console use
    formatter = Formatter('%(name)-12s: %(levelname)-8s %(message)s')

    # tell the handler to use this format
    console.setFormatter(formatter)

    # add the handler to the root logger
    getLogger('').addHandler(console)
