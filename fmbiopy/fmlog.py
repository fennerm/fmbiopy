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
from typing import Tuple

from plumbum import (
        local,
        LocalPath,
        )


def setup_log(logfile: LocalPath) -> None:
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
