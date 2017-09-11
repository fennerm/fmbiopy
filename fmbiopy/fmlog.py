""" Logging utilities """

import logging


def setup_log(logfile: str) -> None:
    """Set up a logging instance for a given logfile path. See:
    https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing\
    -to-a-file-at-the-same-time """
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s \
                                %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=logfile,
                        filemode='w')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')

    # tell the handler to use this format
    console.setFormatter(formatter)

    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
