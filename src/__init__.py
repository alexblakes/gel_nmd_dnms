""""Package docstring."""

import logging


_FORMATTER = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(filename)s %(funcName)s(): %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)


def setup_logger(
    logfile=None,
    mode="w",
    name="src",
    formatter=_FORMATTER,
    stream=False,
    level=logging.DEBUG,
):
    """Retrieve a new or existing logger and optionally add handler(s) to it."""

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Log to file
    if logfile:
        file_handler = logging.FileHandler(logfile, mode=mode)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(level)
        logger.addHandler(file_handler)

    # Log to console
    if stream:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        stream_handler.setLevel(level)
        logger.addHandler(stream_handler)

    return logger


root_logger = setup_logger(name="") # Sponge for logging output from external packages
src_logger = setup_logger(stream=True) # Main logger for the package
