""""Package docstring."""

import logging
from pathlib import Path


def setup_logger(logfile=None, name="src", level=logging.DEBUG, stream=False):
    """Set up a logger."""

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        "{asctime} {levelname} | {module}  {funcName}  {lineno} :: {message}\n",
        style="{",
        datefmt="%d-%m-%y %H:%M:%S",
    )

    if logfile:
        file_handler = logging.FileHandler(logfile, mode="w")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    if stream:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    return logger


def log_file(file):
    return f"data/logs/{'.'.join(Path(file).relative_to(Path.cwd()).with_suffix('.log').parts)}"


root_logger = setup_logger(name="")  # Sponge for logging output from external packages
src_logger = setup_logger(stream=True)  # Main logger for the package
