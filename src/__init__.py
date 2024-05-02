"""Package docstring."""

import logging
from src import constants as C

def setup_logger(name):
    """Start a separate logger for each script."""

    formatter = logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s %(funcName)s() %(message)s", 
        datefmt='%d-%b-%y %H:%M:%S',
        )

    file_handler = logging.FileHandler(f"{C.LOGS_DIR}/{name}.log", "w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)        
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    if logger.hasHandlers(): logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger