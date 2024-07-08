"""Boilerplate code for most modules."""

import logging
from pathlib import Path

import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

logger = logging.getLogger(__name__)


def main():
    """Run as script."""
    pass


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()