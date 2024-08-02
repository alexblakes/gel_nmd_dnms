"""Boilerplate code for most modules."""

import logging

import pandas as pd

import src

logger = logging.getLogger(__name__)


def main():
    """Run as script."""
    return None


if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
