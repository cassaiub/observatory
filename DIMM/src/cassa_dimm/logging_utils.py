"""Shared logging setup for the DIMM monitor (console + optional file)."""

import logging
import os
import sys

_CONSOLE_FMT = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_FILE_FMT = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")


def get_logger(name="cassa_dimm", log_dir=None, level=logging.INFO):
    """Return a configured, idempotent logger.

    Parameters
    ----------
    name : str
        Logger name (handlers are reused, not stacked).
    log_dir : str, optional
        If given, also write to ``<log_dir>/<name>.log``.
    level : int
        Logging level.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
               for h in logger.handlers):
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(_CONSOLE_FMT)
        logger.addHandler(console)

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.abspath(os.path.join(log_dir, f"{name}.log"))
        if not any(isinstance(h, logging.FileHandler)
                   and getattr(h, "baseFilename", None) == log_path for h in logger.handlers):
            fh = logging.FileHandler(log_path)
            fh.setFormatter(_FILE_FMT)
            logger.addHandler(fh)

    return logger
