"""Shared logging setup for every phase of the pipeline.

Historically the three phases each used a different mechanism (``tqdm.write`` +
``print`` in phase 1, the ``logging`` module in phase 2, bare ``print`` in
phase 3). :func:`get_logger` gives all of them one consistent, idempotent logger
that writes to the console and, optionally, to a per-run log file.
"""

import logging
import os
import sys

_CONSOLE_FMT = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_FILE_FMT = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")


def get_logger(name="cassa_photometry", run_dir=None, level=logging.INFO):
    """Return a configured logger.

    Parameters
    ----------
    name : str
        Logger name (shared across a phase so handlers are reused).
    run_dir : str, optional
        If given, a ``<name>.log`` file handler is attached inside this
        directory in addition to the console handler.
    level : int
        Logging level (default :data:`logging.INFO`).

    Notes
    -----
    The function is idempotent: calling it repeatedly with the same ``name``
    will not stack duplicate handlers.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
               for h in logger.handlers):
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(_CONSOLE_FMT)
        logger.addHandler(console)

    if run_dir:
        log_path = os.path.abspath(os.path.join(run_dir, f"{name}.log"))
        already = any(
            isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", None) == log_path
            for h in logger.handlers
        )
        if not already:
            os.makedirs(run_dir, exist_ok=True)
            file_handler = logging.FileHandler(log_path)
            file_handler.setFormatter(_FILE_FMT)
            logger.addHandler(file_handler)

    return logger
