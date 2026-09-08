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

#: Level applied when a caller does not name one. ``load_config`` sets this from
#: ``PipelineConfig.log_level``, which is how a YAML file's ``log_level: DEBUG``
#: reaches the many ``get_logger`` calls that pass no level of their own.
_DEFAULT_LEVEL = logging.INFO


def set_default_level(level):
    """Set the level used by :func:`get_logger` when none is given.

    Accepts a level name (``"DEBUG"``) or a numeric level. An unrecognised name
    is ignored with a warning rather than raised: a typo in a log setting must
    not stop a reduction.
    """
    global _DEFAULT_LEVEL
    if isinstance(level, str):
        resolved = logging.getLevelName(level.strip().upper())
        if not isinstance(resolved, int):
            logging.getLogger(__name__).warning(
                "Unknown log_level %r; keeping %s.", level, logging.getLevelName(_DEFAULT_LEVEL)
            )
            return _DEFAULT_LEVEL
        level = resolved
    _DEFAULT_LEVEL = int(level)
    # Existing loggers were built before the config was read; move them too.
    for name in list(logging.root.manager.loggerDict):
        if name.startswith("cassa"):
            logging.getLogger(name).setLevel(_DEFAULT_LEVEL)
    return _DEFAULT_LEVEL


def get_logger(name="cassa_photometry", run_dir=None, level=None):
    """Return a configured logger.

    Parameters
    ----------
    name : str
        Logger name (shared across a phase so handlers are reused).
    run_dir : str, optional
        If given, a ``<name>.log`` file handler is attached inside this
        directory in addition to the console handler.
    level : int, optional
        Logging level. Defaults to whatever :func:`set_default_level` was last
        given, which ``load_config`` sets from ``PipelineConfig.log_level``.

    Notes
    -----
    The function is idempotent: calling it repeatedly with the same ``name``
    will not stack duplicate handlers.
    """
    logger = logging.getLogger(name)
    logger.setLevel(_DEFAULT_LEVEL if level is None else level)
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
