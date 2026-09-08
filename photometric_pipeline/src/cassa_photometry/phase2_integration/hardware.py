"""CPU and memory budgeting for the integration phase.

Phase 2 is the memory-hungry phase -- it holds a stack of frames plus their
variance planes -- so it is worth telling the user what it is about to need. The
prompt is a convenience, never a requirement: a run with no terminal attached
(a batch job, a SLURM script, ``cassa-run``, a notebook, a test) takes the
default and says so, rather than blocking on ``input()`` or dying with
``EOFError``.
"""

import math
import os
import sys

import psutil

from cassa_photometry.logging_utils import get_logger

#: Fraction of the machine's cores per menu choice.
CPU_LEVELS = {"1": 0.25, "2": 0.50, "3": 0.75, "4": 1.0}
DEFAULT_LEVEL = "2"

#: Thread-count variables the allocation is applied to.
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


class HardwareManager:
    """Chooses a core count and reports the expected resource envelope."""

    def __init__(self, logger=None):
        self.total_cores = os.cpu_count() or 1
        self.sys_ram_gb = psutil.virtual_memory().total / (1024**3)
        self.allocated_cores = 1
        self.logger = logger or get_logger("cassa_integrate")

    def _interactive(self):
        """True only when there is a real terminal to prompt on."""
        try:
            return sys.stdin is not None and sys.stdin.isatty()
        except (AttributeError, ValueError):  # closed or replaced stream
            return False

    def allocate_resources(self, assume_yes=False, cpu_level=None):
        """Pick a core count, prompting only when that is possible and wanted.

        Parameters
        ----------
        assume_yes : bool
            Skip the prompt and take the default. Set by ``-y/--yes``.
        cpu_level : str, optional
            An explicit level ("1".."4"), bypassing the prompt entirely.
        """
        self.logger.info("Hardware: %d cores, %.1f GB RAM", self.total_cores, self.sys_ram_gb)

        choice = cpu_level
        if choice is None:
            if assume_yes or not self._interactive():
                choice = DEFAULT_LEVEL
                if not assume_yes:
                    self.logger.info(
                        "Not attached to a terminal; using the default CPU level "
                        "(%d%%). Pass cpu_level to choose another.",
                        int(CPU_LEVELS[DEFAULT_LEVEL] * 100),
                    )
            else:
                try:
                    choice = input(
                        "\nSelect CPU Resource Level "
                        "(1: 25%, 2: 50%, 3: 75%, 4: 100%) [Default: 2]: "
                    ).strip()
                except (EOFError, KeyboardInterrupt):
                    choice = DEFAULT_LEVEL

        multiplier = CPU_LEVELS.get(str(choice), CPU_LEVELS[DEFAULT_LEVEL])
        self.allocated_cores = max(1, int(self.total_cores * multiplier))

        # Scoped to this process and whatever it spawns. Recorded, because a
        # reduction that quietly reconfigures BLAS threading is worth a log line.
        for variable in THREAD_VARIABLES:
            os.environ[variable] = str(self.allocated_cores)
        self.logger.info(
            "Using %d of %d cores (%d%%).",
            self.allocated_cores, self.total_cores, int(multiplier * 100),
        )
        return self.allocated_cores

    def print_estimations(self, num_files, max_stack, total_size_mb):
        """Log the expected peak memory and runtime.

        A rough envelope, not a measurement: peak memory is driven by holding the
        largest stack's data and variance planes at once, and the time estimate
        assumes sub-linear scaling with cores.
        """
        avg_file_mb = total_size_mb / num_files if num_files > 0 else 0
        peak_ram_gb = ((avg_file_mb * 2) * max_stack * 5) / 1024 + 1.5

        base_time = 3.5
        thread_eff = (
            math.sqrt(self.allocated_cores / 4) if self.allocated_cores > 4
            else (self.allocated_cores / 4)
        )
        est_time_s = (num_files * base_time) / max(thread_eff, 1e-6)

        self.logger.info(
            "Resource estimate: %d frames, largest stack %d, peak RAM ~%.1f GB, "
            "runtime ~%dm %02ds.",
            num_files, max_stack, peak_ram_gb,
            int(est_time_s // 60), int(est_time_s % 60),
        )
        if peak_ram_gb > self.sys_ram_gb:
            self.logger.warning(
                "Estimated peak RAM (%.1f GB) exceeds system memory (%.1f GB); "
                "the run may swap heavily or be killed.",
                peak_ram_gb, self.sys_ram_gb,
            )
