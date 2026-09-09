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

    def __init__(self, logger=None, config=None):
        self.total_cores = os.cpu_count() or 1
        self.sys_ram_gb = psutil.virtual_memory().total / (1024**3)
        self.available_ram_gb = psutil.virtual_memory().available / (1024**3)
        self.allocated_cores = 1
        self.config = config
        self.logger = logger or get_logger("cassa_integrate")
        #: Filled in by print_estimations, so a caller (a notebook, a report)
        #: can show the same numbers the log line carries.
        self.estimate = {}

    # -- the configured budget -------------------------------------------------
    def _budget(self):
        """``(cpu_fraction, max_cores, max_memory_gb)`` from the config."""
        phase2 = getattr(self.config, "phase2", None)
        return (getattr(phase2, "cpu_fraction", None),
                getattr(phase2, "max_cores", None),
                getattr(phase2, "max_memory_gb", None))

    def _cap(self, cores):
        """Apply ``max_cores``, and say so when it bites."""
        _, max_cores, _ = self._budget()
        if max_cores and cores > int(max_cores):
            self.logger.info(
                "Capping %d cores to %d (phase2.max_cores).", cores, int(max_cores)
            )
            return max(1, int(max_cores))
        return cores

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

        # An explicitly configured budget is an answer: honour it and do not
        # prompt. Being asked again after setting it in a config file would be
        # the surprising behaviour, and it makes the setting useless in a
        # notebook or a batch job, which is where it is most wanted.
        fraction, _, _ = self._budget()
        if cpu_level is None and fraction is not None:
            fraction = min(max(float(fraction), 0.0), 1.0)
            self.allocated_cores = self._cap(max(1, round(self.total_cores * fraction)))
            self._apply_thread_limits()
            self.logger.info(
                "Using %d of %d cores (%d%%, from phase2.cpu_fraction).",
                self.allocated_cores, self.total_cores,
                int(round(100 * self.allocated_cores / self.total_cores)),
            )
            return self.allocated_cores

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
        self.allocated_cores = self._cap(max(1, int(self.total_cores * multiplier)))
        self._apply_thread_limits()
        self.logger.info(
            "Using %d of %d cores (%d%%).",
            self.allocated_cores, self.total_cores,
            int(round(100 * self.allocated_cores / self.total_cores)),
        )
        return self.allocated_cores

    def _apply_thread_limits(self):
        """Scoped to this process and whatever it spawns.

        Recorded in the log, because a reduction that quietly reconfigures BLAS
        threading is worth a line.
        """
        for variable in THREAD_VARIABLES:
            os.environ[variable] = str(self.allocated_cores)

    def print_estimations(self, num_files, max_stack, total_size_mb):
        """Log the expected resource envelope, and what of it will be used.

        A rough envelope, not a measurement: peak memory is driven by holding
        the largest stack's data and variance planes at once, and the time
        estimate assumes sub-linear scaling with cores.

        The numbers are also kept on ``self.estimate`` so a notebook or a report
        can show the same figures the log line carries, rather than recomputing
        them and drifting.
        """
        estimate = self.estimate_resources(num_files, max_stack, total_size_mb)

        self.logger.info(
            "Resource estimate: %d frames, largest stack %d, peak RAM ~%.1f GB, "
            "runtime ~%dm %02ds.",
            num_files, max_stack, estimate["peak_ram_gb"],
            int(estimate["runtime_s"] // 60), int(estimate["runtime_s"] % 60),
        )
        self.logger.info(
            "Using %d of %d cores (%d%%) and ~%.1f of %.1f GB RAM (%d%%).",
            estimate["cores_used"], estimate["cores_total"],
            round(estimate["cpu_percent"]),
            estimate["peak_ram_gb"], estimate["ram_total_gb"],
            round(estimate["ram_percent"]),
        )

        limit = estimate["ram_limit_gb"]
        if estimate["peak_ram_gb"] > limit:
            self.logger.warning(
                "Estimated peak RAM (%.1f GB) exceeds %s (%.1f GB); the run may "
                "swap heavily or be killed. Reduce phase2.max_cores, or stack "
                "fewer frames per group.",
                estimate["peak_ram_gb"], estimate["ram_limit_label"], limit,
            )
        return estimate

    def estimate_resources(self, num_files, max_stack, total_size_mb):
        """The resource envelope as a dict, without logging anything.

        Separated from :meth:`print_estimations` so the workshop notebooks can
        show a participant what the run will need *before* starting it, and how
        much of their own machine that is.
        """
        avg_file_mb = total_size_mb / num_files if num_files > 0 else 0
        peak_ram_gb = ((avg_file_mb * 2) * max_stack * 5) / 1024 + 1.5

        base_time = 3.5
        cores = max(self.allocated_cores, 1)
        thread_eff = math.sqrt(cores / 4) if cores > 4 else (cores / 4)
        runtime_s = (num_files * base_time) / max(thread_eff, 1e-6)

        _, _, max_memory_gb = self._budget()
        if max_memory_gb:
            ram_limit_gb = float(max_memory_gb)
            ram_limit_label = "the configured budget (phase2.max_memory_gb)"
        else:
            ram_limit_gb = self.sys_ram_gb
            ram_limit_label = "system memory"

        self.estimate = {
            "frames": num_files,
            "largest_stack": max_stack,
            "total_input_mb": total_size_mb,
            "cores_used": cores,
            "cores_total": self.total_cores,
            "cpu_percent": 100.0 * cores / max(self.total_cores, 1),
            "peak_ram_gb": peak_ram_gb,
            "ram_total_gb": self.sys_ram_gb,
            "ram_available_gb": self.available_ram_gb,
            "ram_percent": 100.0 * peak_ram_gb / max(self.sys_ram_gb, 1e-6),
            "ram_limit_gb": ram_limit_gb,
            "ram_limit_label": ram_limit_label,
            "fits_in_budget": peak_ram_gb <= ram_limit_gb,
            "runtime_s": runtime_s,
        }
        return self.estimate

    def describe(self):
        """A short human-readable summary of the budget, for notebooks."""
        e = self.estimate
        if not e:
            return (f"{self.total_cores} cores, {self.sys_ram_gb:.1f} GB RAM "
                    f"({self.available_ram_gb:.1f} GB free). No estimate yet.")
        verdict = "fits" if e["fits_in_budget"] else "DOES NOT FIT"
        return (
            f"This machine : {e['cores_total']} cores, {e['ram_total_gb']:.1f} GB RAM "
            f"({e['ram_available_gb']:.1f} GB free)\n"
            f"This run will : {e['cores_used']} cores ({e['cpu_percent']:.0f}%), "
            f"~{e['peak_ram_gb']:.1f} GB peak ({e['ram_percent']:.0f}%)\n"
            f"Workload      : {e['frames']} frames, largest stack {e['largest_stack']}, "
            f"~{e['runtime_s'] / 60:.1f} min\n"
            f"Against {e['ram_limit_label']} ({e['ram_limit_gb']:.1f} GB): {verdict}"
        )
