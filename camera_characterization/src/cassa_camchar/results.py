"""Characterization results container with CSV / JSON persistence.

Persisting the measured numbers (not just a dashboard image) was the main
missing piece of the original script.
"""

import csv
import json
import math
from dataclasses import dataclass, field, asdict


@dataclass
class CharacterizationResults:
    """All measured sensor parameters."""

    driver_gains: list = field(default_factory=list)
    system_gain: list = field(default_factory=list)     # e-/ADU per driver gain
    read_noise_e: list = field(default_factory=list)
    full_well_e: list = field(default_factory=list)
    dynamic_range: list = field(default_factory=list)

    temperatures_c: list = field(default_factory=list)
    dark_current: list = field(default_factory=list)    # e-/pixel/s per temperature
    linearity: dict = field(default_factory=dict)       # slope / nonlinearity_pct / ...

    qe_peak_pct: float = float("nan")
    qe_peak_nm: float = float("nan")
    qe_synthetic: bool = False
    filter_metrics: dict = field(default_factory=dict)

    def write_json(self, path):
        with open(path, "w") as fh:
            json.dump(_json_safe(asdict(self)), fh, indent=2)

    def write_csv(self, path):
        """Write the per-gain headline table (full detail is in the JSON)."""
        cols = ["driver_gain", "system_gain_e_per_adu", "read_noise_e",
                "full_well_e", "dynamic_range_stops"]
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for i, g in enumerate(self.driver_gains):
                w.writerow([g, _fmt(self.system_gain[i]), _fmt(self.read_noise_e[i]),
                            _fmt(self.full_well_e[i]), _fmt(self.dynamic_range[i])])


def _fmt(v):
    return "" if v is None or (isinstance(v, float) and math.isnan(v)) else v


def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj
