"""cassa_dimm: a professional-grade Differential Image Motion Monitor (DIMM).

Continuously watches a folder of incoming frames, detects and pairs every star's
two prism spots, measures the differential image motion, and reports airmass-
corrected atmospheric seeing as a live time series.
"""

__version__ = "0.1.0"

from cassa_dimm.config import DimmConfig, load_config

__all__ = ["DimmConfig", "load_config", "__version__"]
