"""cassa_camchar: CMOS sensor characterization via the Photon Transfer Curve.

Measures system gain, read noise, full-well capacity, dynamic range, dark
current, dark linearity, quantum efficiency, and filter transmission, and
writes a results table (CSV/JSON) plus a 9-panel dashboard.
"""

__version__ = "0.1.0"

from cassa_camchar.config import CamCharConfig, load_config

__all__ = ["CamCharConfig", "load_config", "__version__"]
