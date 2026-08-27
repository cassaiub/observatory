"""Console entry points.

* ``cassa-camchar-analyze`` -- characterize a data directory (dashboard + CSV/JSON).
* ``cassa-camchar-sim``     -- generate a synthetic calibration campaign.
"""

import argparse

from cassa_camchar import __version__
from cassa_camchar.config import load_config
from cassa_camchar.logging_utils import get_logger


def analyze():
    p = argparse.ArgumentParser(prog="cassa-camchar-analyze",
                                description="Characterize a CMOS sensor from a calibration data directory.")
    p.add_argument("data_dir", nargs="?", default=None, help="Directory with bias/flat/dark/spectra (default: config).")
    p.add_argument("-o", "--outdir", default=None, help="Output directory for results + dashboard.")
    p.add_argument("--show", action="store_true", help="Also display the dashboard window.")
    p.add_argument("-c", "--config", default=None, help="YAML config file.")
    p.add_argument("--version", action="version", version=f"cassa-camchar {__version__}")
    args = p.parse_args()

    from cassa_camchar.pipeline import analyze as run
    run(data_dir=args.data_dir, config=load_config(args.config), out_dir=args.outdir, show=args.show)


def simulate():
    p = argparse.ArgumentParser(prog="cassa-camchar-sim",
                                description="Generate a synthetic sensor-characterization campaign.")
    p.add_argument("--out", required=True, help="Output data directory (bias/flat/dark/spectra created inside).")
    p.add_argument("--seed", type=int, default=None, help="RNG seed.")
    p.add_argument("-c", "--config", default=None, help="YAML config file.")
    args = p.parse_args()

    from cassa_camchar.simulate import simulate_campaign
    simulate_campaign(args.out, config=load_config(args.config), seed=args.seed,
                      logger=get_logger("cassa_camchar_sim"))
