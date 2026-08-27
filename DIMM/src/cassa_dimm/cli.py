"""Console entry points for the DIMM pipeline.

* ``cassa-dimm-monitor`` -- continuous watch-folder seeing monitor.
* ``cassa-dimm-batch``   -- one-shot reduction of a cube/folder.
* ``cassa-dimm-sim``     -- generate synthetic frames with a known seeing.
"""

import argparse

from cassa_dimm import __version__
from cassa_dimm.config import load_config
from cassa_dimm.logging_utils import get_logger


def monitor():
    p = argparse.ArgumentParser(prog="cassa-dimm-monitor",
                                description="Continuous DIMM seeing monitor (watch folder).")
    p.add_argument("-c", "--config", default=None, help="YAML config file.")
    p.add_argument("--folder", default=None, help="Override watch folder.")
    p.add_argument("--outdir", default=None, help="Override output/log directory.")
    p.add_argument("--version", action="version", version=f"cassa-dimm {__version__}")
    args = p.parse_args()

    config = load_config(args.config)
    if args.folder:
        config.watch.folder = args.folder
    if args.outdir:
        config.output.log_dir = args.outdir

    from cassa_dimm.monitor import DimmMonitor
    DimmMonitor(config=config).run()


def batch():
    p = argparse.ArgumentParser(prog="cassa-dimm-batch",
                                description="One-shot DIMM reduction of a cube or frame directory.")
    p.add_argument("path", help="A FITS cube file, or a directory of frames.")
    p.add_argument("-c", "--config", default=None, help="YAML config file.")
    args = p.parse_args()

    from cassa_dimm.monitor import run_batch
    run_batch(args.path, config=load_config(args.config))


def simulate():
    p = argparse.ArgumentParser(prog="cassa-dimm-sim",
                                description="Generate synthetic DIMM data with a known seeing.")
    p.add_argument("--out", required=True, help="Output directory (or cube dir with --cube).")
    p.add_argument("--seeing", type=float, default=1.5, help="Injected seeing (arcsec). Default 1.5.")
    p.add_argument("--nstars", type=int, default=3, help="Number of doublets. Default 3.")
    p.add_argument("--frames", type=int, default=500, help="Number of frames. Default 500.")
    p.add_argument("--altitude", type=float, default=90.0, help="Target altitude (deg). Default 90.")
    p.add_argument("--cube", action="store_true", help="Write one cube instead of per-frame files.")
    p.add_argument("--seed", type=int, default=None, help="RNG seed.")
    p.add_argument("-c", "--config", default=None, help="YAML config file.")
    args = p.parse_args()

    from cassa_dimm.simulate import simulate_frames
    simulate_frames(args.out, seeing_arcsec=args.seeing, n_stars=args.nstars,
                    n_frames=args.frames, altitude_deg=args.altitude, as_cube=args.cube,
                    seed=args.seed, config=load_config(args.config),
                    logger=get_logger("cassa_dimm_sim"))
