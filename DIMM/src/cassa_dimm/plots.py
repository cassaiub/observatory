"""Live seeing time-series plot for the monitor."""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def plot_timeseries(times, seeing, seeing_err, out_path, title="DIMM Seeing Monitor",
                    max_points=720):
    """Render a seeing-vs-time plot (with error band and running median)."""
    if not times:
        return
    times = list(times)[-max_points:]
    seeing = np.asarray(seeing[-max_points:], dtype=float)
    err = np.asarray(seeing_err[-max_points:], dtype=float)
    err = np.where(np.isfinite(err), err, 0.0)

    fig, ax = plt.subplots(figsize=(11, 4.5), dpi=110)
    ax.plot(times, seeing, "o-", ms=3, lw=1.0, color="steelblue", label="seeing")
    ax.fill_between(times, seeing - err, seeing + err, color="steelblue", alpha=0.2)

    finite = seeing[np.isfinite(seeing)]
    if finite.size:
        med = float(np.median(finite))
        ax.axhline(med, color="firebrick", ls="--", lw=0.9, label=f"median {med:.2f}\"")

    ax.set_ylabel("Seeing (arcsec, zenith)")
    ax.set_xlabel("Time (UTC)")
    ax.set_title(title, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper right")
    try:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        fig.autofmt_xdate()
    except Exception:
        pass
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
