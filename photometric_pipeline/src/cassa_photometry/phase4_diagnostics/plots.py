"""Matplotlib helpers for the diagnostics report.

Follows the visual style of ``phase2_integration/visuals.py`` (ZScale display,
PdfPages multi-page output) and adds reusable panels for histograms, radial
profiles, scatter plots and text summaries.
"""

import matplotlib
import numpy as np

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from matplotlib.backends.backend_pdf import PdfPages

_ZSCALE = ZScaleInterval()


class Report:
    """Collects figures into a multi-page PDF (and optional per-page PNGs)."""

    def __init__(self, pdf_path):
        self.pdf_path = pdf_path
        self._pdf = PdfPages(pdf_path)

    def add(self, fig, png_path=None):
        self._pdf.savefig(fig, bbox_inches="tight")
        if png_path:
            fig.savefig(png_path, dpi=120, bbox_inches="tight")
        plt.close(fig)

    def close(self):
        self._pdf.close()


def zscale_panel(ax, data, title, cmap="gray"):
    """Display an image with ZScale limits (FITS orientation)."""
    finite = data[np.isfinite(data)]
    if finite.size:
        vmin, vmax = _ZSCALE.get_limits(finite)
    else:
        vmin, vmax = 0, 1
    im = ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    return im


def image_panel(ax, data, title, cmap="viridis", vmin=None, vmax=None):
    im = ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    return im


def overlay_positions(ax, positions, color="cyan", s=60):
    """Circle detected/used star positions on an image panel."""
    if not len(positions):
        return
    xy = np.asarray(positions)
    ax.scatter(xy[:, 0], xy[:, 1], s=s, facecolors="none", edgecolors=color, linewidths=0.8)


def histogram_panel(ax, data, title, bins=200, log=True):
    finite = data[np.isfinite(data)]
    if finite.size:
        lo, hi = np.percentile(finite, [0.5, 99.5])
        ax.hist(finite, bins=bins, range=(lo, hi), color="steelblue", log=log)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Pixel value")
    ax.set_ylabel("Count")


def radial_profile_panel(ax, radii, profile, fwhm_px=None, title="Radial profile"):
    if radii is None or profile is None:
        _text_only(ax, "No radial profile\n(insufficient stars)", title)
        return
    ax.plot(radii, profile, "o-", color="darkorange", ms=3)
    ax.axhline(0.5, color="gray", ls=":", lw=0.8)
    if fwhm_px and np.isfinite(fwhm_px):
        ax.axvline(fwhm_px / 2.0, color="red", ls="--", lw=0.8, label=f"HWHM={fwhm_px/2:.2f}px")
        ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Radius (px)")
    ax.set_ylabel("Normalised flux")
    ax.set_ylim(-0.05, 1.05)


def scatter_panel(ax, x, y, xlabel, ylabel, title, logy=False, c=None, cmap="viridis", s=8):
    sc = ax.scatter(x, y, c=c, cmap=cmap, s=s, alpha=0.7)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    return sc


def bar_panel(ax, labels, values, title, color="indianred"):
    ax.bar(range(len(labels)), values, color=color)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("Pixels")


def text_panel(ax, lines, title="Summary"):
    """Render a metrics summary as monospaced text."""
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left",
            family="monospace", fontsize=9, transform=ax.transAxes)


def _text_only(ax, message, title):
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10,
            transform=ax.transAxes, color="gray")


def new_page(nrows, ncols, suptitle, figsize=None):
    """Create a figure/axes grid with a bold page title."""
    figsize = figsize or (5.2 * ncols, 4.6 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    fig.suptitle(suptitle, fontsize=14, fontweight="bold", y=0.99)
    return fig, np.atleast_1d(axes).ravel()
