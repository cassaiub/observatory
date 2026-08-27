"""The 9-panel sensor-characterization dashboard (headless)."""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

_FILTER_COLORS = {
    "Luminance": "gray", "Red": "#e74c3c", "Green": "#2ecc71", "Blue": "#3498db",
    "H-Alpha": "darkred", "OIII": "teal", "SII": "purple",
}


def generate_dashboard(master_ptc, results, spectra, linearity, out_path, show=False):
    """Render the 3x3 characterization dashboard to ``out_path`` (PNG)."""
    style = "seaborn-v0_8-whitegrid"
    plt.style.use(style if style in plt.style.available else "default")
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    dg = np.asarray(results.driver_gains)

    # 1: PTC
    if master_ptc is not None and len(master_ptc["means"]):
        axes[0, 0].loglog(master_ptc["means"], master_ptc["variances"], "bo-", markersize=4)
    axes[0, 0].set(title="Photon Transfer Curve", xlabel="Mean Signal (ADU)", ylabel="Signal Variance (ADU$^2$)")
    axes[0, 0].set_title("Photon Transfer Curve", fontweight="bold")

    # 2-5: gain-sweep metrics
    _line(axes[0, 1], dg, results.system_gain, "System Gain", "e-/ADU", "go-")
    _line(axes[0, 2], dg, results.read_noise_e, "Readout Noise", "e-", "ro-")
    _line(axes[1, 0], dg, results.full_well_e, "Full Well Capacity", "e-", "mo-")
    _line(axes[1, 1], dg, results.dynamic_range, "Dynamic Range", "Stops", "co-")
    for ax in (axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]):
        ax.set_xlabel("Camera Driver Gain Setting")

    # 6: dark current vs temperature
    t, dc = np.asarray(results.temperatures_c), np.asarray(results.dark_current)
    if t.size and np.any(dc > 0):
        axes[1, 2].semilogy(t, np.clip(dc, 1e-6, None), "ko-", linewidth=2)
    axes[1, 2].set_title("Dark Current", fontweight="bold")
    axes[1, 2].set(xlabel="Sensor Temperature (°C)", ylabel="Dark Current (e-/pixel/s)")

    # 7: QE
    qe = spectra["qe"]
    axes[2, 0].plot(qe["wavelengths"], qe["qe_percentage"], "b-", linewidth=2.5)
    axes[2, 0].fill_between(qe["wavelengths"], qe["qe_percentage"], alpha=0.2, color="blue")
    axes[2, 0].set_title("Absolute Quantum Efficiency" + (" (synthetic)" if spectra["qe_synthetic"] else ""),
                        fontweight="bold")
    axes[2, 0].set(xlabel="Wavelength (nm)", ylabel="Quantum Efficiency (%)", ylim=(0, 100), xlim=(350, 1000))

    # 8: filter transmission
    ax = axes[2, 1]
    for name, f in spectra["filters"].items():
        c = _FILTER_COLORS.get(name, "black")
        ax.plot(f["wavelengths"], f["transmission"], color=c, linewidth=2, label=name)
        ax.fill_between(f["wavelengths"], f["transmission"], alpha=0.15, color=c)
    if spectra["filters"]:
        ax.legend(loc="lower center", fontsize=8, ncol=4, bbox_to_anchor=(0.5, -0.35))
    ax.set_title("Filter Transmission Curves", fontweight="bold")
    ax.set(xlabel="Wavelength (nm)", ylabel="Transmission (%)", ylim=(0, 100), xlim=(350, 850))

    # 9: dark linearity
    ax = axes[2, 2]
    if linearity is not None:
        e, y = linearity["exposures"], linearity["thermal_electrons"]
        ax.plot(e, y, "ko-", linewidth=2, label="Measured")
        ax.plot(e, linearity["slope"] * e + linearity["intercept"], "r--", alpha=0.8,
                label=f"Linear fit ({linearity['nonlinearity_pct']:.1f}% dev)")
        ax.legend(loc="upper left", fontsize=9)
    else:
        ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                transform=ax.transAxes, color="gray", style="italic")
    ax.set_title("Dark Linearity (Amp Glow Check)", fontweight="bold")
    ax.set(xlabel="Exposure Time (s)", ylabel="Total Thermal Signal (e-)")

    plt.tight_layout(pad=3.0)
    fig.savefig(out_path, dpi=200)
    if show:
        plt.show()
    plt.close(fig)


def _line(ax, x, y, title, ylabel, style):
    ax.plot(x, y, style, linewidth=2)
    ax.set_title(title, fontweight="bold")
    ax.set_ylabel(ylabel)
