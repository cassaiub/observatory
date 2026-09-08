"""The detector half of the forward model: photons in, ADU out.

Deliberately the exact inverse of what phase 1 undoes, in the same order, so a
correct reduction recovers the input and an incorrect one shows up as a
measurable bias rather than a plausible-looking image:

    (sky + sources) x flat x t   ->  electrons
      + dark x t + amp glow      ->  electrons
      Poisson                    ->  electrons
      / gain + bias + read noise ->  ADU
      non-linearity              ->  ADU
      clip at full well          ->  ADU
      quantise                   ->  uint16

Two details are what make simulated data worth testing against:

* **Dark current is added after the flat**, because it is generated in the
  silicon and not attenuated by the optics. Flat-fielding it as if it were sky
  is a real and common reduction error, and a simulator that also gets it wrong
  cannot detect it.
* **Non-linearity is applied last, in ADU**, growing with signal. That is what
  makes it tilt a magnitude scale rather than offset it -- the defect a
  saturation flag does not catch.
"""

import numpy as np


class DetectorModel:
    """A CMOS detector: gain, noise, dark current, flat field, non-linearity.

    Defaults describe the CASSA 8-inch's QHY miniCAM8M (IMX585) closely enough
    to be representative; they are not measurements of that camera, and the
    simulator says so in the headers it writes.
    """

    def __init__(
        self,
        shape,
        gain=0.15,
        read_noise_e=1.3,
        bias_level=500.0,
        dark_rate_e_per_s=0.005,
        full_well_adu=58000.0,
        nonlinearity=0.02,
        hot_pixel_fraction=2e-5,
        dead_pixel_fraction=5e-6,
        vignetting=0.18,
        rng=None,
    ):
        self.shape = tuple(shape)
        self.gain = float(gain)
        self.read_noise_e = float(read_noise_e)
        self.bias_level = float(bias_level)
        self.dark_rate = float(dark_rate_e_per_s)
        self.full_well_adu = float(full_well_adu)
        self.nonlinearity = float(nonlinearity)
        self.rng = rng if rng is not None else np.random.default_rng(0)

        self.flat = self._build_flat(vignetting)
        self.dark_current = self._build_dark(hot_pixel_fraction)
        self.bias_structure = self._build_bias()
        self.dead_mask = self.rng.random(self.shape) < dead_pixel_fraction
        self.flat[self.dead_mask] = 0.02

    # --- Fixed-pattern structure ---------------------------------------------
    def _build_flat(self, vignetting):
        """Relative response: vignetting x dust donuts x pixel-to-pixel scatter.

        The large-scale term is what a *sky* flat corrects and a dome flat only
        partly does, so it is what makes the flat-type question testable.
        """
        ny, nx = self.shape
        y, x = np.mgrid[0:ny, 0:nx]
        cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
        radius = np.hypot(y - cy, x - cx) / np.hypot(cy, cx)

        flat = 1.0 - vignetting * radius**2
        for _ in range(6):  # dust donuts on the window
            dy = self.rng.uniform(0, ny)
            dx = self.rng.uniform(0, nx)
            r0 = self.rng.uniform(0.02, 0.05) * max(ny, nx)
            depth = self.rng.uniform(0.05, 0.15)
            d = np.hypot(y - dy, x - dx)
            flat *= 1.0 - depth * np.exp(-0.5 * ((d - r0) / (0.25 * r0)) ** 2)
        flat *= self.rng.normal(1.0, 0.01, self.shape)   # pixel response
        return np.clip(flat, 0.02, None)

    def _build_dark(self, hot_pixel_fraction):
        """Dark current per pixel in e-/s, with a hot-pixel tail."""
        dark = np.full(self.shape, self.dark_rate, dtype=float)
        n_hot = int(hot_pixel_fraction * dark.size)
        if n_hot:
            idx = self.rng.integers(0, dark.size, n_hot)
            dark.flat[idx] = self.rng.uniform(50, 500, n_hot) * self.dark_rate
        return dark

    def _build_bias(self):
        """Bias structure in ADU: a column pattern plus a gentle gradient."""
        ny, nx = self.shape
        columns = self.rng.normal(0.0, 1.5, nx)
        gradient = np.linspace(0.0, 2.0, ny)[:, None]
        return columns[None, :] + gradient

    # --- The forward model ----------------------------------------------------
    def _apply_nonlinearity(self, adu):
        """Compress the response toward full well.

        A simple quadratic: response falls short by ``nonlinearity`` (a few
        percent) at full well and by nothing at zero. Brightness-dependent by
        construction, which is the point.
        """
        if self.nonlinearity <= 0:
            return adu
        fraction = np.clip(adu / self.full_well_adu, 0.0, None)
        return adu * (1.0 - self.nonlinearity * fraction)

    def _digitise(self, electrons, exposure, add_dark=True):
        """Shared tail of every frame type: noise, gain, bias, clip, quantise."""
        if add_dark:
            electrons = electrons + self.dark_current * exposure

        noisy = self.rng.poisson(np.clip(electrons, 0, None)).astype(float)
        adu = noisy / self.gain
        adu += self.bias_level + self.bias_structure
        adu += self.rng.normal(0.0, self.read_noise_e / self.gain, self.shape)
        adu = self._apply_nonlinearity(adu)
        adu = np.clip(adu, 0.0, self.full_well_adu)
        return np.rint(adu).astype(np.uint16)

    def expose_science(self, photon_rate, exposure):
        """A science frame from a scene in e-/s (before the flat)."""
        return self._digitise(photon_rate * self.flat * exposure, exposure)

    def expose_bias(self):
        return self._digitise(np.zeros(self.shape), 0.0, add_dark=False)

    def expose_dark(self, exposure):
        return self._digitise(np.zeros(self.shape), exposure)

    def expose_flat(self, level_adu, exposure):
        """A flat frame illuminated to roughly ``level_adu`` at the field centre.

        Taking the level in ADU rather than electrons matches how an observer
        actually judges a flat -- "expose to about half full well" -- and keeps
        the caller from having to know the gain.
        """
        level_e = float(level_adu) * self.gain
        return self._digitise(np.full(self.shape, level_e) * self.flat, exposure)
