"""Star/galaxy separation, calibrated from each frame's own stellar locus.

The pipeline used to classify on ``ELLIPTICITY < 0.15``. Ellipticity is *shape*,
not *concentration*, so a face-on elliptical galaxy is round and passes as a star
while a slightly trailed star is elongated and fails. That mattered beyond the
catalog column: a galaxy leaking into the zero point biases every magnitude in
the frame.

The discriminant here is **``MAG_PSF - MAG_AUTO``**, which is what Pan-STARRS
uses (``iPSFMag - iKronMag``). A point source is fitted well by the PSF and its
Kron aperture adds nothing, so the difference is near zero; an extended source
has more light outside the PSF core, so the difference is negative. The
threshold is **calibrated per frame from the stellar locus itself** -- bright,
high-SNR sources define where point sources sit tonight, at tonight's seeing --
so there is no training data, no runtime catalog dependency, and no fixed cut
that silently goes wrong when the seeing changes.

Two honesty requirements shape the rest:

* Below the SNR where the locus separation falls under its own scatter, the
  answer is ``AMBIGUOUS``. The magnitude at which that happens is reported as
  ``CLASSLIM``, rather than pretending to classify to the detection limit.
* A field too sparse to define a locus returns ``AMBIGUOUS`` for everything with
  a loud warning -- never a silent fixed cut.
"""

import numpy as np
from astropy.stats import sigma_clipped_stats

STAR = "STAR"
EXTENDED = "EXTENDED"
AMBIGUOUS = "AMBIGUOUS"
BLENDED = "BLENDED"
SATURATED = "SATURATED"
EDGE = "EDGE"

#: Fewest locus stars before a per-frame calibration is trusted.
MIN_LOCUS_STARS = 8

#: SNR above which a source is bright enough to help define the locus.
LOCUS_MIN_SNR = 20.0

#: SNR below which no class is asserted at all: the PSF-minus-Kron difference is
#: then dominated by measurement noise rather than by source structure.
MIN_CLASSIFY_SNR = 10.0


class Classification:
    """The outcome of classifying one frame's sources."""

    def __init__(self, classes, stellarity, limit_mag=None, locus_centre=None,
                 locus_scatter=None, method="locus", n_locus=0):
        self.classes = np.asarray(classes)
        #: Continuous stellarity, 1 = point-like, 0 = clearly extended.
        self.stellarity = np.asarray(stellarity, dtype=float)
        #: Magnitude below which the classes stop being meaningful.
        self.limit_mag = limit_mag
        self.locus_centre = locus_centre
        self.locus_scatter = locus_scatter
        self.method = method
        self.n_locus = n_locus


def classify(psf_minus_auto, magnitude, snr, flags=None, saturated=None,
             edge=None, blended=None, logger=None):
    """Classify sources from their PSF-minus-Kron difference.

    Parameters
    ----------
    psf_minus_auto : ndarray
        ``MAG_PSF - MAG_AUTO`` per source.
    magnitude : ndarray
        A calibrated or instrumental magnitude, used only to report the limit.
    snr : ndarray
        Signal-to-noise, used to pick the locus sources and to decide where
        classification stops being meaningful.
    """
    difference = np.asarray(psf_minus_auto, dtype=float)
    magnitude = np.asarray(magnitude, dtype=float)
    snr = np.asarray(snr, dtype=float)
    n = difference.size

    classes = np.full(n, AMBIGUOUS, dtype=object)
    stellarity = np.full(n, np.nan)

    usable = np.isfinite(difference)
    bright = usable & np.isfinite(snr) & (snr >= LOCUS_MIN_SNR)

    if bright.sum() < MIN_LOCUS_STARS:
        if logger is not None:
            logger.warning(
                "Only %d source(s) are bright enough to define a stellar locus "
                "(need %d). Every source is reported AMBIGUOUS rather than "
                "classified against a fixed cut that may not suit this frame.",
                int(bright.sum()), MIN_LOCUS_STARS,
            )
        _apply_overrides(classes, saturated, edge, blended)
        return Classification(classes, stellarity, method="insufficient-locus",
                              n_locus=int(bright.sum()))

    centre, scatter = _stellar_locus(difference[bright])

    # Extended sources have MORE light outside the PSF core, so MAG_PSF is
    # fainter than MAG_AUTO and the difference runs positive.
    deviation = np.where(usable, (difference - centre) / scatter, np.nan)
    with np.errstate(invalid="ignore"):
        stellarity = np.clip(1.0 - np.abs(deviation) / 6.0, 0.0, 1.0)

    limit_mag = _classification_limit(magnitude, snr)
    too_faint = np.isfinite(magnitude) & np.isfinite(limit_mag) & (magnitude > limit_mag)
    # An SNR floor as well as a magnitude limit. Below it the PSF-minus-Kron
    # difference is dominated by measurement noise, so a confident class would
    # be read off nothing -- and it is the *faint* end where a spurious
    # "EXTENDED" would be most misleading.
    too_noisy = np.isfinite(snr) & (snr < MIN_CLASSIFY_SNR)

    classes[usable & (np.abs(deviation) <= 3.0)] = STAR
    classes[usable & (deviation > 3.0)] = EXTENDED
    # A source far on the *compact* side is not a galaxy; it is more likely a
    # cosmic ray or an artefact, so it is not called a star either.
    classes[usable & (deviation < -5.0)] = AMBIGUOUS
    classes[too_faint | too_noisy] = AMBIGUOUS

    _apply_overrides(classes, saturated, edge, blended)

    if logger is not None:
        counts = {c: int(np.count_nonzero(classes == c)) for c in set(classes.tolist())}
        logger.info(
            "    Classification: locus at %+.3f +/- %.3f mag from %d bright "
            "sources; %s. Meaningful to mag %.2f.",
            centre, scatter, int(bright.sum()),
            ", ".join(f"{k} {v}" for k, v in sorted(counts.items())),
            limit_mag if np.isfinite(limit_mag) else float("nan"),
        )
    return Classification(classes, stellarity, limit_mag=limit_mag,
                          locus_centre=centre, locus_scatter=scatter,
                          n_locus=int(bright.sum()))


def _stellar_locus(values):
    """Where the point sources sit, and how tightly.

    Found as the **mode**, not the sigma-clipped mean. Sigma-clipping locates
    the bulk of a distribution, so in a field where galaxies outnumber stars it
    settles on the galaxies and the classification inverts. Point sources form a
    narrow peak at the compact end, and the mode finds that peak whatever the
    mixture -- verified against a population of 20 stars among 80 galaxies,
    where sigma-clipping put the locus at +0.54 instead of 0.0.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < MIN_LOCUS_STARS:
        return 0.0, 0.02

    # Histogram peak, on a range robust to a long extended-source tail.
    low, high = np.percentile(values, [1, 99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return float(np.median(values)), 0.02
    counts, edges = np.histogram(values, bins=max(int(np.sqrt(values.size) * 2), 10),
                                 range=(low, high))
    peak = 0.5 * (edges[np.argmax(counts)] + edges[np.argmax(counts) + 1])

    # Refine on the sources near that peak, so the width is the locus's own and
    # not the whole population's. Iterated, because a single fixed window
    # truncates a broad locus and reports it as narrower than it is -- which
    # then makes every genuine star look like an outlier.
    centre = peak
    width = max((high - low) / 10.0, 1e-3)
    scatter = width / 2.0
    for _ in range(4):
        near = np.abs(values - centre) < width
        if near.sum() < MIN_LOCUS_STARS:
            break
        _, centre, measured = sigma_clipped_stats(values[near], sigma=2.5, maxiters=5)
        if not np.isfinite(measured) or measured <= 1e-4:
            break
        scatter = float(measured)
        new_width = 3.0 * scatter
        if abs(new_width - width) < 0.05 * width:
            width = new_width
            break
        width = new_width

    scatter = float(scatter) if np.isfinite(scatter) and scatter > 1e-4 else 0.02
    return float(centre), scatter


def _apply_overrides(classes, saturated, edge, blended):
    """States that override any measurement, because the measurement is void."""
    for mask, label in ((blended, BLENDED), (edge, EDGE), (saturated, SATURATED)):
        if mask is not None:
            classes[np.asarray(mask, dtype=bool)] = label


def _classification_limit(magnitude, snr, target_snr=None):
    """The magnitude below which classification stops being meaningful.

    Taken as the magnitude at which SNR falls to the level where the locus
    separation is comparable with its own scatter. Reported rather than applied
    silently, so a catalog says how deep its classes can be trusted.
    """
    target_snr = target_snr or LOCUS_MIN_SNR / 2.0
    good = np.isfinite(magnitude) & np.isfinite(snr) & (snr > 0)
    if good.sum() < MIN_LOCUS_STARS:
        return np.nan
    near = good & (snr > target_snr * 0.7) & (snr < target_snr * 1.4)
    if near.sum() >= 3:
        return float(np.median(magnitude[near]))
    # Fall back to an interpolation through the sorted SNR/magnitude relation.
    order = np.argsort(-snr[good])
    sorted_snr = snr[good][order]
    sorted_mag = magnitude[good][order]
    below = np.where(sorted_snr <= target_snr)[0]
    return float(sorted_mag[below[0]]) if below.size else float(np.max(sorted_mag))
