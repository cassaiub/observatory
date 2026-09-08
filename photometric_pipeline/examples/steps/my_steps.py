"""Example custom pipeline steps.

Point a config file at one of these and it joins the plan like any built-in::

    phase3:
      steps:
        custom:
          write_bright_list: examples/steps/my_steps.py:write_bright_list
        order: [aperture_correction, zero_point, flux_calibration,
                catalog, psf_photometry, classification, write_bright_list]

The point of this seam is that your code lives in a file upstream does not have,
so ``git pull`` keeps merging cleanly. Nothing here is imported by the pipeline
unless a config asks for it.

**Arguments.** A step is called with keyword arguments only, so declare the ones
you use and swallow the rest with ``**_``. That is what stops your step breaking
when the call site later grows another argument.

Phase 3 passes: ``engine``, ``path``, ``outdir``, ``base``, ``logger``.
Phase 1 passes: ``ccd``, ``meta``, ``state``, ``logger``.

**Ordering.** Set ``requires`` / ``provides`` on the function to join the
dependency graph. A step that declares neither is unconstrained and simply sits
where ``order`` puts it -- which is fine for something that only reads.
"""

import os


def write_bright_list(*, engine, outdir, base, logger, **_):
    """Write a short list of the brightest calibrated sources.

    Illustrates the common case: a step that runs *after* the catalog and reads
    what earlier steps produced.
    """
    import csv

    catalog_path = os.path.join(outdir, f"{base}_catalog.csv")
    if not os.path.exists(catalog_path):
        logger.warning("%s: no catalog to summarise; skipping bright list.", base)
        return

    from astropy.table import Table

    table = Table.read(catalog_path)

    # Pick the best column that actually has values. MAG_BEST is preferred, but
    # it is all-NaN when there is no zero point (an uncalibratable filter, or a
    # run with no network), and a step that silently produced an empty file
    # would be a poor example of anything. Falling back to the instrumental
    # magnitude keeps the step useful on an uncalibrated catalog.
    column = None
    for candidate in ("MAG_BEST", "MAG_AUTO", "MAG_INST"):
        if candidate in table.colnames and any(
                row[candidate] == row[candidate] for row in table):
            column = candidate
            break
    if column is None:
        logger.warning("%s: no usable magnitude column; skipping bright list.", base)
        return

    rows = [r for r in table if r[column] == r[column]]  # drop NaN
    rows.sort(key=lambda r: r[column])

    out_path = os.path.join(outdir, f"{base}_bright.csv")
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "x", "y", column])
        for rank, row in enumerate(rows[:20], 1):
            writer.writerow([rank, row["X_IMAGE"], row["Y_IMAGE"], row[column]])

    logger.info("%s: wrote %s (%d sources).", base, os.path.basename(out_path),
                min(20, len(rows)))


# This step reads the catalog, so it must not be ordered before the step that
# produces one. Declaring it makes that a load-time error rather than an
# empty file discovered later.
write_bright_list.requires = ("catalog",)


def clip_negatives(*, ccd, logger, **_):
    """A phase 1 example: clamp negative pixels to zero, in place.

    Deliberately something you should think twice about -- negative pixels after
    bias subtraction are real noise, and clipping them biases the background
    high. It is here because a *plausible but questionable* step is the honest
    illustration of what this seam is for: the pipeline lets you do it, records
    that you did, and does not pretend it was the default.
    """
    import numpy as np

    clipped = int(np.count_nonzero(ccd.data < 0))
    if clipped:
        ccd.data = np.clip(ccd.data, 0, None)
        logger.info("clip_negatives: clamped %d pixel(s).", clipped)


clip_negatives.requires = ("bias_subtracted",)
