"""Nominal PDA-23 honeycomb layout (rows 5-4-5-4-5, channel order = routing 9..31 row-major,
as in SPADarray_AberrationPSF/detectorLayout.m). Units: detector pitch."""
import numpy as np


def honeycomb23():
    rows = [5, 4, 5, 4, 5]
    yp = np.sqrt(3) / 2
    xy = []
    for r, n in enumerate(rows):
        y = (2 - r) * yp
        xs = (np.arange(n) - (n - 1) / 2)
        xy += [(x, y) for x in xs]
    return np.array(xy)
