"""Step 2b: nanosecond cross-correlation of every channel pair (crosstalk check).

With 40 MHz pulsed excitation the pair coincidence histogram is a comb of
peaks at multiples of 25 ns. For independent pixels the zero-delay peak is
smaller than the side peaks by the antibunching factor (1 - 1/N_ij) of the
pair. Optical/electrical crosstalk adds a narrow extra peak at zero delay.
Excess counts over the antibunching expectation, divided by the counts of
the triggering pixel, give the crosstalk probability.
"""
import sys, numpy as np
from pathlib import Path

CACHE = Path(sys.argv[1]); OUT = Path("out")
NCH = 23; R0 = 9
PER = 25.0          # ns
WIN = 4.0           # ns half-width of a comb peak window
ns = np.load(CACHE / "atto655_nsync.npy", mmap_mode="r")
ch = np.load(CACHE / "atto655_chan.npy", mmap_mode="r")
dt = np.load(CACHE / "atto655_dtime.npy", mmap_mode="r")
step1 = np.load(OUT / "step1_tcspc.npz")
offs = step1["peak"] - np.median(step1["peak"])   # per-channel timing offset (ns)

NPEAK = 3                                  # peaks -3..+3 x 25 ns
edges = np.arange(-(NPEAK + 0.5) * PER, (NPEAK + 0.5) * PER + 1e-9, 0.1)
hist_all = np.zeros(edges.size - 1)
peakcnt = np.zeros((2 * NPEAK + 1, NCH, NCH))
chunk = 10_000_000
K = 4                                      # neighbours ahead considered
for i0 in range(0, ns.size, chunk):
    i1 = min(i0 + chunk + K, ns.size)
    c = np.asarray(ch[i0:i1]).astype(np.int64) - R0
    t = np.asarray(ns[i0:i1]) * PER + np.asarray(dt[i0:i1]) * 0.02 - offs[np.clip(c, 0, NCH - 1)]
    for k in range(1, K + 1):
        d = t[k:] - t[:-k]
        a = c[:-k]; b = c[k:]
        m = d < (NPEAK + 0.5) * PER
        d, a, b = d[m], a[m], b[m]
        # symmetric: (a->b at +d) and (b->a at -d)
        hist_all += np.histogram(np.concatenate([d, -d]), edges)[0]
        for p in range(-NPEAK, NPEAK + 1):
            sel = np.abs(d - p * PER) < WIN
            np.add.at(peakcnt[NPEAK + p], (a[sel], b[sel]), 1)
            if p != 0:
                np.add.at(peakcnt[NPEAK - p], (b[sel], a[sel]), 1)
            else:
                np.add.at(peakcnt[NPEAK], (b[sel], a[sel]), 1)
np.savez(OUT / "step2b_crosstalk.npz", edges=edges, hist=hist_all, peakcnt=peakcnt)
side = 0.5 * (peakcnt[NPEAK - 1] + peakcnt[NPEAK + 1])
cen = peakcnt[NPEAK]
ratio = cen / np.maximum(side, 1)
np.fill_diagonal(ratio, np.nan)
print("zero-delay / side-peak ratio, off-diagonal: median %.4f, min %.4f, max %.4f" % (
    np.nanmedian(ratio), np.nanmin(ratio), np.nanmax(ratio)))
np.save(OUT / "step2b_ratio.npy", ratio)
