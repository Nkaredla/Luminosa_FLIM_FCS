import numpy as np, json, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from layout import honeycomb23
z = np.load("out/step3_gauss_T24.npz"); t, Y, E, Mo, iu = z["t"], z["Y"], z["E"], z["model"], z["iu"]
xy = honeycomb23(); dist = np.linalg.norm(xy[iu[0]] - xy[iu[1]], axis=1)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
cen = 11
for j, col in [(11, "k"), (12, "C0"), (16, "C1"), (13, "C2"), (22, "C3"), (0, "C4")]:
    a, b = min(cen, j), max(cen, j)
    k = np.flatnonzero((iu[0] == a) & (iu[1] == b))[0]
    ax[0].errorbar(t, Y[:, k], E[:, k], fmt=".", ms=3, color=col, label=f"20 x {j + 9} (d={dist[k]:.2f})")
    ax[0].plot(t, Mo[:, k], "-", color=col, lw=1)
ax[0].set_xscale("log"); ax[0].set_xlabel("tau (s)"); ax[0].set_ylabel("G"); ax[0].legend(fontsize=7); ax[0].set_title("centre pixel vs others: data + Gaussian fit")
R = (Y - Mo) / E
for lo, hi, col in [(-0.1, 0.1, "k"), (0.9, 1.1, "C0"), (1.6, 2.1, "C1"), (2.5, 3.1, "C2"), (3.4, 5, "C3")]:
    m = (dist > lo) & (dist < hi)
    ax[1].plot(t, R[:, m].mean(1), color=col, label=f"pairs d {lo:.1f}-{hi:.1f}: n={m.sum()}")
ax[1].set_xscale("log"); ax[1].axhline(0, color="gray"); ax[1].set_title("mean normalised residual by pair distance"); ax[1].legend(fontsize=7)
# normalised curve shapes: auto (center) / cross (opposite corners) ratio
k1 = np.flatnonzero((iu[0] == 11) & (iu[1] == 11))[0]
for j in [12, 13, 22, 0]:
    a, b = min(cen, j), max(cen, j); k = np.flatnonzero((iu[0] == a) & (iu[1] == b))[0]
    ax[2].plot(t, Y[:, k] / Y[:, k1], ".", ms=3, label=f"G(20,{j + 9})/G(20,20)")
    ax[2].plot(t, Mo[:, k] / Mo[:, k1], "-", lw=1, color="gray")
k = np.flatnonzero((iu[0] == 0) & (iu[1] == 22))[0]
ax[2].plot(t, Y[:, k] / Y[:, k1], ".", ms=3, label="G(9,31)/G(20,20)"); ax[2].plot(t, Mo[:, k] / Mo[:, k1], "-", lw=1, color="gray")
ax[2].set_xscale("log"); ax[2].set_ylim(0.4, 1.1); ax[2].legend(fontsize=7); ax[2].set_title("cross/auto ratio (carries the absolute scale)")
fig.tight_layout(); fig.savefig("out/check_gauss_fit.png", dpi=120)
