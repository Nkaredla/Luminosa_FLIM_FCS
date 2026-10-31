"""Shared helpers: correlation reduction, segment rejection, reference D for Atto 655."""
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent / "out"
NCH = 23


def water_viscosity_mPas(T_C):
    """Vogel equation for water viscosity (mPa s), good to ~0.5% at 0-100 C."""
    T = T_C + 273.15
    return np.exp(-3.7188 + 578.919 / (T - 137.546))


def d_atto655(T_C, D25=426.0):
    """Atto 655-COOH, 4.26e-6 cm^2/s at 25 C (2fFCS), scaled by T/eta(T). um^2/s."""
    return D25 * (T_C + 273.15) / 298.15 * water_viscosity_mPas(25.0) / water_viscosity_mPas(T_C)


def load_segment_G(weight=1):
    """Return tau (s), per-segment G[s, l, i, j], channel weights W[s, i]."""
    z = np.load(OUT / "step2_corr.npz")
    S = z["S"][:, weight]; W = z["W"][:, weight]
    lags, casc = z["lags"], z["casc"]
    M = z["seg_len"] / 2.0 ** casc
    tau = lags * 2.0 ** casc / float(z["sync"])
    norm = (M ** 2 / (M - lags))[None, :, None, None]
    G = S * norm / (W[:, None, :, None] * W[:, None, None, :]) - 1.0
    return tau, G, W, z


def reject_segments(tau, G, nmad=5.0):
    """Flag segments whose mean slow-time correlation (1-100 ms) is anomalous (aggregates, spikes)."""
    sl = (tau > 1e-3) & (tau < 0.1)
    fa = (tau > 5e-6) & (tau < 2e-5)
    slow = G[:, sl].mean(axis=(1, 2, 3))
    fast = G[:, fa].mean(axis=(1, 2, 3))
    keep = np.ones(G.shape[0], bool)
    for q in (slow, fast):
        med = np.median(q); mad = 1.4826 * np.median(np.abs(q - med))
        keep &= np.abs(q - med) < nmad * mad
    return keep, slow, fast


def average_symmetric(G, keep):
    """Mean over kept segments, symmetrised G_ij = G_ji; SEM from segment scatter."""
    Gk = G[keep]
    Gk = 0.5 * (Gk + np.swapaxes(Gk, 2, 3))
    n = Gk.shape[0]
    return Gk.mean(0), Gk.std(0, ddof=1) / np.sqrt(n), n
