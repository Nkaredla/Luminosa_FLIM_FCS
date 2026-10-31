"""3-D Gaussian-MDF ISM-FCS model and fit, reusable for data and simulations (see step3)."""
import numpy as np
from scipy.optimize import least_squares
from layout import honeycomb23

NCH = 23
CEN, RIGHT = 11, 12


def unpack(p):
    A, T, tT = p[0], p[1], p[2]
    k = 3
    s = np.zeros((NCH, 2))
    for i in range(NCH):
        if i == CEN:
            continue
        if i == RIGHT:
            s[i, 0] = p[k]; k += 1
        else:
            s[i] = p[k:k + 2]; k += 2
    sig = p[k:k + NCH]; k += NCH
    sz = p[k:k + NCH]
    return A, T, tT, s, sig, sz


def model(p, t, D, iu):
    A, T, tT, s, sig, sz = unpack(p)
    pi_, pj_ = iu
    a = (sig[pi_] ** 2 + sig[pj_] ** 2)[None] + 2 * D * t[:, None]
    b = (sz[pi_] ** 2 + sz[pj_] ** 2)[None] + 2 * D * t[:, None]
    d2 = np.sum((s[pi_] - s[pj_]) ** 2, axis=1)
    trip = 1 + T / (1 - T) * np.exp(-t / tT)
    return A / (a * np.sqrt(b)) * np.exp(-d2[None] / (2 * a)) * trip[:, None]


def fit(t, Y, E, D, iu, shift0=0.06, sig0=0.12, sz0=0.6, fix_trip=None):
    xy0 = honeycomb23()
    p0 = [0.01, 0.05, 2e-6]
    s_init = (xy0 - xy0[CEN]) * shift0
    for i in range(NCH):
        if i == CEN: continue
        p0 += [s_init[i, 0]] if i == RIGHT else list(s_init[i])
    p0 += [sig0] * NCH + [sz0] * NCH
    p0 = np.array(p0)
    p0[0] = np.median(Y[0]) * (2 * sig0 ** 2) * np.sqrt(2 * sz0 ** 2)
    lb = np.full(p0.size, -np.inf); ub = np.full(p0.size, np.inf)
    lb[0] = 0; lb[1], ub[1] = 0, 0.5; lb[2], ub[2] = 1e-7, 2e-5
    nsig = 3 + 2 * (NCH - 1) - 1
    lb[nsig:] = 0.01; ub[nsig:] = 20
    if fix_trip is not None:
        p0[1] = fix_trip; lb[1] = fix_trip - 1e-9; ub[1] = fix_trip + 1e-9
    r = least_squares(lambda p: ((model(p, t, D, iu) - Y) / E).ravel(), p0, bounds=(lb, ub),
                      x_scale="jac", max_nfev=20000)
    return r


def similarity(s):
    """Shift per pitch (um), rotation (deg), mirror flag, rms residual (um), aligned coordinates."""
    xy0 = honeycomb23()
    best = None
    for mirror in (False, True):
        X = s.copy()
        if mirror: X[:, 1] *= -1
        Xc = X - X.mean(0); N = xy0 - xy0.mean(0)
        U, S_, Vt = np.linalg.svd(N.T @ Xc); R = U @ Vt; sc = S_.sum() / np.sum(N ** 2)
        rms = np.sqrt(np.mean(np.sum((N @ R * sc - Xc) ** 2, 1)))
        if best is None or rms < best[0]:
            best = (rms, mirror, R, sc, Xc @ R.T)
    rms, mirror, R, sc, Xal = best
    return sc, np.degrees(np.arctan2(R[0, 1], R[0, 0])), mirror, rms, Xal
