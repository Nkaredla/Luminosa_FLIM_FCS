"""Gaussian MDF model with a common lateral ellipticity (e, theta): is the 3-D MDF elliptical on average?"""
import numpy as np, json
from scipy.optimize import least_squares
import gauss_model as gm
from fcs_common import OUT, d_atto655
z = np.load(OUT / "step3_gauss_T24.npz")
t, Y, E, iu, p_iso = z["t"], z["Y"], z["E"], z["iu"], z["p"]
D = d_atto655(24.0)
pi_, pj_ = iu


def model(q):
    p, e, th = q[:-2], q[-2], q[-1]
    A, T, tT, s, sig, sz = gm.unpack(p)
    c, sn = np.cos(th), np.sin(th)
    R = np.array([[c, -sn], [sn, c]])
    base = (sig[pi_] ** 2 + sig[pj_] ** 2)
    lu = base[None] * (1 + e) + 2 * D * t[:, None]     # along u
    lv = base[None] * (1 - e) + 2 * D * t[:, None]
    b = (sz[pi_] ** 2 + sz[pj_] ** 2)[None] + 2 * D * t[:, None]
    ds = (s[pi_] - s[pj_]) @ R                          # components along u, v
    arg = ds[:, 0][None] ** 2 / lu + ds[:, 1][None] ** 2 / lv
    trip = 1 + T / (1 - T) * np.exp(-t / tT)
    return A / (np.sqrt(lu * lv) * np.sqrt(b)) * np.exp(-arg / 2) * trip[:, None]


r0 = np.sum(((gm.model(p_iso, t, D, iu) - Y) / E) ** 2)
best = None
for th0 in np.linspace(0, np.pi, 4, endpoint=False):
    q0 = np.r_[p_iso, 0.02, th0]
    lb = np.r_[np.full(p_iso.size, -np.inf), -0.5, -10]; ub = np.r_[np.full(p_iso.size, np.inf), 0.5, 10]
    lb[0] = 0; lb[1], ub[1] = 0, 0.5; lb[2], ub[2] = 1e-7, 2e-5
    nsig = 3 + 2 * 22 - 1; lb[nsig:-2] = 0.01; ub[nsig:-2] = 20
    q0 = np.clip(q0, lb + 1e-9, ub - 1e-9)
    f = least_squares(lambda q: ((model(q) - Y) / E).ravel(), q0, bounds=(lb, ub), x_scale="jac", max_nfev=5000)
    if best is None or f.cost < best.cost: best = f
q = best.x; r1 = 2 * best.cost
e, th = q[-2], np.degrees(q[-1]) % 180
if e < 0: e, th = -e, (th + 90) % 180
J = best.jac; cov = np.linalg.pinv(J.T @ J) * r1 / (best.fun.size - q.size)
print(f"isotropic chi2 {r0:.1f}  ->  elliptical chi2 {r1:.1f}  (delta {r0 - r1:.1f} for 2 extra params; N = {best.fun.size})")
print(f"MDF lateral variance ellipticity e = {e:.4f} +/- {np.sqrt(cov[-2, -2]):.4f}  (sigma ratio {np.sqrt((1 + e) / (1 - e)):.4f}), major axis {th:.0f} deg in the chip frame")
json.dump(dict(dchi2=r0 - r1, e=e, e_err=float(np.sqrt(cov[-2, -2])), theta_deg=th), open(OUT / "step3c_anisotropy.json", "w"))
