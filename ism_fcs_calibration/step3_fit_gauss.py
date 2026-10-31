"""Step 3: global 3-D Gaussian-MDF fit of all 23x23 FLCS correlation curves.

Model (pixel MDF_i = exp(-|r_xy - s_i|^2 / 2 sig_i^2 - z^2 / 2 sz_i^2)):

  G_ij(tau) = A / [(S_ij + 2 D tau) sqrt(Z_ij + 2 D tau)]
              * exp(-|s_i - s_j|^2 / (2 (S_ij + 2 D tau))) * (1 + T/(1-T) exp(-tau/tau_T))

  S_ij = sig_i^2 + sig_j^2,  Z_ij = sz_i^2 + sz_j^2,  A = 1 / ((2 pi)^1.5 c).

D is fixed to the Atto 655 reference, so every length comes out in um.
Gauge: the configuration {s_i} is fixed only up to translation/rotation/
mirror; the centre pixel is pinned at the origin and the right-hand
neighbour on the y = 0 axis, then the result is Procrustes-aligned to the
nominal honeycomb.
"""
import sys, json, numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fcs_common import OUT, NCH, d_atto655, load_segment_G, reject_segments, average_symmetric
from layout import honeycomb23

T_C = float(sys.argv[1]) if len(sys.argv) > 1 else 24.0
D = d_atto655(T_C)
TMIN, TMAX = 1.0e-6, 0.1
CEN = 11  # routing 20, centre of the honeycomb
RIGHT = 12

tau, G, W, z = load_segment_G(weight=1)
keep, slow, fast = reject_segments(tau, G)
Gm, Ge, nkeep = average_symmetric(G, keep)
print(f"T = {T_C} C, D = {D:.1f} um^2/s; segments kept {nkeep}/{G.shape[0]}")
sel = (tau >= TMIN) & (tau <= TMAX)
t = tau[sel]
iu = np.triu_indices(NCH)            # includes autocorrelations
Y = Gm[sel][:, iu[0], iu[1]]         # [nt, npair]
E = Ge[sel][:, iu[0], iu[1]]
E = np.maximum(E, 1e-4 * np.abs(Y).max())
pi_, pj_ = iu
xy0 = honeycomb23()


def unpack(p):
    A, T, tT = p[0], p[1], p[2]
    k = 3
    s = np.zeros((NCH, 2))
    free = [i for i in range(NCH) if i != CEN]
    for i in free:
        if i == RIGHT:
            s[i, 0] = p[k]; k += 1
        else:
            s[i] = p[k:k + 2]; k += 2
    sig = p[k:k + NCH]; k += NCH
    sz = p[k:k + NCH]; k += NCH
    return A, T, tT, s, sig, sz


def model(p, t=t):
    A, T, tT, s, sig, sz = unpack(p)
    Sxy = sig[pi_] ** 2 + sig[pj_] ** 2
    Sz = sz[pi_] ** 2 + sz[pj_] ** 2
    d2 = np.sum((s[pi_] - s[pj_]) ** 2, axis=1)
    a = Sxy[None] + 2 * D * t[:, None]
    b = Sz[None] + 2 * D * t[:, None]
    trip = 1 + T / (1 - T) * np.exp(-t / tT)
    return A / (a * np.sqrt(b)) * np.exp(-d2[None] / (2 * a)) * trip[:, None]


def resid(p):
    return ((model(p) - Y) / E).ravel()


# initial guess: nominal layout with 0.06 um shift per pitch, sig 0.12, sz 0.6
p0 = [0.01, 0.05, 2e-6]
s_init = (xy0 - xy0[CEN]) * 0.06
for i in range(NCH):
    if i == CEN: continue
    if i == RIGHT: p0.append(s_init[i, 0])
    else: p0 += list(s_init[i])
p0 += [0.12] * NCH + [0.6] * NCH
p0 = np.array(p0)
# amplitude from data
p0[0] = np.median(Y[0]) * (2 * 0.12 ** 2) * np.sqrt(2 * 0.6 ** 2)
lb = np.full(p0.size, -np.inf); ub = np.full(p0.size, np.inf)
lb[0], ub[0] = 0, np.inf
lb[1], ub[1] = 0, 0.5
lb[2], ub[2] = 1e-7, 2e-5
nsig = 3 + 2 * (NCH - 1) - 1
lb[nsig:] = 0.01; ub[nsig:] = 20
fit = least_squares(resid, p0, bounds=(lb, ub), x_scale="jac", max_nfev=20000, verbose=0)
p = fit.x
chi2r = np.sum(fit.fun ** 2) / (fit.fun.size - p.size)
J = fit.jac
cov = np.linalg.pinv(J.T @ J) * chi2r
perr = np.sqrt(np.diag(cov))
A, T, tT, s, sig, sz = unpack(p)
c = 1 / (A * (2 * np.pi) ** 1.5)
print(f"fit: chi2_red = {chi2r:.3f}, c = {c:.3f} /um^3 = {c / 0.6022:.3f} nM, T = {T:.3f}, tau_T = {tT * 1e6:.2f} us")

# Procrustes: rotation (and optional mirror) of s onto the scaled nominal layout
def procrustes(X, Y):
    Xc = X - X.mean(0); Yc = Y - Y.mean(0)
    U, S_, Vt = np.linalg.svd(Xc.T @ Yc)
    R = U @ Vt
    scale = S_.sum() / np.sum(Xc ** 2)
    return R, scale
best = None
for mirror in (False, True):
    X = s.copy()
    if mirror: X[:, 1] *= -1
    R, sc = procrustes(xy0, X)   # maps nominal -> measured
    pred = (xy0 - xy0.mean(0)) @ R * sc + X.mean(0)
    rms = np.sqrt(np.mean(np.sum((pred - X) ** 2, 1)))
    if best is None or rms < best[0]:
        best = (rms, mirror, R, sc)
rms, mirror, R, sc = best
X = s.copy()
if mirror: X[:, 1] *= -1
Xal = (X - X.mean(0)) @ R.T     # measured, in nominal orientation, centroid at 0
ang = np.degrees(np.arctan2(R[0, 1], R[0, 0]))
print(f"similarity fit to nominal honeycomb: shift per pitch = {sc * 1e3:.1f} nm, rotation {ang:.2f} deg, mirror={mirror}, rms residual {rms * 1e3:.1f} nm")

# full affine on aligned coordinates: Xal = xy0c @ Aff + b
xy0c = xy0 - xy0.mean(0)
Xa = np.hstack([xy0c, np.ones((NCH, 1))])
coef, *_ = np.linalg.lstsq(Xa, Xal, rcond=None)
Aff = coef[:2]
U_, Sv, Vt_ = np.linalg.svd(Aff)
res_aff = Xal - Xa @ coef
print(f"affine: principal scales {Sv[0] * 1e3:.1f} / {Sv[1] * 1e3:.1f} nm per pitch (anisotropy {Sv[0] / Sv[1] - 1:+.3%}), rms residual {np.sqrt(np.mean(np.sum(res_aff ** 2, 1))) * 1e3:.1f} nm")

out = dict(T_C=T_C, D=D, chi2r=chi2r, c_per_um3=c, c_nM=c / 0.6022, trip=T, tauT=tT,
           s_um=s.tolist(), s_aligned_um=Xal.tolist(), sig_um=sig.tolist(), sz_um=sz.tolist(),
           shift_per_pitch_um=sc, rot_deg=ang, mirror=bool(mirror), rms_similarity_um=rms,
           affine=Aff.tolist(), affine_scales_um=Sv.tolist(), nkeep=int(nkeep), keep=keep.tolist(),
           perr=perr.tolist())
json.dump(out, open(OUT / f"step3_gauss_T{T_C:.0f}.json", "w"), indent=1)
np.savez(OUT / f"step3_gauss_T{T_C:.0f}.npz", p=p, t=t, Y=Y, E=E, model=model(p), iu=np.array(iu),
         tau=tau, Gm=Gm, Ge=Ge, keep=keep, slow=slow, fast=fast)
for i in range(NCH):
    print(f"rout {i + 9:2d}: s_aligned = ({Xal[i, 0] * 1e3:7.1f}, {Xal[i, 1] * 1e3:7.1f}) nm   "
          f"sig = {sig[i] * 1e3:6.1f} nm   sz = {sz[i] * 1e3:7.1f} nm   K = {sz[i] / sig[i]:.2f}")
