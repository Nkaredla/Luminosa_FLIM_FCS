"""Step 4b: same physical fit as step4_fit_wave.py, with the Jacobian columns computed in parallel.

Parameters: d, pitch, x0, y0, ap, fill (fixed), dzd, sat | A, T, tauT. Logs every evaluation.
"""
import sys, json, time, numpy as np
from multiprocessing import Pool
from scipy.optimize import least_squares

NAMES = ["d", "pitch", "x0", "y0", "ap", "fill", "dzd", "sat"]
NOPT = len(NAMES)
STEP = np.array([0.15, 0.002, 0.004, 0.004, 0.03, 0.05, 0.05, 0.1])
FIXED = {5}
_W = None


def _init(tsel, D):
    global _W, _t, _D, _iu
    from fcs_wave_model import WaveFCS
    _W = WaveFCS(); _t = tsel; _D = D; _iu = np.triu_indices(23)


def unit_G(phi):
    """Expensive part, run inside a worker. Returns (Gu [nt, npair], integ)."""
    from psf_mismatch import psf_rz
    from fcs_wave_model import saturate
    from scipy.special import j1
    import scipy.fft as sfft
    m = _W
    d, pitch, x0, y0, ap, fill, dzd, sat = phi
    key = ("psf", round(d, 5), round(dzd, 5))
    if key not in m._psf_cache:
        zs = np.linspace(max(0.05, 0.6 * d - 3), d + 3, 600)
        ax = psf_rz(m.lam_ex, np.array([0.0]), zs, d, m.NA, m.n1, m.n2, None)[:, 0]
        zf = zs[np.argmax(ax)]
        z = np.clip(zf + m.zrel, 0.02, None)
        hex_ = psf_rz(m.lam_ex, m.rho, z, d, m.NA, m.n1, m.n2, None)
        hdet = psf_rz(m.lam_em, m.rho, np.clip(z - dzd, 0.02, None), d, m.NA, m.n1, m.n2, None)
        if len(m._psf_cache) > 3:
            m._psf_cache.pop(next(iter(m._psf_cache)))
        m._psf_cache[key] = (hex_ / hex_.max(), hdet / hdet.max(), zf)
    hex_, hdet, zf = m._psf_cache[key]
    Hex = saturate(m._to_grid(hex_), sat).astype(np.float32)
    Fdet = sfft.fft2(m._to_grid(hdet).astype(np.float32), axes=(0, 1))
    kr = np.hypot(m.KX, m.KY); ra = 0.5 * ap * pitch
    with np.errstate(invalid="ignore", divide="ignore"):
        Fap = np.where(kr > 0, 2 * j1(kr * ra) / (kr * ra), 1.0)
    cen = m.xy_nom * pitch + np.array([x0, y0])[None, :]
    F = []
    for cx, cy in cen:
        ph = (np.exp(-1j * (m.KX * cx + m.KY * cy)) * Fap).astype(np.complex64)
        hd = sfft.ifft2(Fdet * ph[:, :, None], axes=(0, 1)).real
        F.append(sfft.rfftn(Hex * np.clip(hd, 0, None), axes=(0, 1, 2)))
    F = np.array(F)
    n = F.shape[0]
    F = F.reshape(n, -1, m.kz.size)[:, m.order]
    P = np.zeros((n, n, m.nkbin, m.kz.size))
    for b in range(m.nkbin):
        s0, s1 = m.bstart[b], m.bstart[b + 1]
        if s1 > s0:
            Fb = F[:, s0:s1].transpose(2, 0, 1).astype(np.complex128)
            P[:, :, b] = np.real(Fb @ np.conj(Fb).transpose(0, 2, 1)).transpose(1, 2, 0)
    P *= m.wz[None, None, None, :]
    F0 = F[:, 0, 0].astype(np.complex128)
    P[:, :, 0, 0] -= np.real(F0[:, None] * np.conj(F0[None, :]))
    integ = np.real(F0)
    G = m.correlations(P, integ, _t, _D)[:, _iu[0], _iu[1]]
    return G, integ


def main():
    from fcs_common import OUT, NCH, d_atto655
    T_C = float(sys.argv[1]) if len(sys.argv) > 1 else 24.0
    TMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 3e-3
    TAG = sys.argv[3] if len(sys.argv) > 3 else "par"
    D = d_atto655(T_C)
    g = np.load(OUT / f"step3_gauss_T{T_C:.0f}.npz")
    tau, Gm, Ge = g["tau"], g["Gm"], g["Ge"]
    sel = (tau >= 1e-6) & (tau <= TMAX); t = tau[sel]
    iu = np.triu_indices(NCH)
    Y = Gm[sel][:, iu[0], iu[1]]
    E = np.maximum(Ge[sel][:, iu[0], iu[1]], 1e-4 * np.abs(Y).max())
    pool = Pool(7, initializer=_init, initargs=(t, D))
    cache = {}
    log = open(OUT / f"step4b_{TAG}_progress.log", "w")

    def many(phis):
        keys = [tuple(np.round(p, 6)) for p in phis]
        todo = [(k, p) for k, p in zip(keys, phis) if k not in cache]
        if todo:
            for (k, _), res in zip(todo, pool.map(unit_G, [p for _, p in todo])):
                cache[k] = res
        return [cache[k] for k in keys]

    def full(p, Gu):
        A, T, tT = p[NOPT:]
        return A * Gu * (1 + T / (1 - T) * np.exp(-t / tT))[:, None]

    def resid(p):
        Gu, _ = many([p[:NOPT]])[0]
        r = ((full(p, Gu) - Y) / E).ravel()
        log.write(f"{time.strftime('%H:%M:%S')} chi2r {np.sum(r ** 2) / (r.size - 10):.4f} " +
                  " ".join(f"{n}={v:.4g}" for n, v in zip(NAMES + ['A', 'T', 'tauT'], p)) + "\n"); log.flush()
        return r

    def jac(p):
        base = p[:NOPT]
        phis = [base]
        cols = []
        for k in range(NOPT):
            if k in FIXED: continue
            q = base.copy(); q[k] += STEP[k]; phis.append(q); cols.append(k)
        res = many(phis)
        Gu0 = res[0][0]
        r0 = ((full(p, Gu0) - Y) / E).ravel()
        J = np.zeros((r0.size, p.size))
        for (Gu, _), k in zip(res[1:], cols):
            J[:, k] = (((full(p, Gu) - Y) / E).ravel() - r0) / STEP[k]
        for k in range(NOPT, p.size):
            h = max(1e-7, 1e-4 * abs(p[k])); q = p.copy(); q[k] += h
            J[:, k] = (((full(q, Gu0) - Y) / E).ravel() - r0) / h
        return J

    p0 = np.array([18.27, 0.124, 0.0, 0.0, 0.75, 5.0, 0.0, 3.0, 0.1, 0.05, 5e-6])
    if len(sys.argv) > 4:
        p0[:NOPT] = json.loads(sys.argv[4])
    Gu0, _ = many([p0[:NOPT]])[0]
    ic = np.flatnonzero((iu[0] == 11) & (iu[1] == 11))[0]
    p0[NOPT] = Y[2, ic] / Gu0[2, ic]
    lb = [3.0, 0.06, -0.15, -0.15, 0.2, 4.99, -1.5, 0.0, 0, 0, 2e-7]
    ub = [35.0, 0.25, 0.15, 0.15, 1.5, 5.01, 1.5, 30.0, np.inf, 0.4, 3e-5]
    t0 = time.time()
    fit = least_squares(resid, p0, jac=jac, bounds=(lb, ub), x_scale=np.r_[STEP * 10, p0[NOPT], 0.05, 1e-6],
                        max_nfev=40, verbose=2, xtol=1e-6, ftol=1e-6)
    p = fit.x
    chi2r = np.sum(fit.fun ** 2) / (fit.fun.size - 10)
    free = [k for k in range(p.size) if k not in FIXED]
    J = fit.jac; cov = np.zeros((p.size, p.size))
    cov[np.ix_(free, free)] = np.linalg.pinv(J[:, free].T @ J[:, free]) * chi2r
    err = np.sqrt(np.diag(cov))
    Gu, integ = many([p[:NOPT]])[0]
    res = dict(T_C=T_C, D=D, tmax=TMAX, chi2r=chi2r, params=dict(zip(NAMES + ["A", "T", "tauT"], p.tolist())),
               errors=dict(zip(NAMES + ["A", "T", "tauT"], err.tolist())), c_nM=1 / p[NOPT] / 0.6022,
               rel_rates=(integ / integ.max()).tolist(), runtime_s=time.time() - t0, nfev=fit.nfev)
    print(json.dumps(res, indent=1))
    json.dump(res, open(OUT / f"step4b_{TAG}_T{T_C:.0f}.json", "w"), indent=1)
    np.savez(OUT / f"step4b_{TAG}_T{T_C:.0f}.npz", p=p, t=t, Y=Y, E=E, model=full(p, Gu), iu=np.array(iu), cov=cov)
    pool.close()


if __name__ == "__main__":
    main()
