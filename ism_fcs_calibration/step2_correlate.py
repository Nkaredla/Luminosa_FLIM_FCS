"""Step 2: FLCS-filtered multi-tau 23x23 correlation matrix, per time segment.

Every photon carries two weights: 1 (raw) and its FLCS fluorescence-filter
value, built per channel from the measured TCSPC histogram decomposed into
a fluorescence pattern and a flat pattern (dark counts + afterpulsing).
The FLCS-weighted correlation is the correlation of the fluorescence alone,
so it needs no background amplitude correction and has no afterpulsing peak
in the autocorrelations.

Output: out/step2_corr.npz with, per segment s, lag l and channel pair (i, j),
the lag sums S[s, l, i, j] = sum_t w_i(t) w_j(t + tau_l) and the per-segment
channel weight totals, from which G is formed in step 3.
"""
import sys, time, numpy as np
from pathlib import Path
from multiprocessing import Pool

CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else None
OUT = Path("out")
SEG_S = 20.0           # segment length (s)
SYNC = 40e6            # macro-time unit = 1 sync period = 25 ns
NSUB, NCASC = 10, 20   # multi-tau: max lag 10*(2^20-1)*25 ns = 0.26 s
ROUT = np.arange(9, 32)
NCH = 23


def flcs_filters(H, bg_window=(850, 1050)):
    """Fluorescence filter per channel from patterns [fluo, flat]."""
    nb = H.shape[1]
    F = np.zeros_like(H, dtype=float)
    bgfrac = np.zeros(H.shape[0])
    for k, h in enumerate(H.astype(float)):
        b = np.mean(h[bg_window[0]:bg_window[1]])
        pf = np.clip(h - b, 0, None); pf /= pf.sum()
        pb = np.full(nb, 1.0 / nb)
        Mt = np.stack([pf, pb])               # 2 x nb
        Dinv = 1.0 / np.maximum(h, 1.0)
        A = (Mt * Dinv) @ Mt.T
        f = np.linalg.solve(A, Mt * Dinv)      # 2 x nb
        F[k] = f[0]
        bgfrac[k] = b * nb / h.sum()
    return F, bgfrac


def lag_list():
    lags, casc = [], []
    shift, delta = 0, 1
    for j in range(NCASC):
        for _ in range(NSUB):
            shift += delta
            lags.append(shift // delta); casc.append(j)
        delta *= 2
    return np.array(lags), np.array(casc)


def correlate_segment(args):
    s, i0, i1, t0 = args
    ns = np.load(CACHE / "atto655_nsync.npy", mmap_mode="r")
    ch = np.load(CACHE / "atto655_chan.npy", mmap_mode="r")
    dt = np.load(CACHE / "atto655_dtime.npy", mmap_mode="r")
    F = np.load(OUT / "flcs_filters.npy")
    t = np.asarray(ns[i0:i1]) - t0
    c = np.asarray(ch[i0:i1]).astype(np.int64) - ROUT[0]
    d = np.asarray(dt[i0:i1]).astype(np.int64)
    keep = (c >= 0) & (c < NCH) & (d < F.shape[1])
    t, c, d = t[keep], c[keep], d[keep]
    W = np.stack([np.ones(t.size), F[c, d]], axis=1)     # raw, flcs
    Wtot = np.stack([np.bincount(c, weights=W[:, m], minlength=NCH) for m in range(2)])
    lags, casc = lag_list()
    S = np.zeros((2, lags.size, NCH, NCH))
    # merge identical (t, c) at level 0 (two photons of one channel in one sync period are impossible
    # with the SPAD dead time, but keep the code general)
    li = 0
    for j in range(NCASC):
        key = t * 32 + c
        order = np.argsort(key, kind="stable")
        key = key[order]; W = W[order]
        uk, inv = np.unique(key, return_inverse=True)
        if uk.size < key.size:
            W = np.stack([np.bincount(inv, weights=W[:, m]) for m in range(2)], axis=1)
        t = uk // 32; c = uk % 32
        brk = np.flatnonzero(np.diff(t)) + 1
        starts = np.concatenate([[0], brk]); counts = np.diff(np.concatenate([starts, [t.size]]))
        u = t[starts]
        single = bool(np.all(counts == 1))
        for _ in range(NSUB):
            L = lags[li]
            idx = np.searchsorted(u, u + L)
            ok = idx < u.size
            a = np.flatnonzero(ok); b = idx[ok]
            m = u[b] == u[a] + L
            a, b = a[m], b[m]
            if single:
                ia, ib = starts[a], starts[b]
            else:
                na, nb_ = counts[a], counts[b]
                tot = na * nb_
                P = int(tot.sum())
                pid = np.repeat(np.arange(a.size), tot)
                k = np.arange(P) - np.repeat(np.cumsum(tot) - tot, tot)
                ia = starts[a][pid] + k // nb_[pid]
                ib = starts[b][pid] + k % nb_[pid]
            code = c[ia] * NCH + c[ib]
            for mm in range(2):
                S[mm, li] = np.bincount(code, weights=W[ia, mm] * W[ib, mm], minlength=NCH * NCH).reshape(NCH, NCH)
            li += 1
        t = t >> 1
    return s, S, Wtot


def main():
    H = np.load(OUT / "step1_tcspc.npz")["H"]
    F, bgfrac = flcs_filters(H)
    np.save(OUT / "flcs_filters.npy", F)
    print("flat-background fraction per channel:", np.round(bgfrac, 4))
    ns = np.load(CACHE / "atto655_nsync.npy", mmap_mode="r")
    Tsync = int(ns[-1]) + 1
    seg_len = int(SEG_S * SYNC)
    nseg = Tsync // seg_len
    edges = np.searchsorted(ns, np.arange(nseg + 1) * seg_len)
    jobs = [(s, int(edges[s]), int(edges[s + 1]), s * seg_len) for s in range(nseg)]
    lags, casc = lag_list()
    S = np.zeros((nseg, 2, lags.size, NCH, NCH)); Wt = np.zeros((nseg, 2, NCH))
    t0 = time.time()
    with Pool(10) as pool:
        for s, Ss, Ws in pool.imap_unordered(correlate_segment, jobs):
            S[s] = Ss; Wt[s] = Ws
            print(f"segment {s:3d} done  ({time.time() - t0:.0f} s)", flush=True)
    np.savez_compressed(OUT / "step2_corr.npz", S=S, W=Wt, lags=lags, casc=casc,
                        seg_len=seg_len, sync=SYNC, rout=ROUT, bgfrac=bgfrac)


if __name__ == "__main__":
    main()
