"""Stability: Gaussian fit on each 5-min quarter of the acquisition."""
import numpy as np, json
from fcs_common import OUT, NCH, d_atto655, load_segment_G, average_symmetric
import gauss_model as gm
D = d_atto655(24.0)
tau, G, W, z = load_segment_G(weight=1)
sel = (tau >= 1e-6) & (tau <= 0.03); t = tau[sel]
iu = np.triu_indices(NCH)
rows = []
for q in range(4):
    keep = np.zeros(G.shape[0], bool); keep[q * 15:(q + 1) * 15] = True
    Gm, Ge, n = average_symmetric(G, keep)
    Y = Gm[sel][:, iu[0], iu[1]]; E = np.maximum(Ge[sel][:, iu[0], iu[1]], 1e-4)
    r = gm.fit(t, Y, E, D, iu)
    A, T, tT, s, sig, sz = gm.unpack(r.x)
    sc, rot, mir, rms, Xal = gm.similarity(s)
    chi = np.sum(r.fun ** 2) / (r.fun.size - r.x.size)
    c = 1 / (A * (2 * np.pi) ** 1.5)
    rate = W[keep].sum(0).sum() / (n * 20.0)
    print(f"quarter {q}: chi2r {chi:.3f}  shift/pitch {sc * 1e3:.1f} nm  rms {rms * 1e3:.1f} nm  sig(centre) {sig[11] * 1e3:.0f} nm  "
          f"median sig {np.median(sig) * 1e3:.0f}  sz {np.median(sz) * 1e3:.0f} nm  c {c / 0.6022:.2f} nM  T {T:.3f}  tauT {tT * 1e6:.1f} us  FLCS rate {rate / 1e3:.1f} kcps")
    rows.append(dict(q=q, shift=sc, rms=rms, sig=sig.tolist(), sz=sz.tolist(), c=c, s=Xal.tolist()))
json.dump(rows, open(OUT / "step3b_quarters.json", "w"))
