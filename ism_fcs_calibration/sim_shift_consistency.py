"""Does one pitch explain both the bead APR shift (in focus, coverslip) and the FCS Gaussian shift (18 um deep)?

1. In-focus bead image per pixel at the coverslip: I_i(x) = h_ex(x, z0) * [h_det (x) ap](x - c_i, z0);
   APR shift = position of the cross-correlation peak against the centre pixel image.
2. Solution FCS at depth d: noiseless G_ij(tau) from the wave model, then the same Gaussian fit as the data.
Both reported as shift per detector pitch / pitch.
"""
import sys, json, numpy as np
from psf_mismatch import psf_rz
from fcs_wave_model import WaveFCS, saturate
from fcs_common import OUT, d_atto655
import gauss_model as gm
from scipy.special import j1

D = d_atto655(24.0)
pitch = float(sys.argv[1]) if len(sys.argv) > 1 else 0.121
fill = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
ap = float(sys.argv[3]) if len(sys.argv) > 3 else 0.75
depth = float(sys.argv[4]) if len(sys.argv) > 4 else 18.27
sat = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
sat_bead = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
fl = None if fill > 4.9 else fill

# ---- 1. bead at the coverslip (focus at the interface, emitter 10 nm into water)
m = WaveFCS(nxy=128, dx=0.01, nz=8, dz=0.01)
zb = np.array([0.01])
rho = np.linspace(0, 1.0, 400)
hex_ = psf_rz(0.640, rho, zb, 0.0, fill=fl)[0]
hdet = psf_rz(0.690, rho, zb, 0.0)[0]
Rg = m.R
Hex = saturate(np.interp(Rg, rho, hex_), sat_bead); Hdet = np.interp(Rg, rho, hdet)
Fdet = np.fft.fft2(Hdet)
kr = np.hypot(m.KX, m.KY); ra = 0.5 * ap * pitch
with np.errstate(invalid="ignore", divide="ignore"):
    Fap = np.where(kr > 0, 2 * j1(kr * ra) / (kr * ra), 1.0)
ims = []
for cx, cy in m.xy_nom * pitch:
    hd = np.real(np.fft.ifft2(Fdet * Fap * np.exp(-1j * (m.KX * cx + m.KY * cy))))
    ims.append(Hex * np.clip(hd, 0, None))
ims = np.array(ims)
# image centroid relative to the centre-pixel image centroid (equivalent to the APR peak for these smooth images)
X, Y = m.X, m.Y
cen = np.array([[np.sum(im * X) / im.sum(), np.sum(im * Y) / im.sum()] for im in ims])
sc_b, rot_b, mir_b, rms_b, _ = gm.similarity(cen)
print(f"bead in focus: pitch {pitch * 1e3:.0f} nm, fill {fill}, ap {ap}: APR shift per pitch {sc_b * 1e3:.1f} nm (alpha = {sc_b / pitch:.3f})")

# ---- 2. FCS at depth, Gaussian-equivalent shift
w = WaveFCS()
M, zf = w.mdfs(depth, pitch, ap=ap, fill=fl, sat=sat)
P, integ = w.spectra(M)
g = np.load(OUT / "step3_gauss_T24.npz")
t = g["t"]; t = t[t <= 0.03]
iu = np.triu_indices(23)
Gs = w.correlations(P, integ, t, D)[:, iu[0], iu[1]] / 3.0     # c = 3 /um^3
E = np.full_like(Gs, 1e-3)
r = gm.fit(t, Gs, E, D, iu)
A, T, tT, s, sig, sz = gm.unpack(r.x)
sc_f, rot_f, mir_f, rms_f, _ = gm.similarity(s)
print(f"FCS at d={depth}, sat={sat}: Gaussian-equivalent shift per pitch {sc_f * 1e3:.1f} nm (beta = {sc_f / pitch:.3f}), "
      f"sig centre {sig[11] * 1e3:.0f} nm, sz centre {sz[11] * 1e3:.0f} nm, K {sz[11] / sig[11]:.1f}, T {T:.3f}, tauT {tT * 1e6:.1f} us, "
      f"rates min/max {integ.min() / integ.max():.2f}")
print(f"ratio FCS/bead shift = {sc_f / sc_b:.3f}")
