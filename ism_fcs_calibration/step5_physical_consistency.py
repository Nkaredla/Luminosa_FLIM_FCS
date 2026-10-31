"""Step 5: consequences of the fitted physical model.

1. Simulate the FCS data with the fitted optical parameters and fit the same Gaussian model
   -> must reproduce the measured Gaussian description (77 nm/pitch, sigma 130 nm).
2. Switch off saturation (and/or the depth) and predict the Gaussian-equivalent FCS shift.
3. Predict the in-focus bead APR shift at the coverslip for the fitted pitch, unsaturated and saturated.
"""
import sys, json, numpy as np
import step4b_fit_wave_parallel as w4
import gauss_model as gm
from fcs_common import OUT, d_atto655
from psf_mismatch import psf_rz
from fcs_wave_model import WaveFCS, saturate
from scipy.special import j1

D = d_atto655(24.0)
fit = json.load(open(OUT / "step4b_sat_T24.json"))
P = fit["params"]
phi_fit = np.array([P[k] for k in w4.NAMES])
z = np.load(OUT / "step3_gauss_T24.npz")
t = z["t"]; t = t[t <= 0.03]
iu = np.triu_indices(23)
w4._init(t, D)


def gauss_equiv(phi, label):
    G, integ = w4.unit_G(phi)
    r = gm.fit(t, G / 3.0, np.full_like(G, 1e-3), D, iu)
    A, T, tT, s, sig, sz = gm.unpack(r.x)
    sc, rot, mir, rms, _ = gm.similarity(s)
    out = dict(label=label, shift_per_pitch_nm=sc * 1e3, sigma_centre_nm=sig[11] * 1e3, K=sz[11] / sig[11],
               rate_min_max=float(integ.min() / integ.max()))
    print(f"{label:38s}: Gaussian-equivalent shift {sc * 1e3:5.1f} nm/pitch, sigma {sig[11] * 1e3:5.0f} nm, K {sz[11] / sig[11]:5.1f}, rates min/max {integ.min() / integ.max():.2f}")
    return out


def bead_shift(pitch, ap, sat_bead):
    m = WaveFCS(nxy=128, dx=0.01, nz=8, dz=0.01)
    rho = np.linspace(0, 1.0, 400)
    hex_ = psf_rz(0.640, rho, np.array([0.01]), 0.0)[0]
    hdet = psf_rz(0.690, rho, np.array([0.01]), 0.0)[0]
    Hex = saturate(np.interp(m.R, rho, hex_), sat_bead)
    Fdet = np.fft.fft2(np.interp(m.R, rho, hdet))
    kr = np.hypot(m.KX, m.KY); ra = 0.5 * ap * pitch
    with np.errstate(invalid="ignore", divide="ignore"):
        Fap = np.where(kr > 0, 2 * j1(kr * ra) / (kr * ra), 1.0)
    cen = []
    for cx, cy in m.xy_nom * pitch:
        im = Hex * np.clip(np.real(np.fft.ifft2(Fdet * Fap * np.exp(-1j * (m.KX * cx + m.KY * cy)))), 0, None)
        cen.append([np.sum(im * m.X) / im.sum(), np.sum(im * m.Y) / im.sum()])
    return gm.similarity(np.array(cen))[0] * 1e3


res = dict(fit=P)
res["fitted"] = gauss_equiv(phi_fit, "fitted model (as measured)")
nosat = phi_fit.copy(); nosat[7] = 0.0
res["no_saturation"] = gauss_equiv(nosat, "same optics, no saturation")
shallow = nosat.copy(); shallow[0] = 3.0
res["no_sat_3um"] = gauss_equiv(shallow, "no saturation, 3 um nominal depth")
res["bead_unsat"] = bead_shift(P["pitch"], P["ap"], 0.0)
res["bead_sat1"] = bead_shift(P["pitch"], P["ap"], 1.0)
res["bead_sat3"] = bead_shift(P["pitch"], P["ap"], 3.4)
print(f"bead APR at coverslip, fitted pitch {P['pitch'] * 1e3:.1f} nm: unsaturated {res['bead_unsat']:.1f}, "
      f"sat 1 {res['bead_sat1']:.1f}, sat 3.4 {res['bead_sat3']:.1f} nm/pitch")
json.dump(res, open(OUT / "step5_physical_consistency.json", "w"), indent=1)
