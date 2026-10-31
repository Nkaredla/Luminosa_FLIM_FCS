"""Vectorial PSF for an oil objective focusing through a glass/water interface.

Torok/Hell formulation, rotationally symmetric case (circular excitation or
an isotropic emitter in detection):

    I(rho, z) = |I0|^2 + 2 |I1|^2 + |I2|^2
    Im(rho, z) = int_0^alpha sqrt(cos t1) sin t1 a(t1) g_m(t1) J_m(k_perp rho)
                 exp(i k0 [n2 cos t2 z - n1 cos t1 d]) dt1

with g_0 = ts + tp cos t2, g_1 = tp sin t2, g_2 = ts - tp cos t2 and the
Fresnel coefficients of the interface. z is the distance from the
interface into the water, d the nominal focus depth (stage displacement
times nothing: the depth a matched-index focus would have). Rays beyond the
critical angle are evanescent and drop out a few hundred nm into water.
`a(t1)` is the pupil apodisation (Gaussian beam fill) and `w_extra` adds a
symmetric phase in radians as a function of the normalised pupil radius.
"""
import numpy as np
from scipy.special import jv


def psf_rz(lam, rho, z, d, NA=1.45, n1=1.518, n2=1.333, fill=None, w_extra=None, ntheta=6000):
    k0 = 2 * np.pi / lam
    alpha = np.arcsin(NA / n1)
    t1 = (np.arange(ntheta) + 0.5) / ntheta * alpha
    dt = alpha / ntheta
    s1, c1 = np.sin(t1), np.cos(t1)
    s2 = n1 * s1 / n2
    c2 = np.sqrt((1 - s2 ** 2).astype(complex))
    c2 = np.where(np.imag(c2) < 0, np.conj(c2), c2)
    ts = 2 * n1 * c1 / (n1 * c1 + n2 * c2)
    tp = 2 * n1 * c1 / (n2 * c1 + n1 * c2)
    amp = np.sqrt(c1) * s1 * dt
    r = s1 / np.sin(alpha)                          # normalised pupil radius
    if fill is not None:                            # Gaussian beam, 1/e^2 radius = fill * pupil radius
        amp = amp * np.exp(-(r / fill) ** 2)
    phase0 = -k0 * n1 * c1 * d
    if w_extra is not None:
        phase0 = phase0 + w_extra(r)
    g = [ts + tp * c2, tp * s2, ts - tp * c2]
    kperp = k0 * n1 * s1
    Z = np.exp(1j * (k0 * n2 * np.outer(z, c2) + phase0[None, :]))   # [nz, nt]
    out = 0
    for m, w in zip((0, 1, 2), (1, 2, 1)):
        B = jv(m, np.outer(kperp, rho))                                # [nt, nrho]
        Im = (Z * (amp * g[m])[None, :]) @ B
        out = out + w * np.abs(Im) ** 2
    return out                                                         # [nz, nrho]
