"""Wave-optics forward model of the 23x23 ISM-FCS correlation matrix in solution.

MDF_i(r) = h_ex(r) * [h_det (x)_xy aperture](r_xy - c_i, z)

with h_ex, h_det from psf_mismatch.psf_rz (oil objective -> water at
nominal depth d), c_i = pixel centre projected into the sample plane
(nominal honeycomb * pitch + array offset), a disc aperture of diameter
`ap * pitch`. For free 3-D diffusion,

    G_ij(tau) = sum_k Re[M_i(k) M_j*(k)] exp(-D k^2 tau)
                / (c * Ntot * dV * M_i(0) M_j(0))      (k = 0 excluded)

with M_i the DFT of MDF_i on a periodic grid. The lateral |k|^2 is binned
finely so the tau loop is cheap.
"""
import numpy as np
from psf_mismatch import psf_rz
from layout import honeycomb23


def saturate(h, sat):
    """Pulsed excitation: excitation probability per pulse 1 - exp(-sat * h / max(h)).

    `sat` is the mean number of photons absorbed per pulse by a molecule at
    the excitation maximum (pulse << lifetime). sat -> 0 is the linear MDF.
    """
    if sat <= 1e-4:
        return h
    return -np.expm1(-sat * h / h.max()) / sat


class WaveFCS:
    def __init__(self, nxy=96, dx=0.06, nz=240, dz=0.08, nkbin=160, lam_ex=0.640, lam_em=0.690,
                 NA=1.45, n1=1.518, n2=1.333):
        self.nxy, self.dx, self.nz, self.dz = nxy, dx, nz, dz
        self.lam_ex, self.lam_em, self.NA, self.n1, self.n2 = lam_ex, lam_em, NA, n1, n2
        x = (np.arange(nxy) - nxy // 2) * dx
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")
        self.R = np.hypot(self.X, self.Y)
        self.zrel = (np.arange(nz) - nz // 2) * dz
        kx = 2 * np.pi * np.fft.fftfreq(nxy, dx)
        self.KX, self.KY = np.meshgrid(kx, kx, indexing="ij")
        kp2 = (self.KX ** 2 + self.KY ** 2).ravel()
        self.kz = 2 * np.pi * np.fft.rfftfreq(nz, dz)
        kmax2 = kp2.max()
        # bins equally spaced in |k_perp|
        edges = np.linspace(0, np.sqrt(kmax2) * 1.0001, nkbin + 1)
        b = np.digitize(np.sqrt(kp2), edges) - 1
        self.order = np.argsort(b, kind="stable")
        self.bstart = np.searchsorted(b[self.order], np.arange(nkbin + 1))
        self.kp2_bin = np.array([kp2[self.order[self.bstart[i]:self.bstart[i + 1]]].mean()
                                 if self.bstart[i + 1] > self.bstart[i] else 0.0 for i in range(nkbin)])
        self.nkbin = nkbin
        self.wz = np.full(self.kz.size, 2.0); self.wz[0] = 1.0
        if nz % 2 == 0: self.wz[-1] = 1.0
        self.rho = np.linspace(0, self.R.max() * 1.001, 400)
        self._psf_cache = {}
        self.xy_nom = honeycomb23()

    def psfs(self, d, fill):
        key = (round(d, 5), None if fill is None else round(fill, 5))
        if key in self._psf_cache:
            return self._psf_cache[key]
        # locate the excitation focus first (on-axis scan), then sample around it
        zs = np.linspace(max(0.05, 0.6 * d - 3), d + 3, 600)
        ax = psf_rz(self.lam_ex, np.array([0.0]), zs, d, self.NA, self.n1, self.n2, fill)[:, 0]
        zf = zs[np.argmax(ax)]
        z = zf + self.zrel
        z = np.where(z > 0.02, z, 0.02)
        hex_ = psf_rz(self.lam_ex, self.rho, z, d, self.NA, self.n1, self.n2, fill)
        hdet = psf_rz(self.lam_em, self.rho, z, d, self.NA, self.n1, self.n2, None)
        hex_ /= hex_.max(); hdet /= hdet.max()
        self._psf_cache = {key: (hex_, hdet, zf)}
        return hex_, hdet, zf

    def _to_grid(self, h):
        """[nz, nrho] -> [nxy, nxy, nz] by linear interpolation in rho."""
        r = self.R.ravel()
        idx = np.clip(np.searchsorted(self.rho, r) - 1, 0, self.rho.size - 2)
        w = (r - self.rho[idx]) / (self.rho[idx + 1] - self.rho[idx])
        g = h[:, idx] * (1 - w) + h[:, idx + 1] * w            # [nz, nxy*nxy]
        return g.T.reshape(self.nxy, self.nxy, self.nz)

    def mdfs(self, d, pitch, offset=(0.0, 0.0), ap=0.75, fill=None, xy=None, sat=0.0):
        hex_, hdet, zf = self.psfs(d, fill)
        Hex = saturate(self._to_grid(hex_), sat)
        Fdet = np.fft.fft2(self._to_grid(hdet), axes=(0, 1))
        kr = np.hypot(self.KX, self.KY)
        ra = 0.5 * ap * pitch
        with np.errstate(invalid="ignore", divide="ignore"):
            from scipy.special import j1
            Fap = np.where(kr > 0, 2 * j1(kr * ra) / (kr * ra), 1.0)
        xy = self.xy_nom if xy is None else xy
        cen = xy * pitch + np.asarray(offset)[None, :]
        M = []
        for cx, cy in cen:
            ph = np.exp(-1j * (self.KX * cx + self.KY * cy)) * Fap
            hd = np.real(np.fft.ifft2(Fdet * ph[:, :, None], axes=(0, 1)))
            M.append(Hex * np.clip(hd, 0, None))
        return np.array(M), zf

    def spectra(self, M):
        """Binned cross-spectra P[i, j, kbin, kz] and integrals."""
        F = np.fft.rfftn(M, axes=(1, 2, 3))                     # [n, nxy, nxy, nkz]
        n = F.shape[0]
        F = F.reshape(n, -1, self.kz.size)[:, self.order]
        P = np.zeros((n, n, self.nkbin, self.kz.size))
        for b in range(self.nkbin):
            s0, s1 = self.bstart[b], self.bstart[b + 1]
            if s1 > s0:
                Fb = F[:, s0:s1].transpose(2, 0, 1)                 # [nkz, n, nb]
                P[:, :, b] = np.real(Fb @ np.conj(Fb).transpose(0, 2, 1)).transpose(1, 2, 0)
        P *= self.wz[None, None, None, :]
        P[:, :, 0, 0] -= np.real(F[:, 0, 0][:, None] * np.conj(F[:, 0, 0][None, :]))  # drop k = 0
        integ = np.real(F[:, 0, 0]).copy()   # copy: a view would keep all of F alive
        return P, integ

    def correlations(self, P, integ, tau, D):
        """G_ij(tau) per unit concentration (c = 1 /um^3). Returns [ntau, n, n]."""
        dV = self.dx * self.dx * self.dz
        ntot = self.nxy * self.nxy * self.nz
        Ep = np.exp(-D * np.outer(tau, self.kp2_bin))              # [nt, nkbin]
        Ez = np.exp(-D * np.outer(tau, self.kz ** 2))              # [nt, nkz]
        G = np.einsum("ijbz,tb,tz->tij", P, Ep, Ez, optimize=True)
        return G / (ntot * dV * np.outer(integ, integ))[None]

    def gaussian_beam_fill(self):
        return None
