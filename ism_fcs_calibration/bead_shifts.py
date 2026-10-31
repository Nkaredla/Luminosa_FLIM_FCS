"""Collect the bead-image APR shift tables (run_ism_reconstruction_from_ptu.m output) in nm."""
import numpy as np, csv
from pathlib import Path
from ptu_io import read_header

BEAD_DIR = Path(r"D:\Luminosa\Data\270226\ISM\GattaBeads_25\good")


def load_bead_shifts():
    """Return dict name -> (pix_nm, shift_to_apply_nm [23 x 2] as (dx, dy), channel IDs)."""
    out = {}
    for k in range(6):
        f = BEAD_DIR / f"RawImage{k}_detector_map.csv"
        tags, _ = read_header(BEAD_DIR / f"RawImage{k}.ptu")
        pix_nm = tags["ImgHdr_PixResol"] * 1e3
        rows = list(csv.DictReader(open(f)))
        ids = np.array([int(r["ChannelID"]) for r in rows])
        sh = np.array([[float(r["ShiftToCenterDX"]), float(r["ShiftToCenterDY"])] for r in rows]) * pix_nm
        out[f"RawImage{k}"] = (pix_nm, sh, ids, tags["ImgHdr_PixX"], tags["ImgHdr_Z0"])
    return out


if __name__ == "__main__":
    from layout import honeycomb23
    xy = honeycomb23()
    for name, (pix, sh, ids, npx, z0) in load_bead_shifts().items():
        # least-squares similarity scale of -shift (image displacement) vs nominal layout (allow rotation+mirror)
        X = -sh; Xc = X - X.mean(0); N = xy - xy.mean(0)
        best = None
        for mir in (1, -1):
            Xm = Xc * [1, mir]
            U, S, Vt = np.linalg.svd(N.T @ Xm); R = U @ Vt; sc = S.sum() / np.sum(N ** 2)
            r = np.sqrt(np.mean(np.sum((N @ R * sc - Xm) ** 2, 1)))
            if best is None or r < best[0]: best = (r, mir, sc, np.degrees(np.arctan2(R[0, 1], R[0, 0])))
        print(f"{name}: pix {pix:.2f} nm, {npx} px, Z0 {z0:.2f}: shift/pitch {best[2]:.1f} nm, rot {best[3]:.1f} deg, mirror {best[1]}, rms {best[0]:.1f} nm")
