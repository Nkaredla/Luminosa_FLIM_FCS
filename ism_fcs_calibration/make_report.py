"""Calibration table (CSV/JSON) and summary figures for the PDA-23 Atto 655 FCS calibration."""
import json, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection
from layout import honeycomb23
from fcs_common import OUT, d_atto655

xy = honeycomb23()
ROUT = np.arange(9, 32)
g = json.load(open(OUT / "step3_gauss_T24.json"))
g25 = json.load(open(OUT / "step3_gauss_T25.json"))
st = np.load(OUT / "step1_tcspc.npz")
eps = np.load(OUT / "crosstalk_eps.npy")
S = np.array(g["s_aligned_um"]) * 1e3
sig = np.array(g["sig_um"]) * 1e3
sz = np.array(g["sz_um"]) * 1e3
perr = np.array(g["perr"]) * 1e3
nsig = 3 + 2 * 22 - 1
sig_err = perr[nsig:nsig + 23]
# excitation axis from the sigma^2 paraboloid (see analysis log)
Am = np.c_[np.ones(23), S[:, 0], S[:, 1], (S ** 2).sum(1)]
co, *_ = np.linalg.lstsq(Am, sig ** 2, rcond=None)
axis = np.array([-co[1] / (2 * co[3]), -co[2] / (2 * co[3])])
offs_ps = (st["peak"] - np.median(st["peak"])) * 1e3
dist = np.linalg.norm(xy[:, None] - xy[None], axis=2)
nn_eps = np.array([np.nanmean(eps[i][(dist[i] > 0.9) & (dist[i] < 1.1)]) for i in range(23)])
nn_sum = np.array([np.nansum(eps[i][(dist[i] > 0.1)]) for i in range(23)])
beads = np.load(OUT / "bead_pattern_in_fcs_frame.npy")   # similarity-mapped onto the FCS pattern


def afterpulse_probability():
    """Raw-minus-FLCS autocorrelation excess, integrated 25 ns - 100 us, times the count rate."""
    from fcs_common import load_segment_G
    tau, Graw, _, z = load_segment_G(weight=0)
    _, Gfl, _, _ = load_segment_G(weight=1)
    Gr, Gf = Graw.mean(0), Gfl.mean(0)
    b = st["bg"] / st["rate"]
    dG = np.array([Gr[:, i, i] - (1 - b[i]) ** 2 * Gf[:, i, i] for i in range(23)]).T
    dtau = 2.0 ** z["casc"] / float(z["sync"])
    m = (tau >= 25e-9) & (tau <= 1e-4)
    return (dG[m] * dtau[m, None]).sum(0) * st["rate"]


p_ap = afterpulse_probability()

# ---------------- table
rows = []
for i in range(23):
    rows.append(dict(
        routing=int(ROUT[i]), nominal_x_pitch=round(xy[i, 0], 3), nominal_y_pitch=round(xy[i, 1], 3),
        mdf_x_nm=round(S[i, 0], 1), mdf_y_nm=round(S[i, 1], 1),
        mdf_x_rel_axis_nm=round(S[i, 0] - axis[0], 1), mdf_y_rel_axis_nm=round(S[i, 1] - axis[1], 1),
        sigma_xy_nm=round(sig[i], 1), sigma_xy_err_nm=round(sig_err[i], 1), sigma_z_nm=round(sz[i], 0),
        tcspc_offset_ps=round(offs_ps[i], 0), lifetime_ns=round(float(st["tau_tail"][i]), 3),
        rate_cps=round(float(st["rate"][i]), 0), flat_bg_cps=round(float(st["bg"][i]), 1),
        dark_cps=round(float(st["dark"][i]), 1),
        stray_bg_cps=round(float(st["bg"][i] - st["dark"][i]), 1),
        afterpulse_prob=round(float(p_ap[i]), 5),
        crosstalk_nn_mean=float("%.2e" % nn_eps[i]), crosstalk_total=float("%.2e" % nn_sum[i])))
with open(OUT / "pda23_calibration_table.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


def hexmap(ax, vals, title, cmap="viridis", fmt="%.0f", vmin=None, vmax=None, labels=True):
    patches = [RegularPolygon((x, y), 6, radius=0.56, orientation=0) for x, y in xy]
    pc = PatchCollection(patches, cmap=cmap, edgecolor="w", lw=1)
    pc.set_array(np.asarray(vals)); pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    if labels:
        for (x, y), v, r in zip(xy, vals, ROUT):
            ax.text(x, y + 0.12, fmt % v, ha="center", va="center", fontsize=7, color="w")
            ax.text(x, y - 0.2, str(r), ha="center", va="center", fontsize=5.5, color="w", alpha=0.8)
    ax.set_xlim(-2.7, 2.7); ax.set_ylim(-2.3, 2.3); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=9)
    plt.colorbar(pc, ax=ax, fraction=0.04)


# ---------------- figure 1: detector electronics
fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))
hexmap(ax[0], offs_ps, "TCSPC timing offset (ps, vs median)", "coolwarm", vmin=-90, vmax=90)
hexmap(ax[1], st["dark"], "dark counts (cps, dark file)", "magma", "%.0f", vmin=0, vmax=60)
hexmap(ax[2], p_ap * 100, "afterpulsing probability (%, 25 ns-100 us)", "viridis", "%.2f", vmin=0, vmax=0.2)
# crosstalk graph
for i in range(23):
    for j in range(i + 1, 23):
        e = eps[i, j]
        if e > 1e-3:
            ax[3].plot(*xy[[i, j]].T, "-", lw=1 + 6 * (e - 1e-3) / 3e-3, color=plt.cm.plasma((e - 2e-3) / 1.6e-3))
ax[3].scatter(*xy.T, s=180, c="k", zorder=3)
for (x, y), r in zip(xy, ROUT):
    ax[3].text(x, y, str(r), color="w", ha="center", va="center", fontsize=6, zorder=4)
ax[3].set_aspect("equal"); ax[3].axis("off")
ax[3].set_title("nearest-neighbour crosstalk %.2f-%.2f %% (lines)\nnext-nearest %.3f %%" % (
    np.nanmin(eps[(dist > 0.9) & (dist < 1.1)]) * 100, np.nanmax(eps[(dist > 0.9) & (dist < 1.1)]) * 100,
    np.nanmedian(eps[(dist > 1.6) & (dist < 2.1)]) * 100), fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "fig1_detector_electronics.png", dpi=140); plt.close(fig)

# ---------------- figure 2: shift vectors and MDF widths from FCS
sc = g["shift_per_pitch_um"] * 1e3
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
ideal = xy * sc
a.plot(*ideal.T, "o", mfc="none", mec="gray", ms=11, label=f"ideal lattice, {sc:.1f} nm/pitch")
a.plot(*S.T, "o", color="C0", ms=5, label="FCS MDF centres")
res = S - (ideal - ideal.mean(0) + S.mean(0))
a.quiver(*S.T, *(res * 8).T, angles="xy", scale_units="xy", scale=1, color="C3", width=0.004, label="residual x8")
a.plot(*axis, "k+", ms=14, mew=2, label="excitation axis (width minimum)")
for (x, y), r in zip(S, ROUT):
    a.text(x + 6, y + 6, str(r), fontsize=6)
a.set_aspect("equal"); a.legend(fontsize=7, loc="lower left"); a.set_xlabel("x (nm, chip frame)"); a.set_ylabel("y (nm)")
a.set_title(f"FCS-derived MDF positions (D = {g['D']:.0f} um^2/s, 24 C)\nrms residual {g['rms_similarity_um'] * 1e3:.1f} nm, anisotropy "
            f"{(g['affine_scales_um'][0] / g['affine_scales_um'][1] - 1) * 100:.1f} %", fontsize=9)
hexmap(ax[1], sig, "lateral MDF sigma (nm), Gaussian fit", "viridis", "%.0f")
hexmap(ax[2], sz / sig, "axial/lateral ratio K", "viridis", "%.0f")
fig.tight_layout(); fig.savefig(OUT / "fig2_fcs_shifts_widths.png", dpi=140); plt.close(fig)

# ---------------- figure 3: FCS vs beads
fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))
a = ax[0]
for k, b in enumerate(beads):
    a.plot(*b.T, ".", ms=4, color="C1", alpha=0.6, label="bead APR (5 images), similarity-mapped" if k == 0 else None)
a.plot(*S.T, "o", ms=5, color="C0", label="FCS")
bm = beads.mean(0)
a.quiver(*S.T, *((bm - S) * 3).T, angles="xy", scale_units="xy", scale=1, color="C3", width=0.004, label="bead mean - FCS (x3)")
a.set_aspect("equal"); a.legend(fontsize=7); a.set_title("pattern shape: beads vs FCS (bead scale matched to FCS)", fontsize=9)
a = ax[1]
cats = ["FCS", "bead0", "bead1", "bead2", "bead3", "bead5"]
bead_sc = [61.7, 60.5, 54.7, 55.3, 59.3]
bead_an = [27, 26, 11, 29, 13]
a.bar(np.arange(6) - 0.2, [sc] + bead_sc, 0.4, label="shift per pitch (nm, geometric mean)")
a2 = a.twinx(); a2.bar(np.arange(6) + 0.2, [1.3] + bead_an, 0.4, color="C3", label="anisotropy (%)")
a.set_xticks(range(6)); a.set_xticklabels(cats); a.set_ylabel("nm per pitch"); a2.set_ylabel("anisotropy %", color="C3")
a.set_title("magnitude and anisotropy", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "fig3_fcs_vs_beads.png", dpi=140); plt.close(fig)

summary = dict(
    D_um2s_24C=g["D"], D_um2s_25C=g25["D"],
    shift_per_pitch_nm_24C=sc, shift_per_pitch_nm_25C=g25["shift_per_pitch_um"] * 1e3,
    rms_nonsimilarity_nm=g["rms_similarity_um"] * 1e3,
    affine_scales_nm=[v * 1e3 for v in g["affine_scales_um"]],
    excitation_axis_nm=axis.tolist(), excitation_axis_pitch=(axis / sc).tolist(),
    concentration_nM=g["c_nM"], chi2r=g["chi2r"],
    crosstalk_nn_range=[float(np.nanmin(eps[(dist > 0.9) & (dist < 1.1)])), float(np.nanmax(eps[(dist > 0.9) & (dist < 1.1)]))],
    timing_offset_ps_range=[float(offs_ps.min()), float(offs_ps.max())],
    afterpulse_prob_median=float(np.median(p_ap)), stray_bg_cps_median=float(np.median(st["bg"] - st["dark"])),
)
json.dump(summary, open(OUT / "calibration_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))

# ---------------- figure 4: physical model
phys = json.load(open(OUT / "step4b_sat_T24.json"))
cons = json.load(open(OUT / "step5_physical_consistency.json"))
zw = np.load(OUT / "step4b_sat_T24.npz")
t4, Y4, E4, M4, iu4 = zw["t"], zw["Y"], zw["E"], zw["model"], zw["iu"]
fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
for j, col in [(11, "k"), (12, "C0"), (13, "C2"), (22, "C3")]:
    a_, b_ = min(11, j), max(11, j)
    k = np.flatnonzero((iu4[0] == a_) & (iu4[1] == b_))[0]
    ax[0].errorbar(t4, Y4[:, k], E4[:, k], fmt=".", ms=3, color=col, alpha=0.6, label=f"20 x {j + 9}")
    ax[0].plot(t4, M4[:, k], "-", color=col)
# mean over all pairs at the same nominal distance, to show the fit without single-pair noise
dd = np.round(dist[iu4[0], iu4[1]], 2)
for dv, col in [(0.0, "k"), (1.0, "C0"), (2.0, "C2"), (2.65, "C1"), (4.0, "C3")]:
    m_ = np.abs(dd - dv) < 0.05
    ax[1].plot(t4, Y4[:, m_].mean(1), ".", ms=3, color=col, label=f"d = {dv} pitch (n={m_.sum()})")
    ax[1].plot(t4, M4[:, m_].mean(1), "-", color=col)
for a in ax[:2]:
    a.set_xscale("log"); a.set_xlabel("tau (s)"); a.set_ylabel("G"); a.legend(fontsize=7)
pp = phys["params"]; pe = phys["errors"]
ax[0].set_title("wave-optics fit, single pairs", fontsize=9)
ax[1].set_title(f"pair-averaged by distance; chi2r = {phys['chi2r']:.3f}\n"
                f"d = {pp['d']:.1f}+/-{pe['d']:.1f} um, pitch = {pp['pitch'] * 1e3:.1f}+/-{pe['pitch'] * 1e3:.1f} nm, "
                f"sat = {pp['sat']:.2f}+/-{pe['sat']:.2f}, dz_det = {pp['dzd'] * 1e3:.0f}+/-{pe['dzd'] * 1e3:.0f} nm", fontsize=8)
lab = ["in focus,\nunsaturated\n(model)", "beads,\nmeasured", "FCS 18 um,\nno sat. (model)", "FCS 18 um,\nfitted model", "FCS,\nmeasured"]
val = [cons["bead_unsat"], 58.3, 53.5, 70.1, 75.9]
ax[2].bar(range(5), val, color=["C7", "C1", "C7", "C0", "C3"])
for i, v in enumerate(val):
    ax[2].text(i, v + 1, f"{v:.1f}", ha="center", fontsize=8)
ax[2].set_xticks(range(5)); ax[2].set_xticklabels(lab, fontsize=7); ax[2].set_ylabel("shift per pitch (nm)")
ax[2].set_title(f"one pitch ({pp['pitch'] * 1e3:.0f} nm) explains all: saturation x1.31", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "fig4_physical_model.png", dpi=140); plt.close(fig)
summary.update(physical_fit=dict(params=pp, errors=pe, chi2r=phys["chi2r"]),
               shift_ladder_nm=dict(zip(["in_focus_unsat_model", "beads_measured", "fcs_nosat_model", "fcs_fitted_model", "fcs_measured_1us_3ms"], val)))
json.dump(summary, open(OUT / "calibration_summary.json", "w"), indent=1)
