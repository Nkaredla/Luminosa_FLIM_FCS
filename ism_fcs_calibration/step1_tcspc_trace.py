"""Step 1: per-channel TCSPC, background vs dark file, timing offsets, intensity trace."""
import sys, json, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CACHE = Path(sys.argv[1]); OUT = Path("out")
SYNC = 40e6; RES = 20e-12; T = 1200.0
ch = np.load(CACHE / "atto655_chan.npy"); dt = np.load(CACHE / "atto655_dtime.npy")
ns = np.load(CACHE / "atto655_nsync.npy")
dch = np.load(CACHE / "dark_chan.npy"); dns = np.load(CACHE / "dark_nsync.npy")
Tdark = dns[-1] / 20e6
ROUT = np.arange(9, 32)
nb = 1250  # 25 ns / 20 ps
H = np.zeros((23, nb))
for k, r in enumerate(ROUT):
    H[k] = np.bincount(dt[ch == r], minlength=nb)[:nb]
t_ns = (np.arange(nb) + 0.5) * RES * 1e9
res = {}
peak = np.zeros(23); bg = np.zeros(23); dark = np.zeros(23); rate = np.zeros(23); tau_tail = np.zeros(23)
for k, r in enumerate(ROUT):
    h = H[k]
    sm = np.convolve(h, np.ones(5) / 5, "same")
    ip = int(np.argmax(sm))
    # background from the window before the rise (well away from the tail of the previous pulse)
    pre = slice(max(ip - 400, 0), ip - 100) if ip > 500 else slice(ip + 900, ip + 1150)
    bg[k] = np.median(h[pre]) * nb / T  # cps of a flat background
    rate[k] = h.sum() / T
    dark[k] = np.sum(dch == r) / Tdark
    # half-max rise time for timing offset
    half = 0.5 * (sm[ip] - np.median(h[pre])) + np.median(h[pre])
    j = ip
    while j > 0 and sm[j] > half: j -= 1
    peak[k] = (j + (half - sm[j]) / (sm[j + 1] - sm[j])) * RES * 1e9
    # tail lifetime by log-linear fit 1-8 ns after peak
    sel = slice(ip + 50, ip + 400)
    y = h[sel] - np.median(h[pre]); x = t_ns[sel]
    ok = y > 0
    tau_tail[k] = -1 / np.polyfit(x[ok], np.log(y[ok]), 1, w=np.sqrt(y[ok]))[0]
    print(f"rout {r:2d}: rate {rate[k]:7.1f}  flat bg {bg[k]:6.1f} cps  dark-file {dark[k]:6.1f} cps  rise50 {peak[k]:6.3f} ns  tail tau {tau_tail[k]:.3f} ns")
print("timing offsets rel. to median (ps):", np.round((peak - np.median(peak)) * 1e3).astype(int))
np.savez(OUT / "step1_tcspc.npz", H=H, t_ns=t_ns, peak=peak, bg=bg, dark=dark, rate=rate, tau_tail=tau_tail, rout=ROUT)

# intensity trace in 1 s and 10 ms bins (total) for stability / aggregate check
tsec = ns / SYNC
tr1 = np.bincount((tsec).astype(int), minlength=1200)[:1200]
tr10 = np.bincount((tsec * 100).astype(np.int64), minlength=120000)[:120000]
print("1 s trace: mean %.0f cps, std %.0f, min %.0f max %.0f; first 60 s mean %.0f, last 60 s mean %.0f" % (
    tr1.mean(), tr1.std(), tr1.min(), tr1.max(), tr1[:60].mean(), tr1[-60:].mean()))
mu = tr10.mean(); print("10 ms bins: mean %.1f, poisson sd %.1f, measured sd %.1f, max %d, n>mu+8sd: %d" % (
    mu, np.sqrt(mu), tr10.std(), tr10.max(), np.sum(tr10 > mu + 8 * np.sqrt(mu))))
np.savez(OUT / "step1_trace.npz", tr1=tr1, tr10=tr10)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
for k in range(23):
    ax[0].semilogy(t_ns, H[k] + 1, lw=0.6)
ax[0].set_xlabel("micro time (ns)"); ax[0].set_ylabel("counts"); ax[0].set_title("TCSPC, 23 channels")
ax[1].bar(np.arange(23), bg, label="flat TCSPC background")
ax[1].plot(np.arange(23), dark, "ko", label="dark-count file")
ax[1].set_yscale("log"); ax[1].set_xticks(range(23)); ax[1].set_xticklabels(ROUT, fontsize=7)
ax[1].set_xlabel("routing channel"); ax[1].set_ylabel("cps"); ax[1].legend(); ax[1].set_title("uncorrelated background")
ax[2].plot(np.arange(1200), tr1 / 1e3, lw=0.7); ax[2].set_xlabel("time (s)"); ax[2].set_ylabel("total kcps"); ax[2].set_title("intensity trace (1 s bins)")
fig.tight_layout(); fig.savefig(OUT / "step1_tcspc_trace.png", dpi=130)
