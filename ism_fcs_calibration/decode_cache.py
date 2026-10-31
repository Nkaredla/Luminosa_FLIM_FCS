"""Decode the Atto 655 PTU and the dark PTU once into .npy caches."""
import sys, time, numpy as np
from pathlib import Path
from ptu_io import read_t3

CACHE = Path(sys.argv[1])
FILES = {
    "atto655": r"D:\Luminosa\Data\270226\ISM\Atto655_ISM_26\ISM_IRF_20260227-201810\RawData.ptu",
    "dark": r"D:\Luminosa\Data\ISMdark_counts.ptu",
}
for key, path in FILES.items():
    t0 = time.time()
    d = read_t3(path)
    tags = d["tags"]
    np.save(CACHE / f"{key}_nsync.npy", d["nsync"])
    np.save(CACHE / f"{key}_chan.npy", d["chan"])
    np.save(CACHE / f"{key}_dtime.npy", d["dtime"])
    np.save(CACHE / f"{key}_markers.npy", np.stack([d["marker_nsync"], d["marker_chan"].astype(np.int64)]))
    sync = tags["TTResult_SyncRate"]
    T = d["nsync"][-1] / sync
    ch, cnt = np.unique(d["chan"], return_counts=True)
    print(f"== {key}: {len(d['nsync'])} photons, {d['marker_nsync'].size} markers, T={T:.1f} s, sync={sync}, res={tags['MeasDesc_Resolution']}, decode {time.time()-t0:.1f}s")
    for c, n in zip(ch, cnt):
        print(f"   routing {c:2d}: {n:10d}  {n/T:9.1f} cps   dtime max {d['dtime'][d['chan']==c].max()}")
