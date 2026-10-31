# PDA-23 calibration from Atto 655 solution FCS (Luminosa, 2026-02-27)

Data: `D:\Luminosa\Data\270226\ISM\Atto655_ISM_26\ISM_IRF_20260227-201810\RawData.ptu`
(20 min point FCS, 640 nm / 40 MHz, 100x oil NA 1.45, focus ~18 um into water, laser 100 % / attenuator 100 %).
Dark counts: `D:\Luminosa\Data\ISMdark_counts.ptu`. Bead cross-check: `GattaBeads_25\good\RawImage*_detector_map.csv`.

Run order (Python 3.10 + numpy/scipy; the cache directory is any scratch folder):

| step | script | output |
|---|---|---|
| decode | `decode_cache.py <cache>` | `.npy` photon caches |
| 1 | `step1_tcspc_trace.py <cache>` | TCSPC, timing offsets, backgrounds, trace |
| 2 | `step2_correlate.py <cache>` | FLCS-filtered 23x23 multi-tau correlations per 20 s segment (~25 min) |
| 2b | `step2b_crosstalk.py <cache>` | ns coincidences -> crosstalk matrix |
| 3 | `step3_fit_gauss.py 24` | 3-D Gaussian-MDF global fit: MDF positions, widths |
| 3b/3c | `step3b_quarters.py`, `step3c_anisotropy.py` | stability, lateral ellipticity test |
| 4b | `step4b_fit_wave_parallel.py 24 3e-3 sat` | wave-optics fit (depth SA, pitch, offset, saturation; ~45 min) |
| 5 | `step5_physical_consistency.py` | model predictions with/without saturation, bead APR |
| report | `make_report.py` | `out/pda23_calibration_table.csv`, `out/calibration_summary.json`, `out/fig1-4*.png` |

`D` = 4.26e-6 cm^2/s (Atto 655-COOH, 25 C, 2fFCS) scaled by T/eta(T); 24 C -> 415 um^2/s.
All FCS lengths scale with sqrt(D).

Key caveat: the FCS was taken at saturating power (fitted 3.4 photons absorbed per pulse at
the focus) 18 um deep with an oil objective. The measured MDF shifts (77 nm per pitch) are valid
for that condition only; the physical model gives a sample-plane pitch of 109 nm, i.e. ~51 nm
per pitch APR shift for unsaturated in-focus imaging.
