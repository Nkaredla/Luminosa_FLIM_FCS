# B12 tail-only fitting

Run `run_B12_MitoLysoCytosol_TailFLIM.m` in this folder after setting its PTU
path and segmentation parameters. It reuses the main program's channel
extraction and segmentation, then stops before IRF estimation or Bayesian
pixel processing. Tail-only results and figures are saved to `B12_TailFLIM`
beside the PTU file by default.

Alternatively, reuse an existing main-program result without reading the PTU
or repeating segmentation:

```matlab
addpath('b12');
tailResult = PTU_MitoLysoCytosol_TailFLIM(result);
```

Each segmented mitochondrial region and lysosome has its own TCSPC peak and
tail cutoff. Cytosol uses its pooled mask, as in the main program. Global and
pooled organelle fits are diagnostics and do not fill skipped regions.

- `cfg.tailOffsetNs = 1`: retain bins starting at least 1 ns after the peak.
  The offset is rounded up to whole bins (1.024 ns for 0.128 ns bins).
- `cfg.minPhotonsPerRegionFit = 10`: fewer photons means no exponential fit.
- `cfg.minPhotonsForTwoExp = 500`: 10–499 photons permits only a single
  exponential; 500 or more permits the existing one/two-exponential BIC
  selection. These settings can be raised, but not below the 10/500 floors.

Tail-only eligibility counts photons **after** the cutoff; the main program
counts the full gated curve. Too few retained time bins also skips the
corresponding model. Unfitted regions have model order zero and NaN lifetime.
Excluded bins and skipped curves have NaN fitted counts.

Tail fits use unconvolved exponentials with optional constant background.
Component fractions and amplitudes refer to the retained tail; they are not
extrapolated to the original peak. `result.fits` includes each curve's
`peakBin`, `firstFitBin`, `firstFitTimeNs`, `totalPhotonCount`, and
`retainedPhotonCount`. A first-fit bin beyond the histogram means no tail is
available and its first-fit time is NaN.

The main runner `run_B12_MitoLysoCytosol_BayesFLIM.m` also exposes the 10/500
cutoffs. Its Bayesian maps use only eligible regional seeds, so pooled or
global models cannot bypass the cutoffs.

Run `test_B12_tail_fitting` to check photon boundaries, independent cutoff
times, synthetic lifetime recovery, empty tails, and skipped-region maps.
