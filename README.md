# Luminosa_FLIM_FCS

Luminosa_FLIM_FCS is a MATLAB toolbox for fluorescence lifetime imaging microscopy (FLIM) and fluorescence correlation spectroscopy (FCS) workflows built around TTTR and TCSPC data. It includes PTU and HT3 readers, IRF models and estimation, multi-exponential reconvolution fitting, global pattern matching, image scanning microscopy (APR + ACO-ISM), and line-scan FCS utilities.

## Quickstart

Read a PTU MultiFrame scan, reconstruct APR + ACO-ISM, and build reassigned FLIM:

```matlab
name = 'path/to/data.ptu';

ptuOut = PTU_MultiFrameScanReadFast(name, 1e6, false);

params = struct();
params.imageSource = 'tags';
params.smoothSigma = 1;
params.useWindow = true;
params.normalizeImages = true;
params.upsampleReg = 20;

params.nIter = 1000;
params.checkEvery = 25;
params.stopTol = 1e-7;
params.minIter = 50;

params.pixelSize = ptuOut.head.ImgHdr_PixResol * 1e3; % nm
params.lambda = 690; % nm
params.NA = 1.45;

params.showPlots = true;

ismRes = run_ism_reconstruction_from_ptu(ptuOut, params);

flim = reassigned_flim(ptuOut, ismRes, params);
```

Estimate a global IRF and lifetimes from a TCSPC cube:

```matlab
tcspc_pix = ptuOut.tcspc_pix;
outIRF = Calc_mIRF_Global_GammaShifted(ptuOut.head, tcspc_pix, [0.35 1.5 5]);
irf = outIRF.IRF;
```

Global multi-exponential pattern matching:

```matlab
pulsePeriodNs = ptuOut.head.MeasDesc_GlobalResolution * 1e9;
dtNs = ptuOut.head.MeasDesc_Resolution * 1e9;
tau0 = [0.4 2.0];

out = GlobalMultiExpPatternMatchFromTCSPC(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0);
```

## Requirements

MATLAB is required.

Optional toolboxes and features:
- Parallel Computing Toolbox is required for `parpool`, `parfeval`, and GPU paths (`gpuArray`, `gpuDeviceCount`, `pagemtimes`).
- Image Processing Toolbox is required for `imgaussfilt`, `regionfill`, and `montage` used in ISM reconstruction.
- Signal Processing Toolbox may be required for `hann` depending on MATLAB version.

External files not included in this repository (required by some functions):
- `HT3_Read.m`
- `mHist.m`
- `mHist3.m`
- `PTU_LineScanRead.m`
- `PTU_LineScanCorr.m`
- `tttr2xfcs.m`
- `AutodetectTimeGates.m`
- `mim.m`
- `mexUtil.h` (required to build `cIntersect.cpp`)

If you only use the PTU MultiFrame readers and the ISM and FLIM pipeline, the missing line-scan and HT3 helpers are not required.

## Units and Conventions

- Time inputs and outputs are usually in nanoseconds (ns) unless stated otherwise.
- PTU headers store `MeasDesc_Resolution` in seconds; readers convert to ns internally where needed.
- TCSPC cubes are typically `[nx, ny, nTime, nCh]` and stored as `uint16`.
- Photon-level arrays use 1-based indices consistent with MATLAB.

## Core Data Structures

### ptuOut

`ptuOut` is returned by `PTU_MultiFrameScanRead` or `PTU_MultiFrameScanReadFast`.

Fields (most common):
- `head` PTU header struct from `PTU_Read_Head`.
- `tag` and `tau` per-frame maps `[nx, ny, nCh, nFrames]`.
- `tags` and `taus` global (frame-summed) maps `[nx, ny, nCh]`.
- `time` cell array of per-channel photon times.
- `tcspc_pix` optional TCSPC cube `[nx, ny, nTime, nCh]`.
- `im_sync`, `im_tcspc`, `im_chan`, `im_line`, `im_col`, `im_frame` photon-level arrays.
- `dind` active channel indices.

### ismRes

`ismRes` is returned by `run_ism_reconstruction_from_ptu`.

Key fields:
- `imgStack`, `channelIDs`, `shiftsToCenter`.
- `rawSum`, `aprImage`, `acoAverage`, `acoImage`.
- `convHistory`, `otfWF`, `otfISMexp`, `otfISMideal`.
- `paramsUsed`.

### flim

`flim` is returned by `reassigned_flim`.

Key fields:
- `coordMode`, `tAxisNs`.
- `unassigned.total` and `reassigned.total` stats: `tag`, `meanArrival`, `tauMean`, `tauRMS`, `globalDecay`, `tAxisNs`, `t0Bin`, and optional `xyT`.
- `unassigned.frames` and `reassigned.frames` per-frame stats and optional cubes.
- `total` and `frames` are aliases of `reassigned.total` and `reassigned.frames`.

### res (line-scan FCS)

`res` is returned by `lsFCS`.

Key fields:
- `tau`, `tcspc`, `autotime`, `auto`, `automean`, `rate`, `time`, `head`, `line_idx`.

## Workflows

### ISM reconstruction and reassigned FLIM

1) Read PTU and build detector images.
2) Reconstruct APR + ACO-ISM.
3) Reassign photons and compute lifetime maps.

```matlab
ptuOut = PTU_MultiFrameScanReadFast(name, 1e6, false);
ismRes = run_ism_reconstruction_from_ptu(ptuOut, params);
flim = reassigned_flim(ptuOut, ismRes, params);
```

### Global IRF estimation

Choose one of the global IRF models:

```matlab
outIRF = Calc_mIRF_Global_ExGauss(head, tcspc_pix, [0.35 1.8 4.2]);
irf = outIRF.IRF;
```

```matlab
outIRF = Calc_mIRF_Global_GammaShifted(head, tcspc_pix, [0.35 1.5 5]);
irf = outIRF.IRF;
```

### Global multi-exponential pattern matching

```matlab
out = GlobalMultiExpPatternMatchFromTCSPC(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0);
tauMean = out.tauMeanArithmetic;
```

### Line-scan FCS

```matlab
res = lsFCS('data.ptu', 1, 10, [], 0);
[G, Gcross, Gcarp, GcarpCross, t, xxi] = lsCrossRead(res, 25);
```

## Function Reference

### Readers and I/O

`PTU_Read_Head.m`
Signature: `head = PTU_Read_Head(name)`
Summary: Reads PicoQuant PTU headers into a struct with all tag fields and `head.length`.

`PTU_Read.m`
Signature: `[sync, tcspc, chan, special, num, loc, head] = PTU_Read(name, cnts, head)`
Summary: Reads TTTR records and decodes sync index, TCSPC bin, channel, and marker events.
Inputs: `cnts` can be a scalar count or `[start count]`.
Outputs: `special` is nonzero for markers, `loc` is overcount information used for chunk stitching.

`PTU_MultiFrameScanRead.m`
Signature: `out = PTU_MultiFrameScanRead(name, photonsPerChunk)`
Summary: MultiFrame PTU reader that builds full `tcspc_pix` cubes and per-frame maps.
Notes: Assumes `head.ImgHdr_Ident == 9` and monodirectional scanning.

`PTU_MultiFrameScanReadFast.m`
Signature: `out = PTU_MultiFrameScanReadFast(name, photonsPerChunk, storeTcspcPix)`
Summary: Faster, lower-memory MultiFrame reader. Uses single-precision tag and tau maps and optionally omits `tcspc_pix`.
Notes: Assumes `head.ImgHdr_Ident == 9` and monodirectional scanning.

`LPTU_LineScanRead.m`
Signature: `[head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame, num] = LPTU_LineScanRead(name, cnts, nx)`
Summary: Reads photons from PTU line-scan files and assigns each photon to line and column indices.
Notes: Only monodirectional line scans (`ImgHdr_Ident == 9` and `ImgHdr_BiDirect == 0`).

`Harp_tcspc.m`
Signature: `[bin, tcspcdata, head] = Harp_tcspc(name, resolution, deadtime, photons)`
Summary: Builds TCSPC histograms from `.ht3` or `.ptu` files with optional deadtime filtering.
Notes: Caches results to `*.ht3tcspc`. Requires `HT3_Read.m` and `mHist.m`.

### ISM and FLIM reconstruction

`run_ism_reconstruction_from_ptu.m`
Signature: `results = run_ism_reconstruction_from_ptu(ptuOut, params)`
Summary: APR + ACO-ISM reconstruction with phase-correlation registration and Schultz-Snyder inversion.
Key params: `imageSource`, `smoothSigma`, `upsampleReg`, `nIter`, `stopTol`, `pixelSize`, `lambda`, `NA`, `showPlots`.
Outputs: `results.aprImage`, `results.acoImage`, `results.shiftsToCenter`, and diagnostics.

`reassigned_flim.m`
Signature: `flim = reassigned_flim(ptuOut, ismRes, params)`
Summary: Reassigns photons using APR shifts and computes lifetime statistics on native or oversampled grids.
Key params: `oversampleXY`, `keepSameSize`, `storeTotalCubes`, `storeFrameCubes`, `minCounts`, `useBackground`, `bgBins`, `t0Mode`, `t0Bin`, `overflowAction`.

`run_ISM.m`
Type: Script
Summary: Example end-to-end pipeline for PTU read, ISM reconstruction, FLIM reassignment, and IRF estimation.

### IRF models and estimation

`IRF_Fun.m`
Signature: `z = IRF_Fun(p, t, pic)`
Summary: Walther-style IRF model with Gaussian peak and exponential tail. Returns normalized IRF.

`IRF_ExGauss.m`
Signature: `z = IRF_ExGauss(p, t, pic)`
Summary: Ex-Gaussian IRF model (Gaussian convolved with a one-sided exponential).

`IRF_GammaShifted.m`
Signature: `z = IRF_GammaShifted(p, t, pic)`
Summary: Shifted gamma IRF model with shape and scale parameters.

`Calc_mIRF.m`
Signature: `IRF = Calc_mIRF(head, tcspc)`
Summary: Estimates IRF per channel (and per PIE) using multi-start Simplex on `TCSPC_Fun`.

`Calc_mIRF_Global_ExGauss.m`
Signature: `out = Calc_mIRF_Global_ExGauss(head, tcspc, tau0)`
Summary: Fits a single global ex-Gaussian IRF and lifetimes from the summed decay.

`Calc_mIRF_Global_GammaShifted.m`
Signature: `out = Calc_mIRF_Global_GammaShifted(head, tcspc, tau0)`
Summary: Fits a single global shifted-gamma IRF and lifetimes from the summed decay.

### TCSPC fitting and utilities

`Fluofit.m`
Signature: `[c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Fluofit(irf, y, p, dt, tau, lim, init)`
Summary: Multi-exponential reconvolution fit with optional initialization from `DistFluofit`.
Notes: Uses `Simplex`, `mlfit` or `lsfit`, `Convol`, and `DistFluofit`.

`DistFluofit.m`
Signature: `[cx, tau, offset, csh, z, t, err] = DistFluofit(irf, y, p, dt, shift, flag, bild, N)`
Summary: Distributed lifetime fit using nonnegative least squares and IRF shift search.

`TCSPC_Fun.m`
Signature: `[err, c, zz, z] = TCSPC_Fun(p, t, y, para)`
Summary: Reconvolution model used by `Calc_mIRF` with Walther IRF.

`TCSPC_Fun_ExGauss.m`
Signature: `[err, c, zz, z, IRF] = TCSPC_Fun_ExGauss(p, t, y, para)`
Summary: Reconvolution model using `IRF_ExGauss`.

`TCSPC_Fun_GammaShifted.m`
Signature: `[err, c, zz, z, IRF] = TCSPC_Fun_GammaShifted(p, t, y, para)`
Summary: Reconvolution model using `IRF_GammaShifted` and PIRLS-based nonnegative fitting.

`lsfit.m`
Signature: `[err, A, z] = lsfit(param, irf, y, p)`
Summary: Least-squares error model for reconvolution fits.

`mlfit.m`
Signature: `err = mlfit(param, irf, y, p)`
Summary: Maximum-likelihood error model for reconvolution fits.

`Simplex.m`
Signature: `[x, dx, steps] = Simplex(fname, x, xmin, xmax, tol, steps, varargin)`
Summary: Bounded Nelder-Mead optimizer used for IRF and decay fitting.

`Convol.m`
Signature: `y = Convol(irf, x)`
Summary: Periodic convolution used for reconvolution fitting.

### Pattern matching and PIRLS

`PatternMatchIm_matlab.m`
Signature: `[C, Z, info] = PatternMatchIm_matlab(y, M, mode, useGPU, batchSize)`
Summary: Pixelwise pattern matching for FLIM cubes with least squares, nonnegative, or PIRLS modes.

`PIRLSnonneg.m`
Signature: `[beta, k] = PIRLSnonneg(x, y, max_num_iter)`
Summary: Poisson IRLS solver with nonnegative constraints.

`PIRLSnonneg_batch_gpu_matlab.m`
Signature: `[Beta, info] = PIRLSnonneg_batch_gpu_matlab(M, Y, maxPirlsIter, maxNnlsIter, batchSize, useGPU)`
Summary: GPU-friendly batch PIRLS with projected-gradient NNLS and CPU fallback.

`GlobalMultiExpPatternMatchFromTCSPC.m`
Signature: `out = GlobalMultiExpPatternMatchFromTCSPC(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0, opts)`
Summary: Global multi-exponential fit followed by pixelwise pattern matching.
Outputs: `out.taufit`, `out.Amp`, `out.AmpFrac`, `out.tauMeanArithmetic`, `out.tauMeanHarmonic`, `out.recon`.

### FCS and correlation

`lsFCS.m`
Signature: `res = lsFCS(name, cnum, maxtime, timegates, flagparallel)`
Summary: Line-scan FCS with multi-tau correlation and optional parallel processing.
Outputs: `res.tcspc`, `res.autotime`, `res.auto`, `res.automean`, `res.rate`, `res.time`, `res.head`.

`lsCrossRead.m`
Signature: `[G, Gcross, Gcarp, GcarpCross, t, xxi] = lsCrossRead(res, flag)`
Summary: Reorders line-scan correlation output into correlation vs spatial shift and component pairs.

`tttr2xfcsSym.m`
Signature: `[auto, autotime] = tttr2xfcsSym(y, num, Ncasc, Nsub)`
Summary: Symmetrically normalized auto and cross correlation with multiple-tau binning.

`cIntersect.m`
Signature: `[values, idx1, idx2] = cIntersect(y, lag)`
Summary: MEX-accelerated intersection of `y` and `y + lag`.
Notes: Precompiled `cIntersect.mexw64` and `cIntersect.mexa64` are included. Rebuilding requires `mexUtil.h`.

### Display utilities

`cim.m`
Signature: `handle = cim(x, p1, p2, p3, clmp)`
Summary: Convenience image display wrapper around `imagesc` with optional overlays.
Notes: Calls `mim`, which is not included in this repository.