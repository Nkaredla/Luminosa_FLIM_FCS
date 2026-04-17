# Luminosa_FLIM_FCS

Luminosa_FLIM_FCS is a MATLAB toolbox for fluorescence lifetime imaging microscopy (FLIM) and fluorescence correlation spectroscopy (FCS) workflows built around TTTR and TCSPC data. It includes PTU and HT3 readers, IRF models and estimation, multi-exponential reconvolution fitting, global pattern matching, image scanning microscopy (APR + ACO-ISM), and line-scan FCS utilities.

## Quickstart

Read a PTU MultiFrame scan, reconstruct APR + ACO-ISM, and build reassigned FLIM:

```matlab
name = 'path/to/data.ptu';

ptuOut = PTU_MultiFrameScanReadFast(name, 1e6, true); % store tcspc_pix for IRF fitting

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

If you do not need `tcspc_pix`, set `storeTcspcPix = false` in `PTU_MultiFrameScanReadFast` to save memory.

## Bundled MIET GUI

The MIET GUI code is vendored in `external/miet-gui`, so a normal clone of this repository includes it. No `git submodule update` step is required.

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
- Image Processing Toolbox is required for `imgaussfilt`, `regionfill`, `montage`, `drawline`, and `improfile` (ISM and line-profile tools).
- Signal Processing Toolbox may be required for `hann` depending on MATLAB version.
- `AutodetectTimeGates.m` uses `mCluster` (not included). Provide an equivalent on your MATLAB path or replace the gate detection.

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
- `imgStack`, `channelIDs`.
- `tempReferenceIndex`, `centerDetectorIndex`.
- `shiftsToTempRef`, `shiftsToCenter`, `detectorPositions`.
- `rawSum`, `aprImage`, `acoAverage`, `acoAverageMasked`, `acoImage`.
- `deconvolvedImage`, `deconvFilter` (empty if `params.doISMDeconv = false`).
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
Notes:
- `lsFCS` now plots spatio-temporal correlations (tau vs xi) live in a `cnum x dnum` subplot grid.
- Use `lsCrossRead(res, flag)` to compute `G`, `Gcross`, and carpet correlations.

## Workflows

### ISM reconstruction and reassigned FLIM

1) Read PTU and build detector images.
2) Reconstruct APR + ACO-ISM.
3) Reassign photons and compute lifetime maps.

Reference for ACO-ISM: Ancora et al., *Image scanning microscopy reconstruction by autocorrelation inversion*, Journal of Physics: Photonics 6(4):045003 (2024), DOI: 10.1088/2515-7647/ad68dd.

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
res = lsFCS('data.ptu', 1, 10, []);
[G, Gcross, Gcarp, GcarpCross, t, xxi] = lsCrossRead(res, 25);
```

`lsFCS` displays a live `cnum x dnum` grid of spatio-temporal correlation maps (`tau` vs `xi`) while processing. The time axis uses decade ticks (e.g., 1 ms, 10 ms, 100 ms when in range).

### Batch line-scan FCS (example script)

`run_lsFCS.m` loops over subfolders in a root directory and runs `lsFCS` on the single `.ptu` file in each subfolder. It can optionally fit a single-diffusion model for components 1 and 4 and save PNG figures with `w0` (Gaussian sigma) and diffusion coefficient values in the title.

### Batch ISM on beads (example script)

`run_beads_batch.m` is a batch-processing example that loops through a folder of PTU files, reconstructs APR and ACO-ISM images, writes PNG and TIFF outputs, and measures bead widths.

Notes:
- Requires `CombineImages.m`, `save_tiff.m`, and `Gauss2D.m`.
- The script uses `results.deconvolvedImage`, which is produced when `params.doISMDeconv` is true (default). If you disable deconvolution, swap in `results.acoImage` or `results.aprImage`.

### Batch FLIM quicklook (example script)

`FLIM_read.m` is a simple batch script for reading PTU files in a folder, building intensity and lifetime maps, adding a scale bar, and saving a PNG plus a MAT file per dataset. Edit `folderName` at the top of the script before running.

### Line profile across images

`lineProfileAcrossImages.m` lets you draw one line on a reference image and extract intensity profiles along the same line in all images.

```matlab
titles = {'Raw sum','APR','ACO-ISM'};
ims = cat(3, results.rawSum, results.aprImage, results.acoImage);
[profiles, distPix, pos] = lineProfileAcrossImages(ims, 2, 600, titles);
```

## Function Reference

### Readers and I/O

`PTU_Read_Head.m`
Signature: `head = PTU_Read_Head(name, verbose)`
Summary: Reads PicoQuant PTU headers into a struct with all tag fields and `head.length`.
Notes: `verbose` defaults to `false` (wide-string tag output is suppressed unless enabled).

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

`LPTU_LineScanReadOld.m`
Signature: `[head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame, num] = LPTU_LineScanReadOld(name, cnts, nx)`
Summary: Legacy PTU line-scan reader kept for backwards compatibility (defaults to `nx = 50`).
Notes: Prefer `LPTU_LineScanRead` for current workflows.

`Harp_tcspc.m`
Signature: `[bin, tcspcdata, head] = Harp_tcspc(name, resolution, deadtime, photons)`
Summary: Builds TCSPC histograms from `.ht3` or `.ptu` files with optional deadtime filtering.
Notes: Caches results to `*.ht3tcspc`. Requires `HT3_Read.m` and `mHist.m`.

### ISM and FLIM reconstruction

`run_ism_reconstruction_from_ptu.m`
Signature: `results = run_ism_reconstruction_from_ptu(ptuOut, params)`
Summary: APR + ACO-ISM reconstruction with phase-correlation registration and Schultz-Snyder inversion.
Key params: `imageSource`, `frameIndices`, `frameCombine`, `dropEmptyChannels`, `centerDetectorIndex`, `tempReferenceIndex`, `smoothSigma`, `useWindow`, `normalizeImages`, `upsampleReg`, `nIter`, `checkEvery`, `stopTol`, `minIter`, `supportMask`, `preserveFlux`, `pixelSize`, `lambda`, `NA`, `doISMDeconv`, `deconvLambda`, `deconvClipNegative`, `deconvPreserveFlux`, `showPlots`.
Outputs: `results.aprImage`, `results.acoImage`, `results.deconvolvedImage`, `results.shiftsToCenter`, plus OTF and registration diagnostics.

`reassigned_flim.m`
Signature: `flim = reassigned_flim(ptuOut, ismRes, params)`
Summary: Reassigns photons using APR shifts and computes lifetime statistics on native or oversampled grids.
Key params: `oversampleXY`, `keepSameSize`, `storeTotalCubes`, `storeFrameCubes`, `frameIndices`, `minCounts`, `useBackground`, `bgBins`, `t0Mode`, `t0Bin`, `overflowAction`.

`run_ISM.m`
Type: Script
Summary: Example end-to-end pipeline for PTU read, ISM reconstruction, FLIM reassignment, and IRF estimation.

`FLIM_read.m`
Type: Script
Summary: Batch PTU quicklook that builds intensity/lifetime maps, adds a scale bar, and saves PNG and MAT outputs. Edit `folderName` in the script.

`run_beads_batch.m`
Type: Script
Summary: Batch ISM reconstruction over a folder of PTU files, saves images, and estimates bead widths.
Notes: Uses `CombineImages.m`, `save_tiff.m`, and `Gauss2D.m`, and references `results.deconvolvedImage`.

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

### Histogram utilities

`mHist.m`
Signature: `[z, xv] = mHist(x, xv, weight)`
Summary: 1D histogram with optional weights and custom bin centers/edges.

`mHist2.m`
Signature: `[z, xv, yv] = mHist2(x, y, xv, yv, weight)`
Summary: 2D histogram with optional weights and custom bin vectors.

`mHist3.m`
Signature: `[h, xv, yv, zv] = mHist3(x, y, z, xv, yv, zv)`
Summary: 3D histogram with custom bin vectors or sizes.

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
Signature: `res = lsFCS(name, cnum, maxtime, timegates)`
Summary: Line-scan FCS with multi-tau correlation; live spatio-temporal plotting in a `cnum x dnum` grid.
Outputs: `res.tcspc`, `res.autotime`, `res.auto`, `res.automean`, `res.rate`, `res.time`, `res.head`.

`AutodetectTimeGates.m`
Signature: `[t1, len] = AutodetectTimeGates(tcspcdata, cnum)`
Summary: Heuristic gate detection for multi-pulse TCSPC data (returns gate start indices and a common length).
Notes: Uses `mCluster` (not included).

`lsCrossRead.m`
Signature: `[G, Gcross, Gcarp, GcarpCross, t, xxi] = lsCrossRead(res, flag)`
Summary: Reorders line-scan correlation output into correlation vs spatial shift and component pairs.

`tttr2xfcs.m`
Signature: `[auto, autotime] = tttr2xfcs(y, num, Ncasc, Nsub)`
Summary: Multi-tau auto- and cross-correlation for weighted photon streams.

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
Notes: Calls `mim` (included).

`mim.m`
Signature: `mim(x, p1, p2, p3)`
Summary: Minimal image display helper used by `cim` and `CombineImages`.

`addPTUScaleBar.m`
Signature: `h = addPTUScaleBar(axOrIm, head, imSize, location, varargin)`
Summary: Draws a scale bar on an image using PTU header pixel size metadata.

`lineProfileAcrossImages.m`
Signature: `[profiles, distPix, pos] = lineProfileAcrossImages(imStack, refIdx, nSamples, imgTitles)`
Summary: Interactive line selection on a reference image and profile extraction across a stack.
Notes: Requires Image Processing Toolbox for `drawline` and `improfile`.

`CombineImages.m`
Signature: `im = CombineImages(imraw, n, m, flag, labelx, labely, fsize)`
Summary: Tiles an image stack into an `n`-by-`m` mosaic with optional scaling and labels.

`CombineImagesMultiCmap.m`
Signature: `imRGB = CombineImagesMultiCmap(imraw, n, m, cmaps, flag, labelx, labely, fsize, clims)`
Summary: Tiles an image stack into a mosaic and applies a separate colormap per panel.

`ShowImagesMultiCmapWithColorbars.m`
Signature: `ShowImagesMultiCmapWithColorbars(imraw, cmaps, clims, titlestr, blackbg)`
Summary: Displays a stack with per-panel colormaps and colorbars.

`save_tiff.m`
Signature: `save_tiff(filename, imgdata, datatype, bitdepth, force_overwrite)`
Summary: Saves an image or stack to TIFF with datatype and bitdepth control.
