function out = PTU_MultiFrameScanReadLog(name, photonsPerChunk, storeTcspcPixMT, useGPU, keepNativePhotonTcspc, binsPerOctave, keepPhotonLists)
% PTU_MultiFrameScanReadFast_multiTauTCSPC
%
% Read a PTU multi-frame scan, preserve native photon TCSPC bins, and store
% the image TCSPC as a compressed multi-tau cube instead of a dense linear
% tcspc_pix cube. Full resolution is preserved up to the global peak, then
% multi-tau binning starts after the peak.
%
% Why this function exists
% ------------------------
% A dense tcspc_pix of size [nx, ny, Ngate, nCh] can become enormous.
% This function avoids that by:
%   1) reading the PTU with native photon TCSPC preserved,
%   2) computing tag/tau maps as usual,
%   3) building only a compressed multi-tau TCSPC image cube for display/QC.
%
% Important caveat
% ----------------
% The stored tcspc_pix_mt is NOT a drop-in replacement for algorithms that
% assume a uniformly sampled linear TCSPC axis (for example some exact IRF-
% convolved fitting routines). It is intended for display, ROI extraction,
% and compact storage. For exact reconvolution fits, use im_tcspc_native.
%
% Inputs
% ------
% name                 : PTU filename
% photonsPerChunk      : chunk size for PTU_Read
% storeTcspcPixMT      : store multitau cube (default true)
% useGPU               : forwarded to the base reader (default false)
% keepNativePhotonTcspc: preserve native per-photon TCSPC bins (default true)
% binsPerOctave        : number of bins before doubling width (default 8)
% keepPhotonLists      : keep photon lists in output (default true)
%
% Output additions relative to PTU_MultiFrameScanReadFast_nativeTCSPC
% -------------------------------------------------------------------
% out.tcspc_pix_mt             : [nx, ny, nMt, nCh] uint16 multitau cube
% out.tcspc_mt_centers_ns      : multitau time-bin centers [ns]
% out.tcspc_mt_edges_ns        : multitau bin edges [ns]
% out.tcspc_mt_width_ns        : multitau bin widths [ns]
% out.tcspc_mt_edges_nativebin : multitau bin edges in native-bin units
% out.tcspc_pix                : []  (dense linear cube intentionally omitted)
%
% See also
% --------
% PTU_MultiFrameScanReadFast_nativeTCSPC

if nargin < 2 || isempty(photonsPerChunk)
    photonsPerChunk = 1e6;
end
if nargin < 3 || isempty(storeTcspcPixMT)
    storeTcspcPixMT = true;
end
if nargin < 4 || isempty(useGPU)
    useGPU = false;
end
if nargin < 5 || isempty(keepNativePhotonTcspc)
    keepNativePhotonTcspc = true;
end
if nargin < 6 || isempty(binsPerOctave)
    binsPerOctave = 8;
end
if nargin < 7 || isempty(keepPhotonLists)
    keepPhotonLists = true;
end

% Read with native photon TCSPC preserved, but do NOT build dense tcspc_pix.
out = PTU_MultiFrameScanReadFast_nativeTCSPC( ...
    name, photonsPerChunk, false, useGPU, keepNativePhotonTcspc);

% Remove dense cube if the base reader happened to create one.
out.tcspc_pix = [];

if ~storeTcspcPixMT
    out.tcspc_pix_mt = [];
    out.tcspc_mt_centers_ns = [];
    out.tcspc_mt_edges_ns = [];
    out.tcspc_mt_width_ns = [];
    out.tcspc_mt_edges_nativebin = [];
    if ~keepPhotonLists
        out = stripPhotonLists(out);
    end
    return;
end

nx = double(out.head.ImgHdr_PixX);
ny = double(out.head.ImgHdr_PixY);
dind = out.dind(:).';
nCh = numel(dind);

% Native TCSPC resolution in ns.
if isfield(out, 'tcspc_native_resolution_ns') && ~isempty(out.tcspc_native_resolution_ns)
    dtNativeNs = double(out.tcspc_native_resolution_ns);
elseif isfield(out.head, 'MeasDesc_Resolution_Original') && ~isempty(out.head.MeasDesc_Resolution_Original)
    dtNativeNs = double(out.head.MeasDesc_Resolution_Original) * 1e9;
else
    dtNativeNs = double(out.head.MeasDesc_Resolution) * 1e9;
end

pulsePeriodNs = double(out.head.MeasDesc_GlobalResolution) * 1e9;
NgateNative = max(1, ceil(pulsePeriodNs / dtNativeNs) + 1);

% Need native per-photon TCSPC to preserve the finest recorded timing.
if isfield(out, 'im_tcspc_native') && ~isempty(out.im_tcspc_native)
    tcspcNative = double(out.im_tcspc_native(:));
else
    error(['Native per-photon TCSPC bins are not available. ', ...
           'Call with keepNativePhotonTcspc = true if you want multi-tau ', ...
           'compression while preserving the finest recorded resolution.']);
end

% Determine global peak bin to anchor full-resolution region.
if isempty(tcspcNative)
    peakIdx = 1;
else
    tcspcNative = max(1, min(NgateNative, round(tcspcNative)));
    peakCounts = accumarray(tcspcNative, 1, [NgateNative, 1], @sum, 0);
    [~, peakIdx] = max(peakCounts);
    peakIdx = max(1, min(NgateNative, peakIdx));
end

[edgesNativeBin, centersNs, edgesNs, widthsNs, native2mt] = ...
    buildMultiTauAxis(NgateNative, dtNativeNs, binsPerOctave, peakIdx);

nMt = numel(centersNs);
mtIdx = native2mt(max(1, min(NgateNative, round(tcspcNative))));

imLine = double(out.im_line(:));
imCol  = double(out.im_col(:));
imChan = uint8(out.im_chan(:));

% Safety mask in case any photon coordinates drifted.
valid = (imLine >= 1) & (imLine <= nx) & ...
        (imCol  >= 1) & (imCol  <= ny) & ...
        (mtIdx  >= 1) & (mtIdx  <= nMt);

imLine = imLine(valid);
imCol  = imCol(valid);
imChan = imChan(valid);
mtIdx  = mtIdx(valid);

if keepPhotonLists
    out.im_tcspc_mt = uint16(mtIdx);
else
    out.im_tcspc_mt = uint16([]);
end

% Build compressed TCSPC cube channel by channel.
tcspc_pix_mt = zeros(nx, ny, nMt, nCh, 'uint16');

for ch = 1:nCh
    ind = (imChan == dind(ch));
    if ~any(ind)
        continue;
    end

    linIdx = sub2ind([nx, ny, nMt], imLine(ind), imCol(ind), double(mtIdx(ind)));
    counts = accumarray(linIdx, 1, [nx * ny * nMt, 1], @sum, 0);

    if max(counts) > double(intmax('uint16'))
        warning('PTU_MultiFrameScanReadFast_multiTauTCSPC:uint16clip', ...
            'tcspc_pix_mt exceeds uint16 range in channel %d; clipping.', ch);
        counts = min(counts, double(intmax('uint16')));
    end

    tcspc_pix_mt(:,:,:,ch) = reshape(uint16(counts), [nx, ny, nMt]);
end

out.tcspc_pix_mt = tcspc_pix_mt;
out.tcspc_mt_centers_ns = centersNs(:);
out.tcspc_mt_edges_ns = edgesNs(:);
out.tcspc_mt_width_ns = widthsNs(:);
out.tcspc_mt_edges_nativebin = edgesNativeBin(:);
out.head.TCSPC_MultiTau_BinsPerOctave = binsPerOctave;
out.head.TCSPC_MultiTau_Nbins = nMt;
out.head.TCSPC_MultiTau_IsCompressed = true;
out.head.TCSPC_MultiTau_NativeResolution_ns = dtNativeNs;
out.head.TCSPC_MultiTau_PulsePeriod_ns = pulsePeriodNs;
out.head.TCSPC_MultiTau_PeakBin = peakIdx;

% Dense linear cube intentionally omitted to avoid huge files.
out.tcspc_pix = [];
out.tcspc_pix_is_multitau = false;

% Useful compact global decays for plotting or ROI defaults.
out.tcspc_mt_global_total = squeeze(sum(sum(sum(tcspc_pix_mt, 4), 1), 2));
out.tcspc_mt_global_by_channel = squeeze(sum(sum(tcspc_pix_mt, 1), 2));

if ~keepPhotonLists
    out = stripPhotonLists(out);
end
end


function out = stripPhotonLists(out)
% Remove large photon-level arrays after tcspc_pix_mt has been built.
fields = {'im_sync','im_tcspc','im_tcspc_native','im_tcspc_mt','im_chan','im_line','im_col','im_frame','time'};
for k = 1:numel(fields)
    if isfield(out, fields{k})
        out.(fields{k}) = [];
    end
end
end


function [edgesNativeBin, centersNs, edgesNs, widthsNs, native2mt] = buildMultiTauAxis(NgateNative, dtNs, binsPerOctave, peakIdx)
% Build a multi-tau/log-binned time axis in native-bin units.
% Full resolution is kept up to the peak bin, then multitau after.
%
% Native bins are 1-based. The grouped multi-tau bins cover intervals
% [edgesNativeBin(i)+1, edgesNativeBin(i+1)] in native-bin indexing.

binsPerOctave = max(2, round(binsPerOctave));
if nargin < 4 || isempty(peakIdx)
    peakIdx = 1;
end
peakIdx = max(1, min(NgateNative, round(peakIdx)));

% Linear edges up to the peak bin.
edgesNativeBin = 0:peakIdx;

% Multitau after the peak bin.
remain = NgateNative - peakIdx;
if remain > 0
    widths = [];
    level = 0;
    covered = 0;
    while covered < remain
        w = 2^level;
        add = repmat(w, 1, binsPerOctave);
        widths = [widths add]; %#ok<AGROW>
        covered = sum(widths);
        level = level + 1;
    end
    edgesPost = peakIdx + cumsum(widths);
    edgesPost(edgesPost > NgateNative) = NgateNative;
    edgesNativeBin = unique([edgesNativeBin, edgesPost], 'stable');
end
if edgesNativeBin(end) < NgateNative
    edgesNativeBin(end+1) = NgateNative;
end

nMt = numel(edgesNativeBin) - 1;
widthsBins = diff(edgesNativeBin);

% Time edges in ns: native interval k occupies [(k-1)dt, k dt].
edgesNs = edgesNativeBin(:) * dtNs;
widthsNs = widthsBins(:) * dtNs;

% Arithmetic center of the grouped native-bin centers.
a = edgesNativeBin(1:end-1) + 1;
b = edgesNativeBin(2:end);
centersNs = (((double(a) + double(b)) / 2) - 0.5) * dtNs;
centersNs = centersNs(:);

native2mt = zeros(NgateNative, 1, 'uint16');
for i = 1:nMt
    native2mt((edgesNativeBin(i)+1):edgesNativeBin(i+1)) = uint16(i);
end
end
