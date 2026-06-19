function diag = ism_reassignment_scale_diagnostic(results, opts)
%ISM_REASSIGNMENT_SCALE_DIAGNOSTIC  Test the ISM pixel-reassignment 1/2 factor.
%
% Compares the MEAd SURED per-detector reassignment shifts (results.shiftsToCenter)
% against the KNOWN physical detector-array geometry built from the hex layout.
%
% Why not results.detectorPositions?
%   In run_ism_reconstruction_from_ptu, detectorPositions == shiftsToTempRef,
%   i.e. it is itself a measured (already-halved) image displacement. Fitting
%   shifts against it would be circular and would *always* return a~1, hiding
%   exactly the deviation you want to see. This routine therefore builds the
%   geometry independently from the detector array layout.
%
% Ideal ISM reassignment:  shift_k = (1/2) * g_k
%   g_k = physical detector offset from the centre element, in OBJECT pixels.
%   This routine fits  shift = a * R * g + t  (similarity) and an unconstrained
%   2x2 map, and reports a. a ~ 0.5 => ideal reassignment.
%   a < 0.5 => under-reassignment (PSF mismatch / aberration).
%   a far from 0.5 with isotropic SVD => magnification / pitch error.
%   anisotropic SVD (s1/s2 >> 1) => astigmatism or scan-axis miscalibration.
%
% USAGE
%   diag = ism_reassignment_scale_diagnostic(results, opts)
%
% results : output struct of run_ism_reconstruction_from_ptu
%
% opts fields (all optional)
%   .detectorPitchObjectPixels : nearest-neighbour detector spacing expressed
%       in OBJECT pixels (FULL geometry, before the 1/2). REQUIRED for an
%       absolute 1/2 test. Compute as  (physicalPitch / totalMagnification) /
%       objectPixelSize. If omitted, geometry is left in detector-pitch units
%       and the reported scale is (object px per detector pitch); the absolute
%       1/2 test is then skipped.
%   .detectorPitchPhysical : alternative input (same units as objectPixelSize
%       after dividing by magnification).
%   .magnification         : total sample->detector magnification (with above).
%   .objectPixelSize       : object pixel size; default results.paramsUsed.pixelSize.
%   .channelToLayoutIndex  : permutation mapping channel k -> hex-layout slot k.
%       Default identity, i.e. the SAME ordering assumption your hex image-grid
%       plots already rely on. If that assumption is wrong, those plots are
%       wrong too, and this fit will be garbage -- verify it.
%   .useChannels           : logical mask or index list selecting detectors to fit.
%   .makePlot              : default true.
%   .plotFile              : optional PNG path.
%
% OUTPUT  diag struct
%   .scaleSimilarity    a in shift = a*R*g + t  (headline number; ~0.5 ideal)
%   .reassignmentFactor a expressed dimensionlessly (== scaleSimilarity when
%                       geometry was converted to object px; else NaN)
%   .rotationDeg        fitted array rotation relative to scan axes
%   .reflection         true if the convention reflects (det(A) < 0)
%   .singularValues     SVD of the unconstrained 2x2 map (anisotropy probe)
%   .anisotropy         s1/s2  (1.0 = isotropic)
%   .slopeThroughOrigin |shift| = slope*|g|, rotation-invariant scalar
%   .R2                 R^2 of that magnitude fit
%   .residualRMS        RMS of similarity-fit residual (object px)
%   .geometryUnits      'object px' or 'detector pitch'
%   .perChannel         table: channel, gMag, shiftMag, residual

if nargin < 2 || isempty(opts), opts = struct(); end
opts = setdef(opts, 'detectorPitchObjectPixels', NaN);
opts = setdef(opts, 'detectorPitchPhysical', NaN);
opts = setdef(opts, 'magnification', NaN);
opts = setdef(opts, 'objectPixelSize', NaN);
opts = setdef(opts, 'channelToLayoutIndex', []);
opts = setdef(opts, 'useChannels', []);
opts = setdef(opts, 'makePlot', true);
opts = setdef(opts, 'plotFile', '');

% ---- resolve pitch in object pixels ------------------------------------
pitchObjPx = opts.detectorPitchObjectPixels;
if isnan(pitchObjPx) && ~isnan(opts.detectorPitchPhysical) && ~isnan(opts.magnification)
    objPx = opts.objectPixelSize;
    if isnan(objPx)
        if isfield(results, 'paramsUsed') && isfield(results.paramsUsed, 'pixelSize')
            objPx = results.paramsUsed.pixelSize;
        else
            error('Provide opts.objectPixelSize (or results.paramsUsed.pixelSize).');
        end
    end
    pitchObjPx = (opts.detectorPitchPhysical / opts.magnification) / objPx;
end

haveAbsolute = ~isnan(pitchObjPx);
if haveAbsolute
    geomScale = pitchObjPx;
    geometryUnits = 'object px';
else
    geomScale = 1;                 % geometry stays in detector-pitch units
    geometryUnits = 'detector pitch';
    warning(['detectorPitchObjectPixels not supplied: geometry kept in pitch ', ...
        'units. The absolute 1/2 test is skipped; reported scale is object px ', ...
        'per detector pitch (expected ~0.5*pitchObjPx).']);
end

% ---- build INDEPENDENT physical geometry from the hex layout -----------
nDet = numel(results.channelIDs);
[layoutXY, useHex] = local_detector_layout(nDet);
if ~useHex
    warning(['No hex layout for nDet=%d; using a circular fallback that is ', ...
        'NOT real geometry. Supply a real layout for a meaningful fit.'], nDet);
end

map = opts.channelToLayoutIndex;
if isempty(map), map = (1:nDet).'; end
assert(numel(map) == nDet, 'channelToLayoutIndex length must equal nDet.');
layoutXY = layoutXY(map, :);

centerIdx = results.centerDetectorIndex;
layoutXY = layoutXY - layoutXY(centerIdx, :);     % centre on centre detector
geomXY   = layoutXY * geomScale;                  % full offset g_k, [x y]

% ---- measured shifts in the same y-up display frame --------------------
% shiftsToCenter is stored [dy dx]; convert to [x y] = [dx, -dy] (matches the
% existing save_detector_shift_vectors convention, so no spurious reflection).
shiftXY = [results.shiftsToCenter(:,2), -results.shiftsToCenter(:,1)];

% ---- detector subset ---------------------------------------------------
sel = true(nDet, 1);
if ~isempty(opts.useChannels)
    if islogical(opts.useChannels)
        sel = opts.useChannels(:);
    else
        sel = false(nDet,1); sel(opts.useChannels) = true;
    end
end
sel(centerIdx) = true;                            % keep centre as anchor
G = geomXY(sel, :);
S = shiftXY(sel, :);
chSel = results.channelIDs(sel);

% ---- fits --------------------------------------------------------------
[a, R, t, residRMS, res] = fit_similarity(G, S);  % shift = a*R*g + t
rotationDeg = atan2d(R(2,1), R(1,1));

% unconstrained 2x2 (centred) for anisotropy / reflection
Gc = G - mean(G,1); Sc = S - mean(S,1);
A  = Gc \ Sc;                                     % g_row * A ~ s_row
sv = svd(A);
anisotropy = sv(1) / max(sv(2), eps);
reflection = det(A) < 0;

% rotation-invariant magnitude fit through origin
gMag = sqrt(sum(geomXY.^2, 2));
sMag = sqrt(sum(shiftXY.^2, 2));
gm = gMag(sel); sm = sMag(sel);
slope = sum(gm .* sm) / max(sum(gm.^2), eps);
ssTot = sum((sm - mean(sm)).^2);
ssRes = sum((sm - slope*gm).^2);
R2 = 1 - ssRes / max(ssTot, eps);

% ---- pack --------------------------------------------------------------
diag = struct();
diag.scaleSimilarity   = a;
diag.reassignmentFactor = ternary(haveAbsolute, a, NaN);
diag.rotationDeg       = rotationDeg;
diag.reflection        = reflection;
diag.singularValues    = sv(:).';
diag.anisotropy        = anisotropy;
diag.slopeThroughOrigin = slope;
diag.R2                = R2;
diag.residualRMS       = residRMS;
diag.geometryUnits     = geometryUnits;
diag.geometryXY        = geomXY;
diag.shiftXY           = shiftXY;
diag.selection         = sel;

resFull = nan(nDet,1);
resFull(sel) = sqrt(sum(res.^2, 2));
diag.perChannel = table(results.channelIDs(:), gMag, sMag, resFull, ...
    'VariableNames', {'ChannelID','gMag','shiftMag','residual'});

% ---- report ------------------------------------------------------------
fprintf('\nISM reassignment-scale diagnostic\n');
fprintf('  geometry units      : %s\n', geometryUnits);
if haveAbsolute
    fprintf('  reassignment factor : %.4f   (ideal 0.5000)\n', a);
else
    fprintf('  scale (px / pitch)  : %.4f   (ideal 0.5 * pitchObjPx)\n', a);
end
fprintf('  array rotation      : %.2f deg\n', rotationDeg);
fprintf('  reflection in conv. : %d\n', reflection);
fprintf('  SVD (2x2)           : [%.4f %.4f]  anisotropy %.3f\n', sv(1), sv(2), anisotropy);
fprintf('  |shift| vs |g| slope: %.4f   R^2 %.4f\n', slope, R2);
fprintf('  similarity resid RMS: %.4f object px\n\n', residRMS);

% ---- plot --------------------------------------------------------------
if opts.makePlot
    plot_diagnostic(diag, haveAbsolute, char_or_empty(opts.plotFile));
end
end

% =========================================================================
function [a, R, t, residRMS, res] = fit_similarity(X, Y)
% Umeyama (1991) similarity fit: Y ~ a*R*X + t, rotation only (no reflection).
N   = size(X,1);
muX = mean(X,1); muY = mean(Y,1);
Xc  = X - muX;   Yc = Y - muY;
Sig = (Yc.' * Xc) / N;                 % 2x2
[U,D,V] = svd(Sig);
W = eye(2);
if det(U)*det(V) < 0, W(2,2) = -1; end % force proper rotation
R = U*W*V.';
varX = mean(sum(Xc.^2, 2));
a = trace(D*W) / max(varX, eps);
t = muY.' - a*R*muX.';
pred = (a*R*X.' + t).';
res  = Y - pred;
residRMS = sqrt(mean(sum(res.^2, 2)));
end

% =========================================================================
function plot_diagnostic(diag, haveAbsolute, plotFile)
geomXY = diag.geometryXY; shiftXY = diag.shiftXY; sel = diag.selection;
gMag = sqrt(sum(geomXY.^2,2)); sMag = sqrt(sum(shiftXY.^2,2));

vis = 'on'; if ~isempty(plotFile), vis = 'off'; end
fig = figure('Color','w','Visible',vis,'Name','ISM reassignment diagnostic', ...
    'Position',[100 100 1180 520]);
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% Panel 1: magnitude scatter + slope (rotation-invariant)
ax1 = nexttile(tl,1); hold(ax1,'on');
plot(ax1, gMag(sel), sMag(sel), 'o', 'MarkerFaceColor',[0.1 0.35 0.85], ...
    'MarkerEdgeColor','k', 'MarkerSize',6);
plot(ax1, gMag(~sel), sMag(~sel), 'o', 'Color',[0.6 0.6 0.6], 'MarkerSize',5);
xl = [0, max(gMag)*1.05];
plot(ax1, xl, diag.slopeThroughOrigin*xl, '-', 'Color',[0.1 0.35 0.85], 'LineWidth',1.4);
if haveAbsolute
    plot(ax1, xl, 0.5*xl, '--', 'Color',[0.85 0.1 0.1], 'LineWidth',1.4);
    legend(ax1, {'fit detectors','excluded',sprintf('fit slope %.3f',diag.slopeThroughOrigin), ...
        'ideal 0.5'}, 'Location','northwest');
else
    legend(ax1, {'fit detectors','excluded',sprintf('fit slope %.3f',diag.slopeThroughOrigin)}, ...
        'Location','northwest');
end
xlabel(ax1,'|g| physical detector offset'); ylabel(ax1,'|shift| measured reassignment');
title(ax1, sprintf('Magnitude fit  (R^2 = %.3f)', diag.R2));
grid(ax1,'on'); axis(ax1,'tight'); box(ax1,'on');

% Panel 2: 2D vector field, measured vs similarity prediction
ax2 = nexttile(tl,2); hold(ax2,'on');
a = diag.scaleSimilarity; th = deg2rad(diag.rotationDeg);
R = [cos(th) -sin(th); sin(th) cos(th)];
pred = (a*R*geomXY.').';                 % similarity prediction (t ~ 0, centred)
plot(ax2, geomXY(:,1), geomXY(:,2), 'o', 'Color',[0.55 0.55 0.55], ...
    'MarkerFaceColor',[0.88 0.88 0.88], 'MarkerSize',6);
hM = quiver(ax2, geomXY(:,1), geomXY(:,2), shiftXY(:,1)-0, shiftXY(:,2)-0, 0, ...
    'Color',[0.1 0.35 0.85], 'LineWidth',1.0, 'MaxHeadSize',0.25);
hP = quiver(ax2, geomXY(:,1), geomXY(:,2), pred(:,1), pred(:,2), 0, ...
    'Color',[0.85 0.1 0.1], 'LineWidth',1.0, 'LineStyle','--', 'MaxHeadSize',0.25);
plot(ax2, 0, 0, 'rp', 'MarkerSize',14, 'LineWidth',2, 'MarkerFaceColor','r');
axis(ax2,'equal'); grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'\Delta x [object px]'); ylabel(ax2,'\Delta y [object px]');
legend(ax2, [hM hP], {'measured shift','similarity fit'}, 'Location','best');
if haveAbsolute
    ttl = sprintf('factor a = %.3f (ideal 0.5),  rot %.1f\\circ,  aniso %.2f,  refl %d', ...
        a, diag.rotationDeg, diag.anisotropy, diag.reflection);
else
    ttl = sprintf('scale %.3f px/pitch,  rot %.1f\\circ,  aniso %.2f,  refl %d', ...
        a, diag.rotationDeg, diag.anisotropy, diag.reflection);
end
title(ax2, ttl);

title(tl, 'ISM pixel-reassignment scale: measured shift vs known geometry', ...
    'Interpreter','none');

if ~isempty(plotFile)
    print(fig, plotFile, '-dpng', '-r300'); close(fig);
end
end

% =========================================================================
function [detXY, useHex] = local_detector_layout(nDet)
% Same convention as run_ism_reconstruction_from_ptu (hex 23 / circular fallback).
if nDet == 23
    rowCounts = [5 4 5 4 5];
    detPitch = 1;
    yPitch = sqrt(3)/2 * detPitch;
    yGrid = ((numel(rowCounts)-1)/2:-1:-(numel(rowCounts)-1)/2) * yPitch;
    detXY = zeros(sum(rowCounts), 2);
    k = 0;
    for r = 1:numel(rowCounts)
        nThisRow = rowCounts(r);
        if mod(nThisRow,2) == 1
            xRow = (-(nThisRow-1)/2:(nThisRow-1)/2) * detPitch;
        else
            xRow = ((-nThisRow/2):(nThisRow/2-1)) * detPitch + detPitch/2;
        end
        for c = 1:nThisRow
            k = k + 1;
            detXY(k,:) = [xRow(c), yGrid(r)];
        end
    end
    useHex = true;
else
    theta = linspace(0, 2*pi, nDet+1).'; theta(end) = [];
    detXY = [cos(theta), sin(theta)];
    useHex = false;
end
end

% =========================================================================
function s = setdef(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
function t = char_or_empty(t)
if isstring(t), t = char(t); end
if isempty(t), t = ''; end
end