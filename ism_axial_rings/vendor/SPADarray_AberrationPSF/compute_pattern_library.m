function lib = compute_pattern_library(varargin)
%COMPUTE_PATTERN_LIBRARY  Build ISM + SPAD pattern libraries for all orientations.
%
%   lib = compute_pattern_library('Name',Value,...)
%
%   For each (alpha,beta) on a uniform spherical grid, computes:
%     1) the ISM BUCKET image within the 1.25 AU pinhole (RADIAL / p-wave arm)
%     2) the configured SPAD pattern (AZIMUTHAL / s-wave arm)
%     3) the bucket count, SPAD counts, and f_r
%
%   DETECTOR ASSIGNMENT (ideal S-plate):
%     radial  (p-wave)  branch -> ISM bucket (1.25 AU circular pinhole)
%     azimuthal (s-wave) branch -> configured SPAD array
%
%   ORIENTATION GRID:
%     Default is a Fibonacci spiral on the upper hemisphere giving
%     near-uniform coverage.  Alternatively, 'grid_type','patgen' replicates
%     the PatternGeneration.m sin(alpha)-weighted grid.
%
%   All optics Name/Value pairs are forwarded to dipole_tomography.m.

% ---- separate library params ----
addpath 'D:\MATLAB\server\Luminosa\3D_Orient_simult'

lib_keys = {'Norient','al_res','be_res','grid_type','pic','detector_layout','detXY','detector_shape','detector_fill'};
[la, oa] = local_split(varargin, lib_keys);

lp = inputParser; lp.KeepUnmatched = true;
addParameter(lp,'Norient',   200);       % # orientations (Fibonacci)
addParameter(lp,'al_res',    10);        % deg (patgen mode)
addParameter(lp,'be_res',    10);        % deg (patgen mode)
addParameter(lp,'grid_type', 'fibonacci'); % 'fibonacci' | 'patgen'
addParameter(lp,'pic',       true);
addParameter(lp,'detector_layout', 'honeycomb23');
addParameter(lp,'detXY', []);
addParameter(lp,'detector_shape', 'hex');
addParameter(lp,'detector_fill', 1.0);
parse(lp, la{:});
L = lp.Results;

% ---- one reference call to bootstrap optics ----
ref = dipole_tomography(oa{:}, 'al',pi/2, 'be',0);
S   = ref.params;

if isempty(L.detXY)
    [spadXY, spadGrid, spadLayout] = detectorLayout(L.detector_layout, S.spad_pitch);
else
    spadXY = double(L.detXY);
    spadGrid = detectorIndexGridFromXY(spadXY, false, false, false);
    spadLayout = struct( ...
        'name', 'custom', ...
        'description', 'Custom detector coordinates', ...
        'nDet', size(spadXY,1), ...
        'detPitch', S.spad_pitch, ...
        'gridSize', size(spadGrid), ...
        'idxGrid', spadGrid);
end
nSpad = size(spadXY, 1);
detectorShape = lower(char(L.detector_shape));
detectorFill = L.detector_fill;
switch detectorShape
    case {'hex','hexagon','hexagonal','honeycomb'}
        detectorRadius = detectorFill * S.spad_pitch / sqrt(3);
        detectorHalfWidth = sqrt(3)/2 * detectorRadius;
    case {'square','rect','rectangle'}
        detectorRadius = [];
        detectorHalfWidth = 0.5 * detectorFill * S.spad_pitch;
    otherwise
        error('Unknown detector_shape: %s', L.detector_shape);
end

% ---- build orientation grid ----
switch lower(L.grid_type)
    case 'fibonacci'
        [theta, phi] = fibonacci_hemisphere(L.Norient);
    case 'patgen'
        [theta, phi] = patgen_grid(L.al_res, L.be_res);
    otherwise
        error('Unknown grid_type: %s', L.grid_type);
end
N = numel(theta);
fprintf('Orientation grid: %d (al,be) pairs (%s)\n', N, L.grid_type);

% =====================================================================
% PRECOMPUTE: emission radial profiles (orientation-INDEPENDENT)
% =====================================================================
I_pc    = ref.em.I_pc;    % radial arm, in-plane (J1)
I_v     = ref.em.I_v;     % radial arm, z-dipole (J0)
I_ps    = ref.em.I_ps;    % azimuthal arm, in-plane (J1)
rho_img = ref.em.rho;
AU_img  = ref.AU_image;

% image-plane grid (covers both pinhole and SPAD)
pinR     = 0.5 * S.pinhole_AU * AU_img;
half_det = max(0.8*S.pinhole_AU*AU_img, max(abs(spadXY(:))) + 1.2*detectorHalfWidth);
N1d      = max(2*ceil(half_det/AU_img*80)+1, 161);
xv       = linspace(-half_det, half_det, N1d);
dx       = xv(2)-xv(1);
[XX, YY] = ndgrid(xv, xv);
RR  = sqrt(XX.^2 + YY.^2);
PSI = atan2(YY, XX);

% interpolate profiles to 2D grid (done once)
Ipc_2d = reshape(interp1(rho_img, I_pc.', RR(:),'pchip',0), size(RR));
Iv_2d  = reshape(interp1(rho_img, I_v.',  RR(:),'pchip',0), size(RR));
Ips_2d = reshape(interp1(rho_img, I_ps.', RR(:),'pchip',0), size(RR));

% detector masks
mask_pin = RR <= pinR;
spad_mask = false([size(RR), nSpad]);
for k = 1:nSpad
    cx = spadXY(k,1);
    cy = spadXY(k,2);
    switch detectorShape
        case {'hex','hexagon','hexagonal','honeycomb'}
            spad_mask(:,:,k) = detectorHexMask(XX-cx, YY-cy, detectorRadius);
        case {'square','rect','rectangle'}
            spad_mask(:,:,k) = (abs(XX-cx)<=detectorHalfWidth) & ...
                               (abs(YY-cy)<=detectorHalfWidth);
    end
end

% ISM pinhole grid (small, for display)
Npin = 31;
pin_xv = linspace(-pinR, pinR, Npin);
[PXX, PYY] = ndgrid(pin_xv, pin_xv);
PRR = sqrt(PXX.^2 + PYY.^2);
PPSI = atan2(PYY, PXX);
pin_mask_small = PRR <= pinR;
Ipc_pin = reshape(interp1(rho_img, I_pc.', PRR(:),'pchip',0), size(PRR));
Iv_pin  = reshape(interp1(rho_img, I_v.',  PRR(:),'pchip',0), size(PRR));

% =====================================================================
% LOOP OVER ORIENTATIONS
% =====================================================================
spad_lib   = nan(size(spadGrid,1), size(spadGrid,2), N);
spad_counts = zeros(nSpad, N);
bucket     = zeros(1, N);
fr         = zeros(1, N);
ISM_images = zeros(Npin, Npin, N);

tic;
for j = 1:N
    al = theta(j);   be = phi(j);
    sa = sin(al);    ca = cos(al);

    % ---- radial arm (p-wave -> ISM bucket) ----
    U_r = 1i*sa*cos(PSI-be).*Ipc_2d - ca.*Iv_2d;
    Ir  = abs(U_r).^2;

    % ---- azimuthal arm (s-wave -> configured SPAD) ----
    U_p = sa*sin(PSI-be).*Ips_2d;
    Ip  = abs(U_p).^2;

    % counts
    Itot = (sum(Ir(:)) + sum(Ip(:)))*dx^2;
    if Itot < eps, Itot = 1; end
    nf = S.N_em / Itot;

    bucket(j) = sum(Ir(mask_pin))*dx^2 * nf;
    sv = squeeze(sum(sum(bsxfun(@times, spad_mask, Ip), 1), 2))*dx^2 * nf;
    spad_counts(:,j) = sv;
    spadMap = nan(size(spadGrid));
    validGrid = isfinite(spadGrid) & spadGrid >= 1 & spadGrid <= nSpad;
    spadMap(validGrid) = sv(spadGrid(validGrid));
    spad_lib(:,:,j) = spadMap;
    fr(j) = bucket(j) / max(bucket(j) + sum(sv), eps);

    % ISM image within pinhole
    U_r_pin = 1i*sa*cos(PPSI-be).*Ipc_pin - ca.*Iv_pin;
    ISM_images(:,:,j) = abs(U_r_pin).^2 .* pin_mask_small;
end
fprintf('Built %d-orientation library in %.1f s.\n', N, toc);

% =====================================================================
% PACK
% =====================================================================
lib.theta      = theta;
lib.phi        = phi;
lib.spad       = spad_lib;        % detector display grid x N  (azimuthal / s-wave arm)
lib.spadCounts = spad_counts;     % nDetector x N raw detector counts
lib.detXY      = spadXY;
lib.detectorIndexGrid = spadGrid;
lib.detectorLayout = spadLayout;
lib.detectorShape = detectorShape;
lib.detectorFill = detectorFill;
lib.bucket     = bucket;          % 1xN    (radial / p-wave arm, integrated)
lib.fr         = fr;              % 1xN
lib.ISM.images = ISM_images;      % NpinxNpinxN  (radial arm image in pinhole)
lib.ISM.x      = pin_xv;
lib.ISM.y      = pin_xv;
lib.ISM.pinR   = pinR;
lib.optics     = S;
lib.AU_image   = AU_img;
lib.grid_type  = L.grid_type;

if L.pic
    show_library(lib);
end
end


% =========================================================================
function [theta, phi] = fibonacci_hemisphere(N)
% Fibonacci spiral on the upper hemisphere: near-uniform area coverage.
% cos(alpha) is uniformly spaced from 1 (vertical) to 0 (in-plane).
golden = (1 + sqrt(5))/2;
theta = zeros(1, N);
phi   = zeros(1, N);
for i = 1:N
    theta(i) = acos(1 - (i-0.5)/N);
    phi(i)   = mod(2*pi*(i-1)/golden, 2*pi);
end
theta = min(theta, pi/2);
end


% =========================================================================
function [theta, phi] = patgen_grid(al_res, be_res)
theta = [];  phi = [];
for k = 90:-al_res:0
    al = k*pi/180;
    if k == 90
        jj = round(180/be_res);  dbe = pi/jj;
    elseif k == 0
        jj = 1;  dbe = 0;
    else
        jj = round(sin(al)*360/be_res);  dbe = 2*pi/jj;
    end
    theta = [theta, al*ones(1,jj)]; %#ok
    phi   = [phi,   dbe*(0:jj-1)];  %#ok
end
if isempty(theta) || theta(end) ~= 0
    theta(end+1) = 0;  phi(end+1) = 0;
end
end


% =========================================================================
function [a_lib, a_other] = local_split(args, keys)
a_lib = {};  a_other = {};
i = 1;
while i <= numel(args)
    if any(strcmpi(args{i}, keys))
        a_lib(end+1:end+2) = args(i:i+1);
    else
        a_other(end+1:end+2) = args(i:i+1);
    end
    i = i+2;
end
end


% =========================================================================
function show_library(lib)
% Display ISM and SPAD patterns at representative (alpha,beta) pairs.
% Rows = alpha, columns = beta.  This shows the actual orientation-
% dependent structure that the azimuth-averaged display was hiding.

al_show = [0 15 30 45 60 75 90] * pi/180;
be_show = [0 45 90 135 225 315] * pi/180;
n_al = numel(al_show);
n_be = numel(be_show);

% ---- Figure 1: ISM bucket images (radial / p-wave arm) ----
figure('Color','w','Name','ISM bucket images (radial/p-wave arm)');
for ia = 1:n_al
    for ib = 1:n_be
        if al_show(ia) < 1e-3 && ib > 1, continue; end
        j = find_nearest(lib.theta, lib.phi, al_show(ia), be_show(ib));
        subplot(n_al, n_be, (ia-1)*n_be + ib);
        img = lib.ISM.images(:,:,j);
        imagesc(lib.ISM.x*1e6, lib.ISM.y*1e6, img.');
        axis image off; colormap(hot);
        if ia == 1 && ib == 1
            title(sprintf('\\alpha=%.0f^\\circ', lib.theta(j)*180/pi),'FontSize',7);
        else
            title(sprintf('%.0f/%.0f', lib.theta(j)*180/pi, lib.phi(j)*180/pi),'FontSize',7);
        end
    end
end
sgtitle({'ISM bucket (radial/p-wave arm) within 1.25 AU pinhole', ...
         'rows = \alpha,  cols = \beta'});

% ---- Figure 2: SPAD detector map (azimuthal / s-wave arm) ----
figSpad = figure('Color','w','Name','SPAD detector hex map (azimuthal/s-wave arm)', ...
    'Position', [80 80 1200 850]*0.85);
tl = tiledlayout(figSpad, n_al, n_be, 'TileSpacing', 'compact', 'Padding', 'compact');

% Use a GLOBAL color scale so you can see how total brightness changes with alpha.
finiteSpad = lib.spadCounts(isfinite(lib.spadCounts));
if isempty(finiteSpad)
    all_max = 1;
else
    all_max = max(finiteSpad);
end
if all_max <= 0
    all_max = 1;
end

for ia = 1:n_al
    for ib = 1:n_be
        ax = nexttile(tl, (ia-1)*n_be + ib);
        if al_show(ia) < 1e-3 && ib > 1
            axis(ax, 'off');
            continue;
        end
        j = find_nearest(lib.theta, lib.phi, al_show(ia), be_show(ib));
        plotDetectorHexMap(lib.detXY, lib.spadCounts(:,j), ...
            'Parent', ax, ...
            'CLim', [0 all_max*0.8], ...
            'CellScale', 1.01*lib.detectorFill, ...
            'EdgeColor', [0.50 0.50 0.50], ...
            'LineWidth', 0.6);
        colormap(ax, hot);
        if ia == 1 && ib == 1
            title(sprintf('\\alpha=%.0f^\\circ', lib.theta(j)*180/pi),'FontSize',7);
        else
            title(sprintf('%.0f/%.0f', lib.theta(j)*180/pi, lib.phi(j)*180/pi),'FontSize',7);
        end
    end
end
title(tl, {sprintf('SPAD %s hex-cell map (azimuthal/s-wave arm) - GLOBAL color scale', lib.detectorLayout.name), ...
         '\alpha=0 is dark; lobes rotate with \beta'});

% ---- Figure 3: f_r on the hemisphere ----
figure('Color','w','Name','f_r on orientation hemisphere');
subplot(1,2,1);
scatter(lib.theta*180/pi, lib.fr, 30, lib.phi*180/pi, 'filled');
xlabel('\alpha (deg)'); ylabel('f_r');
cb = colorbar; cb.Label.String = '\beta (deg)';
title(sprintf('f_r (%d orientations)', numel(lib.theta))); grid on;

subplot(1,2,2);
px = sin(lib.theta).*cos(lib.phi);
py = sin(lib.theta).*sin(lib.phi);
pz = cos(lib.theta);
scatter3(px, py, pz, 40, lib.fr, 'filled');
xlabel('p_x'); ylabel('p_y'); zlabel('p_z');
cb = colorbar; cb.Label.String = 'f_r'; axis equal;
title('f_r on the hemisphere');
end


% =========================================================================
function j = find_nearest(theta, phi, al, be)
% Great-circle distance on the unit sphere (correct for the hemisphere).
cos_d = sin(al).*sin(theta).*cos(be-phi) + cos(al).*cos(theta);
cos_d = max(-1, min(1, cos_d));
dists = acos(cos_d);
[~,j] = min(dists);
end
