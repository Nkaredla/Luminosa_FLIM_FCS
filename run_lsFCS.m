% run_lsFCS
% Batch lsFCS correlation over subfolders with optional diffusion fitting.
%
% Each subfolder under rootDir is expected to contain exactly one .ptu file.
% The output for each file is saved by lsFCS as *_lsFCS.mat next to the .ptu.

% === User settings ===
rootDir   = 'D:\Luminosa\Data\270226\POPC_DiO_10pc_NTA_B7H6_Atto488_21\'; % <-- change this
cnum      = 2;     % number of pulses in one sync (set to your experiment)
maxtime   = 1;    % max correlation time (seconds)
timegates = [];    % [] = auto-detect; or Nx2 array [start len] per pulse
overwrite = true; % set true to re-run even if output exists
scan_length = 5; %um - match with the measurement
nx = 50; % number of pixels in line (default)

% Prompt once before processing the directory
ansFit = input('Run diffusion fit for components 1 and 4? y/n [n]: ', 's');
doFit = ~isempty(ansFit) && lower(ansFit(1)) == 'y';

% Fit settings (used only if doFit == true)
fitCompIdx = [1 4];
p0   = [5e-3  0.30]; % [tauD (s), w0 (um)] w0 is Gaussian sigma
pmin = [1e-6  0.05];
pmax = [10    5.00];
tol  = 1e-5;
steps = 400;

%%
% === Batch processing ===
dirInfo = dir(rootDir);
for k = 1:numel(dirInfo)
    if ~dirInfo(k).isdir
        continue;
    end
    name = dirInfo(k).name;
    if strcmp(name,'.') || strcmp(name,'..')
        continue;
    end

    subDir = fullfile(rootDir, name);
    ptuFiles = dir(fullfile(subDir, '*.ptu'));

    if isempty(ptuFiles)
        fprintf('No PTU found in %s\n', subDir);
        continue;
    end
    if numel(ptuFiles) > 1
        fprintf('Multiple PTUs found in %s (skipping)\n', subDir);
        continue;
    end

    ptuPath = fullfile(subDir, ptuFiles(1).name);
    outMat = [ptuPath(1:end-4), '_lsFCS.mat'];
    if ~overwrite && exist(outMat, 'file')
        fprintf('Output exists, skipping: %s\n', outMat);
        continue;
    end

    fprintf('Processing %s\n', ptuPath);
    try
        res = lsFCS(ptuPath, cnum, maxtime, timegates);
    catch ME
        fprintf('Failed %s: %s\n', ptuPath, ME.message);
        continue;
    end

    if doFit
        try
            [G, ~, ~, ~, t, xxi] = lsCrossRead(res, 25);
            dnum = size(res.tcspc, 2);
            cnum_here = size(res.tcspc, 3);
            pixel_um = scan_length/nx; %um

            for ii = 1:numel(fitCompIdx)
                comp = fitCompIdx(ii);
                if comp < 1 || comp > dnum * cnum_here
                    warning('Component %d is out of range (max %d). Skipping.', comp, dnum * cnum_here);
                    continue;
                end

                compName = localComponentName(comp, dnum, cnum_here);
                data = G(:,:,comp);

                [pfit, ~] = Simplex(@lsFCS_diff1_fit, p0, pmin, pmax, tol, steps, t, xxi, data, pixel_um);
                [~, c, z] = lsFCS_diff1_fit(pfit, t, xxi, data, pixel_um);
                tauD = pfit(1);
                w0 = pfit(2); % sigma
                D = (w0.^2) ./ tauD; % um^2/s for sigma-defined w0

                xx = xxi * pixel_um;
                figure('Name', compName, 'NumberTitle', 'off');

                subplot(1,3,1);
                h = surf(xx, t, data);
                set(get(h,'Parent'),'YScale','Log');
                title([compName ' raw']);
                shading interp
                view(2)
                axis tight
                box on
                ylabel('lag time / s')
                xlabel('lag space / \mum')
                colorbar
                localSetDecadeTicks(gca, t)

                subplot(1,3,2);
                h = surf(xx, t, z);
                set(get(h,'Parent'),'YScale','Log');
                title([compName ' fit']);
                shading interp
                view(2)
                axis tight
                box on
                ylabel('lag time / s')
                xlabel('lag space / \mum')
                colorbar
                localSetDecadeTicks(gca, t)

                subplot(1,3,3);
                resid = (data - z) ./ max(abs(data), eps);
                h = surf(xx, t, resid);
                set(get(h,'Parent'),'YScale','Log');
                title([compName ' resid']);
                shading interp
                view(2)
                axis tight
                box on
                ylabel('lag time / s')
                xlabel('lag space / \mum')
                colorbar
                localSetDecadeTicks(gca, t)

                sgtitle(sprintf('%s  w0=%.3g um  D=%.3g um^2/s', compName, w0, D))
                print([ptuPath(1:end-4) '_' compName '_diff1_fit.png'], '-dpng', '-r300');
            end
        catch ME
            fprintf('Fit failed for %s: %s\n', ptuPath, ME.message);
        end
    end
end

function name = localComponentName(comp, dnum, cnum)
    det = ceil(comp / cnum);
    pulse = comp - (det-1) * cnum;
    if det > dnum
        det = dnum;
    end
    name = sprintf('D%dL%d', det, pulse);
end

function localSetDecadeTicks(ax, t)
    if isempty(t)
        return;
    end
    tpos = t(isfinite(t) & t > 0);
    if isempty(tpos)
        return;
    end
    tmin = min(tpos);
    tmax = max(tpos);
    nmin = floor(log10(tmin));
    nmax = ceil(log10(tmax));
    ticks = 10.^(nmin:nmax);
    ticks = ticks(ticks >= tmin & ticks <= tmax);
    if numel(ticks) < 3
        ticks = logspace(log10(tmin), log10(tmax), 3);
    end
    set(ax, 'YTick', unique(ticks));
end
