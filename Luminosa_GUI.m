function Luminosa_GUI
% Luminosa_GUI
% Simple GUI for FLIM / ISM-FLIM processing and FCS analysis.

    app = struct();
    app.ptuOut = [];
    app.ismRes = [];
    app.flim = [];
    app.fcs = [];
    app.roi = [];
    app.lastFile = '';

    app.fig = uifigure('Name','Luminosa FLIM / ISM-FLIM / FCS','Position',[100 100 1200 720]);
    app.grid = uigridlayout(app.fig, [1 1]);
    app.tg = uitabgroup(app.grid);

    tabFLIM = uitab(app.tg, 'Title', 'FLIM / ISM-FLIM');
    tabFCS  = uitab(app.tg, 'Title', 'FCS');

    buildFLIMTab(tabFLIM);
    buildFCSTab(tabFCS);

    % ----------------------- UI builders -----------------------
    function buildFLIMTab(parent)
        g = uigridlayout(parent, [1 2]);
        g.ColumnWidth = {320, '1x'};
        g.RowHeight = {'1x'};

        ctl = uipanel(g, 'Title', 'Controls');
        axp = uipanel(g, 'Title', 'Display');

        ctlGrid = uigridlayout(ctl, [22 1]);
        ctlGrid.RowHeight = {26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,'1x'};

        app.btnLoad = uibutton(ctlGrid, 'Text', 'Load PTU', 'ButtonPushedFcn', @onLoadPTU);
        app.chkStoreTCSPC = uicheckbox(ctlGrid, 'Text', 'Store TCSPC cube (for ROI)', 'Value', true);
        app.chkUseGPU = uicheckbox(ctlGrid, 'Text', 'Use GPU (if available)', 'Value', false);

        app.btnQuickFLIM = uibutton(ctlGrid, 'Text', 'Run FLIM (quick)', 'ButtonPushedFcn', @onQuickFLIM);
        app.btnISMFLIM = uibutton(ctlGrid, 'Text', 'Run ISM-FLIM', 'ButtonPushedFcn', @onISMFLIM);

        app.btnShowIntensity = uibutton(ctlGrid, 'Text', 'Show Intensity', 'ButtonPushedFcn', @onShowIntensity);
        app.btnShowTau = uibutton(ctlGrid, 'Text', 'Show TauMean', 'ButtonPushedFcn', @onShowTauMean);

        app.btnSelectROI = uibutton(ctlGrid, 'Text', 'Select ROI', 'ButtonPushedFcn', @onSelectROI);
        uilabel(ctlGrid, 'Text', 'TCSPC log bins');
        app.txtBins = uieditfield(ctlGrid, 'numeric', 'Value', 64, 'Limits', [8 512], ...
            'LowerLimitInclusive','on','UpperLimitInclusive','on');
        app.btnShowTCSPC = uibutton(ctlGrid, 'Text', 'Show ROI TCSPC', 'ButtonPushedFcn', @onShowTCSPC);

        app.btnSaveMAT = uibutton(ctlGrid, 'Text', 'Save MAT (PTU/FLIM)', 'ButtonPushedFcn', @onSaveMAT);
        app.btnSavePNG = uibutton(ctlGrid, 'Text', 'Save PNG (current)', 'ButtonPushedFcn', @onSavePNG);

        app.status = uilabel(ctlGrid, 'Text', 'Ready.');
        app.status.WordWrap = 'on';

        axGrid = uigridlayout(axp, [2 1]);
        axGrid.RowHeight = {'2x','1x'};
        app.axImage = uiaxes(axGrid);
        title(app.axImage, 'Image');
        app.axTCSPC = uiaxes(axGrid);
        title(app.axTCSPC, 'ROI TCSPC (log bins)');

        % GPU availability
        if ~gpuIsAvailable()
            app.chkUseGPU.Enable = 'off';
            app.chkUseGPU.Value = false;
        end
    end

    function buildFCSTab(parent)
        g = uigridlayout(parent, [1 2]);
        g.ColumnWidth = {320, '1x'};
        g.RowHeight = {'1x'};

        ctl = uipanel(g, 'Title', 'Controls');
        axp = uipanel(g, 'Title', 'Display');

        ctlGrid = uigridlayout(ctl, [18 1]);
        ctlGrid.RowHeight = {26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,'1x'};

        app.btnLoadFCS = uibutton(ctlGrid, 'Text', 'Load PTU (FCS)', 'ButtonPushedFcn', @onLoadPTU);
        app.btnUseLoaded = uibutton(ctlGrid, 'Text', 'Use Loaded PTU', 'ButtonPushedFcn', @onUseLoadedPTU);

        uilabel(ctlGrid, 'Text', 'cnum');
        app.editCnum = uieditfield(ctlGrid, 'numeric', 'Value', 1, 'Limits', [1 16]);
        uilabel(ctlGrid, 'Text', 'maxtime (s)');
        app.editMaxtime = uieditfield(ctlGrid, 'numeric', 'Value', 10, 'Limits', [0.01 1e3]);
        uilabel(ctlGrid, 'Text', 'timegates (MATLAB)');
        app.editTimegates = uieditfield(ctlGrid, 'text', 'Value', '[]');

        app.btnPointFCS = uibutton(ctlGrid, 'Text', 'Run Point FCS', 'ButtonPushedFcn', @onPointFCS);
        app.btnLineFCS = uibutton(ctlGrid, 'Text', 'Run LineScan FCS', 'ButtonPushedFcn', @onLineFCS);

        app.btnSaveFCS = uibutton(ctlGrid, 'Text', 'Save FCS MAT', 'ButtonPushedFcn', @onSaveFCS);

        app.statusFCS = uilabel(ctlGrid, 'Text', 'Ready.');
        app.statusFCS.WordWrap = 'on';

        axGrid = uigridlayout(axp, [1 1]);
        app.axFCS = uiaxes(axGrid);
        title(app.axFCS, 'FCS');
    end

    % ----------------------- Callbacks -----------------------
    function onLoadPTU(~, ~)
        [f, p] = uigetfile('*.ptu', 'Select PTU file');
        if isequal(f, 0)
            return;
        end
        app.lastFile = fullfile(p, f);
        addStatus('Loading PTU...');
        drawnow;

        storeTCSPC = app.chkStoreTCSPC.Value;
        try
            app.ptuOut = PTU_MultiFrameScanReadFast(app.lastFile, 1e6, storeTCSPC);
            addStatus('PTU loaded.');
            showIntensityFromPTU();
        catch ME
            addStatus(['Load failed: ' ME.message]);
        end
    end

    function onUseLoadedPTU(~, ~)
        if isempty(app.ptuOut)
            app.statusFCS.Text = 'No PTU loaded.';
            return;
        end
        app.statusFCS.Text = 'Using currently loaded PTU.';
    end

    function onQuickFLIM(~, ~)
        if isempty(app.ptuOut) || ~isfield(app.ptuOut, 'tcspc_pix') || isempty(app.ptuOut.tcspc_pix)
            addStatus('No TCSPC cube available. Enable "Store TCSPC cube" and reload.');
            return;
        end
        addStatus('Computing quick FLIM...');
        drawnow;

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        flim = quickFLIMFromTCSPC(app.ptuOut, useGPU);
        app.flim = flim;
        addStatus('Quick FLIM done.');
        showTauMean();
    end

    function onISMFLIM(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        addStatus('Running ISM reconstruction...');
        drawnow;
        params = struct();
        params.imageSource = 'tags';
        params.smoothSigma = 1;
        params.useWindow = true;
        params.normalizeImages = true;
        params.upsampleReg = 20;
        params.nIter = 500;
        params.checkEvery = 25;
        params.stopTol = 1e-7;
        params.minIter = 50;
        params.pixelSize = app.ptuOut.head.ImgHdr_PixResol * 1e3;
        params.lambda = 690;
        params.NA = 1.45;
        params.showPlots = false;
        params.useGPU = app.chkUseGPU.Value && gpuIsAvailable();

        try
            app.ismRes = run_ism_reconstruction_from_ptu(app.ptuOut, params);
            addStatus('ISM reconstruction done. Running FLIM...');
            drawnow;
            flim = reassigned_flim(app.ptuOut, app.ismRes, struct());
            app.flim = flim;
            addStatus('ISM-FLIM done.');
            showTauMean();
        catch ME
            addStatus(['ISM-FLIM failed: ' ME.message]);
        end
    end

    function onShowIntensity(~, ~)
        showIntensityFromPTU();
    end

    function onShowTauMean(~, ~)
        showTauMean();
    end

    function onSelectROI(~, ~)
        if isempty(app.axImage) || isempty(app.axImage.Children)
            addStatus('Display an image first.');
            return;
        end
        if ~isempty(app.roi) && isvalid(app.roi)
            delete(app.roi);
        end
        app.roi = drawrectangle(app.axImage, 'Color', 'y');
        addStatus('ROI selected.');
    end

    function onShowTCSPC(~, ~)
        if isempty(app.ptuOut) || ~isfield(app.ptuOut, 'tcspc_pix') || isempty(app.ptuOut.tcspc_pix)
            addStatus('No TCSPC cube available.');
            return;
        end
        if isempty(app.roi) || ~isvalid(app.roi)
            addStatus('Select an ROI first.');
            return;
        end
        [counts, tAxisNs] = tcspcFromROI(app.ptuOut, app.roi);
        nBins = app.txtBins.Value;
        [tBinNs, cBin] = logBinTCSPC(counts, tAxisNs, nBins);
        cla(app.axTCSPC);
        semilogx(app.axTCSPC, tBinNs*1e-9, cBin, 'LineWidth', 1.2);
        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Time (s)');
        ylabel(app.axTCSPC, 'Counts');
        title(app.axTCSPC, 'ROI TCSPC (log bins)');
    end

    function onSaveMAT(~, ~)
        if isempty(app.ptuOut)
            addStatus('Nothing to save.');
            return;
        end
        [f, p] = uiputfile('*.mat', 'Save MAT');
        if isequal(f, 0)
            return;
        end
        outFile = fullfile(p, f);
        ptuOut = app.ptuOut; %#ok<NASGU>
        flim = app.flim; %#ok<NASGU>
        ismRes = app.ismRes; %#ok<NASGU>
        save(outFile, 'ptuOut', 'flim', 'ismRes');
        addStatus(['Saved: ' outFile]);
    end

    function onSavePNG(~, ~)
        if isempty(app.axImage) || isempty(app.axImage.Children)
            addStatus('No image to save.');
            return;
        end
        [f, p] = uiputfile('*.png', 'Save PNG');
        if isequal(f, 0)
            return;
        end
        outFile = fullfile(p, f);
        exportgraphics(app.axImage, outFile, 'Resolution', 300);
        addStatus(['Saved: ' outFile]);
    end

    function onPointFCS(~, ~)
        if isempty(app.ptuOut)
            app.statusFCS.Text = 'Load a PTU first.';
            return;
        end
        app.statusFCS.Text = 'Running point FCS...';
        drawnow;

        head = app.ptuOut.head;
        Timeunit = head.MeasDesc_GlobalResolution;
        maxtime = app.editMaxtime.Value;
        Nsub = 10;
        Ncasc = ceil(log2(maxtime/Timeunit/Nsub));

        % Photon selection: ROI if available, otherwise all photons
        idx = true(numel(app.ptuOut.im_sync), 1);
        if ~isempty(app.roi) && isvalid(app.roi)
            idx = roiPhotonMask(app.ptuOut, app.roi);
        end

        tPhoton = double(app.ptuOut.im_sync(idx));
        num = true(numel(tPhoton), 1);

        [auto, autotime] = tttr2xfcsSym(tPhoton, num, Ncasc, Nsub);
        autotime = autotime * Timeunit;

        app.fcs.point.auto = auto;
        app.fcs.point.autotime = autotime;

        cla(app.axFCS);
        semilogx(app.axFCS, autotime, squeeze(auto(:,1,1)));
        grid(app.axFCS, 'on');
        xlabel(app.axFCS, 'Lag time (s)');
        ylabel(app.axFCS, 'G(\tau)');
        title(app.axFCS, 'Point FCS');

        app.statusFCS.Text = 'Point FCS done.';
    end

    function onLineFCS(~, ~)
        if isempty(app.lastFile)
            app.statusFCS.Text = 'Load a PTU first.';
            return;
        end
        app.statusFCS.Text = 'Running line-scan FCS...';
        drawnow;

        cnum = app.editCnum.Value;
        maxtime = app.editMaxtime.Value;
        timegates = [];
        try
            timegates = eval(app.editTimegates.Value);
        catch
            timegates = [];
        end

        try
            res = lsFCS(app.lastFile, cnum, maxtime, timegates);
            app.fcs.linescan = res;
            [G, ~, ~, ~, t, xxi] = lsCrossRead(res, 25);
            cla(app.axFCS);
            imagesc(app.axFCS, xxi, t, G(:,:,1));
            axis(app.axFCS, 'xy');
            set(app.axFCS, 'YScale', 'log');
            colorbar(app.axFCS);
            xlabel(app.axFCS, '\xi');
            ylabel(app.axFCS, '\tau (s)');
            title(app.axFCS, 'Line-scan FCS (component 1)');
            app.statusFCS.Text = 'Line-scan FCS done.';
        catch ME
            app.statusFCS.Text = ['Line-scan FCS failed: ' ME.message];
        end
    end

    function onSaveFCS(~, ~)
        if isempty(app.fcs)
            app.statusFCS.Text = 'No FCS data to save.';
            return;
        end
        [f, p] = uiputfile('*.mat', 'Save FCS MAT');
        if isequal(f, 0)
            return;
        end
        outFile = fullfile(p, f);
        fcs = app.fcs; %#ok<NASGU>
        save(outFile, 'fcs');
        app.statusFCS.Text = ['Saved: ' outFile];
    end

    % ----------------------- Helpers -----------------------
    function addStatus(msg)
        app.status.Text = msg;
    end

    function showIntensityFromPTU()
        if isempty(app.ptuOut)
            return;
        end
        if isfield(app.ptuOut, 'tags') && ~isempty(app.ptuOut.tags)
            img = sum(app.ptuOut.tags, 3);
        elseif isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix)
            img = sum(app.ptuOut.tcspc_pix, 3);
            img = sum(img, 4);
        else
            addStatus('No intensity data available.');
            return;
        end
        imagesc(app.axImage, img);
        axis(app.axImage, 'image');
        colorbar(app.axImage);
        title(app.axImage, 'Intensity');
    end

    function showTauMean()
        if isempty(app.flim)
            addStatus('No FLIM result available.');
            return;
        end
        if isfield(app.flim, 'total') && isfield(app.flim.total, 'tauMean')
            img = app.flim.total.tauMean;
        elseif isfield(app.flim, 'reassigned') && isfield(app.flim.reassigned, 'total')
            img = app.flim.reassigned.total.tauMean;
        else
            addStatus('No tauMean in FLIM result.');
            return;
        end
        imagesc(app.axImage, img);
        axis(app.axImage, 'image');
        colorbar(app.axImage);
        title(app.axImage, 'Tau mean (ns)');
    end
end

% ----------------------- Local functions -----------------------
function ok = gpuIsAvailable()
    ok = false;
    try
        ok = gpuDeviceCount > 0;
    catch
        ok = false;
    end
end

function flim = quickFLIMFromTCSPC(ptuOut, useGPU)
    cube = ptuOut.tcspc_pix;
    dt_ns = ptuOut.head.MeasDesc_Resolution * 1e9;
    Ng = size(cube, 3);
    tAxisNs = ((1:Ng) - 0.5) * dt_ns;

    if useGPU
        cubeG = gpuArray(double(cube));
        globalDecay = gather(squeeze(sum(sum(cubeG,1),2)));
    else
        cubeG = double(cube);
        globalDecay = squeeze(sum(sum(cubeG,1),2));
    end

    [~, t0Bin] = max(globalDecay);
    tShift = tAxisNs - tAxisNs(t0Bin);
    tShift3 = reshape(tShift, 1, 1, []);

    if useGPU
        tag = gather(sum(cubeG, 3));
        meanArrival = gather(sum(cubeG .* tShift3, 3) ./ max(tag, 1));
    else
        tag = sum(cubeG, 3);
        meanArrival = sum(cubeG .* tShift3, 3) ./ max(tag, 1);
    end

    flim = struct();
    flim.total = struct();
    flim.total.tag = tag;
    flim.total.tauMean = meanArrival;
    flim.total.tauRMS = zeros(size(meanArrival));
    flim.total.meanArrival = meanArrival;
    flim.total.globalDecay = globalDecay;
    flim.total.tAxisNs = tAxisNs(:);
    flim.total.t0Bin = t0Bin;
end

function [counts, tAxisNs] = tcspcFromROI(ptuOut, roi)
    cube = ptuOut.tcspc_pix;
    dt_ns = ptuOut.head.MeasDesc_Resolution * 1e9;
    Ng = size(cube, 3);
    tAxisNs = ((1:Ng) - 0.5) * dt_ns;

    pos = roi.Position;
    x1 = max(1, floor(pos(1)));
    y1 = max(1, floor(pos(2)));
    x2 = min(size(cube, 2), ceil(pos(1) + pos(3)));
    y2 = min(size(cube, 1), ceil(pos(2) + pos(4)));

    sub = cube(y1:y2, x1:x2, :, :);
    counts = squeeze(sum(sum(sum(sub,4),2),1));
end

function [tBinNs, cBin] = logBinTCSPC(counts, tAxisNs, nBins)
    t = tAxisNs(:);
    c = counts(:);
    tpos = t(t > 0);
    tmin = min(tpos);
    tmax = max(tpos);
    edges = logspace(log10(tmin), log10(tmax), nBins + 1);
    [~,~,bin] = histcounts(t, edges);
    cBin = accumarray(bin(bin>0), c(bin>0), [nBins 1], @sum, 0);
    tBinNs = sqrt(edges(1:end-1) .* edges(2:end));
end

function idx = roiPhotonMask(ptuOut, roi)
    pos = roi.Position;
    x1 = max(1, floor(pos(1)));
    y1 = max(1, floor(pos(2)));
    x2 = ceil(pos(1) + pos(3));
    y2 = ceil(pos(2) + pos(4));

    col = double(ptuOut.im_col);
    row = double(ptuOut.im_line);
    idx = col >= x1 & col <= x2 & row >= y1 & row <= y2;
end
