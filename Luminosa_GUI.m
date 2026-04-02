function Luminosa_GUI
% Luminosa_GUI
% Simple GUI for FLIM / ISM-FLIM processing and FCS analysis.

    app = struct();
    app.ptuOut = [];
    app.ismRes = [];
    app.flim = [];
    app.pattern = [];
    app.fcs = [];
    app.roi = [];
    app.lastFile = '';
    app.tcspc = [];
    app.tcspcFit = [];
    app.gamma = 1;
    app.displayMode = 'intensity';
    app.cbTau = [];

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
        g.ColumnWidth = {380, '1x'};
        g.RowHeight = {'1x'};

        ctl = uipanel(g, 'Title', 'Controls');
        try
            ctl.Scrollable = 'on';
        catch
        end
        ctlInner = uipanel(ctl, 'BorderType', 'none');
        axp = uipanel(g, 'Title', 'Display');

        rowHeights = [40,34,44,44,44,34,34,34,40,34,34,34,34,120];
        ctlGrid = uigridlayout(ctlInner, [14 1]);
        ctlGrid.RowHeight = num2cell(rowHeights);
        ctlGrid.RowSpacing = 8;
        ctlGrid.Padding = [6 6 6 6];

        innerH = sum(rowHeights) + ctlGrid.RowSpacing * (numel(rowHeights)-1) + ctlGrid.Padding(2) + ctlGrid.Padding(4) + 20;
        ctlInner.Units = 'pixels';
        ctlInner.Position = [0 0 360 innerH];

        app.btnLoad = uibutton(ctlGrid, 'Text', 'Load PTU', 'ButtonPushedFcn', @onLoadPTU, 'FontSize', 12);

        rowChecks = uigridlayout(ctlGrid, [1 2]);
        rowChecks.ColumnWidth = {'1x','1x'};
        app.chkStoreTCSPC = uicheckbox(rowChecks, 'Text', 'Store TCSPC', 'Value', true);
        app.chkUseGPU = uicheckbox(rowChecks, 'Text', 'Use GPU', 'Value', false);

        rowGamma = uigridlayout(ctlGrid, [1 3]);
        rowGamma.ColumnWidth = {70,'1x',60};
        uilabel(rowGamma, 'Text', 'Gamma');
        app.sliderGamma = uislider(rowGamma, 'Limits', [0.3 1.5], 'Value', 1, ...
            'MajorTicks', [0.3 1 1.5], 'ValueChangedFcn', @onGammaChanged);
        app.lblGamma = uilabel(rowGamma, 'Text', '1.00');

        rowRun = uigridlayout(ctlGrid, [1 4]);
        rowRun.ColumnWidth = {'1x','1x','1x','1x'};
        app.btnQuickFLIM = uibutton(rowRun, 'Text', 'FLIM', 'ButtonPushedFcn', @onQuickFLIM, 'FontSize', 11);
        app.btnISMFLIM = uibutton(rowRun, 'Text', 'ISM', 'ButtonPushedFcn', @onISMFLIM, 'FontSize', 11);
        app.btnShowIntensity = uibutton(rowRun, 'Text', 'Intensity', 'ButtonPushedFcn', @onShowIntensity, 'FontSize', 11);
        app.btnShowTau = uibutton(rowRun, 'Text', 'Tau', 'ButtonPushedFcn', @onShowTauMean, 'FontSize', 11);

        rowROI = uigridlayout(ctlGrid, [1 3]);
        rowROI.ColumnWidth = {'1x','1x','1x'};
        app.btnSelectROI = uibutton(rowROI, 'Text', 'ROI', 'ButtonPushedFcn', @onSelectROI, 'FontSize', 11);
        app.btnShowTCSPC = uibutton(rowROI, 'Text', 'Show TCSPC', 'ButtonPushedFcn', @onShowTCSPC, 'FontSize', 11);
        app.btnFitTCSPC = uibutton(rowROI, 'Text', 'Fit TCSPC', 'ButtonPushedFcn', @onFitTCSPC, 'FontSize', 11);

        rowNsub = uigridlayout(ctlGrid, [1 2]);
        rowNsub.ColumnWidth = {110,'1x'};
        uilabel(rowNsub, 'Text', 'TCSPC Nsub');
        app.txtBins = uieditfield(rowNsub, 'numeric', 'Value', 10, 'Limits', [2 32], ...
            'LowerLimitInclusive','on','UpperLimitInclusive','on');

        rowTau0 = uigridlayout(ctlGrid, [1 2]);
        rowTau0.ColumnWidth = {80,'1x'};
        uilabel(rowTau0, 'Text', 'Tau0 (ns)');
        app.editTau0 = uieditfield(rowTau0, 'text', 'Value', '0.35 1.5 5');

        rowFit = uigridlayout(ctlGrid, [1 2]);
        rowFit.ColumnWidth = {'1x','1x'};
        app.chkOptimizeTau = uicheckbox(rowFit, 'Text', 'Optimize tau', 'Value', true);
        app.chkIncludeBG = uicheckbox(rowFit, 'Text', 'Include BG', 'Value', true);

        rowRange = uigridlayout(ctlGrid, [1 5]);
        rowRange.ColumnWidth = {80,70,'1x','1x',70};
        uilabel(rowRange, 'Text', 'Tau range');
        app.chkAutoTauRange = uicheckbox(rowRange, 'Text', 'Auto', 'Value', true);
        app.editTauMin = uieditfield(rowRange, 'numeric', 'Value', 1.0, 'Limits', [0 1e3]);
        app.editTauMax = uieditfield(rowRange, 'numeric', 'Value', 7.0, 'Limits', [0 1e3]);
        app.editTauMin.Tooltip = 'Min tau (ns)';
        app.editTauMax.Tooltip = 'Max tau (ns)';
        app.btnApplyTauRange = uibutton(rowRange, 'Text', 'Apply', 'ButtonPushedFcn', @onApplyTauRange);

        rowPattern = uigridlayout(ctlGrid, [1 2]);
        rowPattern.ColumnWidth = {110,'1x'};
        uilabel(rowPattern, 'Text', 'Pattern source');
        app.dropPatternSource = uidropdown(rowPattern, 'Items', {'Confocal (ptuOut)','ISM reassigned (flim)'}, ...
            'Value', 'Confocal (ptuOut)');

        rowIRF = uigridlayout(ctlGrid, [1 2]);
        rowIRF.ColumnWidth = {80,'1x'};
        uilabel(rowIRF, 'Text', 'IRF');
        app.dropIRF = uidropdown(rowIRF, 'Items', {'GammaShifted','ExGauss','Simple'}, ...
            'Value', 'GammaShifted');

        app.btnPatternMatch = uibutton(ctlGrid, 'Text', 'Global Fit + Pattern Match', 'ButtonPushedFcn', @onPatternMatch, 'FontSize', 12);

        rowSave = uigridlayout(ctlGrid, [1 2]);
        rowSave.ColumnWidth = {'1x','1x'};
        app.btnSaveMAT = uibutton(rowSave, 'Text', 'Save MAT', 'ButtonPushedFcn', @onSaveMAT, 'FontSize', 11);
        app.btnSavePNG = uibutton(rowSave, 'Text', 'Save PNG', 'ButtonPushedFcn', @onSavePNG, 'FontSize', 11);

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
        g.ColumnWidth = {380, '1x'};
        g.RowHeight = {'1x'};

        ctl = uipanel(g, 'Title', 'Controls');
        try
            ctl.Scrollable = 'on';
        catch
        end
        ctlInner = uipanel(ctl, 'BorderType', 'none');
        axp = uipanel(g, 'Title', 'Display');

        rowHeights = [30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,100];
        ctlGrid = uigridlayout(ctlInner, [18 1]);
        ctlGrid.RowHeight = num2cell(rowHeights);
        ctlGrid.RowSpacing = 6;
        ctlGrid.Padding = [6 6 6 6];
        innerH = sum(rowHeights) + ctlGrid.RowSpacing * (numel(rowHeights)-1) + ctlGrid.Padding(2) + ctlGrid.Padding(4) + 20;
        ctlInner.Units = 'pixels';
        ctlInner.Position = [0 0 360 innerH];

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
        useGPU = app.chkUseGPU.Value && gpuIsAvailable() && storeTCSPC;
        try
            app.flim = [];
            app.ismRes = [];
            app.pattern = [];
            app.tcspcFit = [];
            app.ptuOut = PTU_MultiFrameScanReadFast(app.lastFile, 1e6, storeTCSPC, useGPU);
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

    function onGammaChanged(~, ~)
        app.gamma = app.sliderGamma.Value;
        app.lblGamma.Text = sprintf('%.2f', app.gamma);
        switch app.displayMode
            case 'tau'
                showTauMean();
            case 'pattern'
                showPatternOverlay();
            otherwise
                showIntensityFromPTU();
        end
    end

    function onApplyTauRange(~, ~)
        if strcmp(app.displayMode, 'tau')
            showTauMean();
        elseif strcmp(app.displayMode, 'pattern')
            showPatternOverlay();
        end
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
        [counts, dtNs] = tcspcFromROI(app.ptuOut, app.roi);
        nSub = app.txtBins.Value;
        [tBinNs, cBin, edgesNs] = multiTauBinTCSPC(counts, dtNs, nSub);
        app.tcspc = struct('tBinNs', tBinNs, 'cBin', cBin, 'edgesNs', edgesNs, 'dtNs', dtNs);
        cla(app.axTCSPC);
        semilogx(app.axTCSPC, tBinNs*1e-9, cBin, 'LineWidth', 1.2);
        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Time (s)');
        ylabel(app.axTCSPC, 'Counts');
        title(app.axTCSPC, 'ROI TCSPC (multitau)');
    end

    function onFitTCSPC(~, ~)
        if isempty(app.tcspc)
            onShowTCSPC();
            if isempty(app.tcspc)
                return;
            end
        end

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5 5".');
            return;
        end

        includeBG = app.chkIncludeBG.Value;
        optimizeTau = app.chkOptimizeTau.Value;

        [tauFit, coeff, fitCounts] = fitRoiTCSPC(app.tcspc.edgesNs, app.tcspc.cBin, tau0, includeBG, optimizeTau);
        app.tcspcFit = struct('tauFit', tauFit, 'coeff', coeff, 'fitCounts', fitCounts);

        cla(app.axTCSPC);
        semilogx(app.axTCSPC, app.tcspc.tBinNs*1e-9, app.tcspc.cBin, 'k.', 'MarkerSize', 10); hold(app.axTCSPC, 'on');
        semilogx(app.axTCSPC, app.tcspc.tBinNs*1e-9, fitCounts, 'r-', 'LineWidth', 1.2);
        hold(app.axTCSPC, 'off');
        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Time (s)');
        ylabel(app.axTCSPC, 'Counts');
        title(app.axTCSPC, 'ROI TCSPC fit (multitau)');

        addStatus(sprintf('ROI fit tau (ns): %s', num2str(tauFit(:).', '%.3g ')));
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
        pattern = app.pattern; %#ok<NASGU>
        tcspcFit = app.tcspcFit; %#ok<NASGU>
        save(outFile, 'ptuOut', 'flim', 'ismRes', 'pattern', 'tcspcFit');
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

    function onPatternMatch(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end
        if isempty(app.ptuOut.tcspc_pix)
            addStatus('No tcspc_pix. Reload with "Store TCSPC cube".');
            return;
        end

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5 5".');
            return;
        end

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        includeBG = app.chkIncludeBG.Value;

        % Select tcspc source
        src = app.dropPatternSource.Value;
        if contains(src, 'ISM') && ~isempty(app.flim) && isfield(app.flim, 'reassigned') ...
                && isfield(app.flim.reassigned, 'total') && isfield(app.flim.reassigned.total, 'xyT') ...
                && ~isempty(app.flim.reassigned.total.xyT)
            tcspc_pix = app.flim.reassigned.total.xyT;
        else
            tcspc_pix = app.ptuOut.tcspc_pix;
        end

        head = app.ptuOut.head;
        dtNs = head.MeasDesc_Resolution * 1e9;
        pulsePeriodNs = head.MeasDesc_GlobalResolution * 1e9;

        addStatus('Estimating IRF...');
        drawnow;
        irf = [];
        try
            switch app.dropIRF.Value
                case 'GammaShifted'
                    outIRF = Calc_mIRF_Global_GammaShifted(head, tcspc_pix, tau0);
                    irf = outIRF.IRF(:);
                case 'ExGauss'
                    outIRF = Calc_mIRF_Global_ExGauss(head, tcspc_pix, tau0);
                    irf = outIRF.IRF(:);
                otherwise
                    globalDecay = squeeze(sum(sum(tcspc_pix,1),2));
                    head2 = head;
                    head2.Resolution = dtNs;
                    if isfield(head, 'TTResult_SyncRate')
                        head2.SyncRate = head.TTResult_SyncRate;
                    elseif isfield(head, 'SyncRate')
                        head2.SyncRate = head.SyncRate;
                    end
                    irf = Calc_mIRF(head2, globalDecay');
            end
        catch ME
            addStatus(['IRF failed: ' ME.message]);
            return;
        end

        if isempty(irf)
            addStatus('IRF is empty.');
            return;
        end

        addStatus('Running global fit + pattern match...');
        drawnow;

        opts = struct();
        opts.useGPU = useGPU;
        opts.mode = 'PIRLS';
        opts.batchSize = 4000;
        opts.includeBackground = includeBG;
        opts.pieIndex = 1;
        opts.normalizePatterns = true;
        opts.sortLifetimes = true;

        try
            outFLIM = GlobalMultiExpPatternMatchFromTCSPC( ...
                tcspc_pix, irf, pulsePeriodNs, dtNs, tau0, opts);
            app.pattern = outFLIM;
        catch ME
            addStatus(['Pattern match failed: ' ME.message]);
            return;
        end

        % Show amplitude maps
        try
            cmaps = {'gray','red','green','blue','orange','cyan','magenta'};
            nMaps = size(outFLIM.Amp, 3);
            cmaps = cmaps(1:min(nMaps, numel(cmaps)));
            ShowImagesMultiCmapWithColorbars(outFLIM.Amp, cmaps, [], [], true);
        catch
        end

        showPatternOverlay();
        addStatus('Pattern match complete.');
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
        img = getIntensityMap();
        if isempty(img)
            addStatus('No intensity data available.');
            return;
        end
        imgG = applyGamma(img, app.gamma);
        imagesc(app.axImage, imgG);
        axis(app.axImage, 'image');
        colormap(app.axImage, gray);
        colorbar(app.axImage);
        title(app.axImage, 'Intensity (gamma)');
        if ~isempty(app.cbTau) && isvalid(app.cbTau)
            delete(app.cbTau);
            app.cbTau = [];
        end
        app.displayMode = 'intensity';
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
        intensity = getIntensityMap();
        if isempty(intensity)
            imagesc(app.axImage, img);
            axis(app.axImage, 'image');
            colorbar(app.axImage);
            title(app.axImage, 'Tau mean (ns)');
        else
            showTauOverlay(img, intensity, 'Tau mean (ns)');
        end
        app.displayMode = 'tau';
    end

    function showPatternOverlay()
        if isempty(app.pattern)
            addStatus('No pattern-match result available.');
            return;
        end
        tau = app.pattern.tauMeanArithmetic;
        intensity = app.pattern.intensity;
        showTauOverlay(tau, intensity, 'Pattern-match tau (ns)');
        app.displayMode = 'pattern';
    end

    function img = getIntensityMap()
        img = [];
        if ~isempty(app.flim) && isfield(app.flim, 'total') && isfield(app.flim.total, 'tag')
            img = app.flim.total.tag;
            return;
        end
        if ~isempty(app.ptuOut)
            if isfield(app.ptuOut, 'tags') && ~isempty(app.ptuOut.tags)
                img = sum(app.ptuOut.tags, 3);
                return;
            elseif isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix)
                tmp = sum(app.ptuOut.tcspc_pix, 3);
                img = sum(tmp, 4);
                return;
            end
        end
    end

    function imgG = applyGamma(img, gamma)
        img = double(img);
        if isempty(img)
            imgG = img;
            return;
        end
        img = img - min(img(:));
        if max(img(:)) > 0
            img = img ./ max(img(:));
        end
        g = max(0.1, min(3, gamma));
        imgG = img .^ g;
    end

    function showTauOverlay(tau, intensity, titleStr)
        cmap = jet(256);
        cmap = cmap(30:end-30,:);

        tau = double(tau);
        intensity = double(intensity);
        if isempty(tau) || isempty(intensity)
            return;
        end
        trange = getTauRange(tau);
        bright = applyGamma(intensity, app.gamma);
        rgb = flim_rgb(tau, bright, trange, cmap);
        image(app.axImage, rgb);
        axis(app.axImage, 'image');
        axis(app.axImage, 'off');
        title(app.axImage, titleStr);
        colormap(app.axImage, cmap);
        caxis(app.axImage, trange);
        if ~isempty(app.cbTau) && isvalid(app.cbTau)
            delete(app.cbTau);
        end
        app.cbTau = colorbar(app.axImage);
        ylabel(app.cbTau, 'Tau (ns)');
    end

    function trange = getTauRange(tau)
        tvals = tau(isfinite(tau));
        if isempty(tvals)
            trange = [1 5];
            return;
        end
        if app.chkAutoTauRange.Value
            trange = prctile(tvals, [5 95]);
        else
            trange = [app.editTauMin.Value, app.editTauMax.Value];
        end
        if trange(1) == trange(2)
            trange = [min(tvals) max(tvals)];
        end
        if trange(1) < 0
            trange(1) = 0;
        end
        if trange(2) <= trange(1)
            trange(2) = trange(1) + eps;
        end
    end

    function rgb = flim_rgb(tau, inten, tauRange, cmap)
        tmin = tauRange(1);
        tmax = tauRange(2);
        if tmax <= tmin
            tmax = tmin + eps;
        end
        tau = min(max(tau, tmin), tmax);
        tnorm = (tau - tmin) / (tmax - tmin);
        idx = 1 + tnorm * (size(cmap,1) - 1);

        rgb = zeros([size(tau) 3]);
        for c = 1:3
            rgb(:,:,c) = interp1(1:size(cmap,1), cmap(:,c), idx, 'linear', 'extrap') .* inten;
        end
        rgb = max(min(rgb, 1), 0);
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

function [counts, dtNs] = tcspcFromROI(ptuOut, roi)
    cube = ptuOut.tcspc_pix;
    dtNs = ptuOut.head.MeasDesc_Resolution * 1e9;

    pos = roi.Position;
    x1 = max(1, floor(pos(1)));
    y1 = max(1, floor(pos(2)));
    x2 = min(size(cube, 2), ceil(pos(1) + pos(3)));
    y2 = min(size(cube, 1), ceil(pos(2) + pos(4)));

    sub = cube(y1:y2, x1:x2, :, :);
    counts = squeeze(sum(sum(sum(sub,4),2),1));
end

function [tBinNs, cBin, edgesNs] = multiTauBinTCSPC(counts, dtNs, nSub)
% Multi-tau logscale binning (FCS-style).
    if nargin < 3 || isempty(nSub)
        nSub = 10;
    end
    nSub = max(2, round(nSub));

    counts = double(counts(:));
    Nt = numel(counts);
    if Nt < 2
        edgesNs = [0 dtNs];
        tBinNs = dtNs/2;
        cBin = counts;
        return;
    end

    Ncasc = max(1, ceil(log2(Nt / nSub)));
    widths = repelem(2.^(0:Ncasc-1), nSub);
    edgesIdx = [0 cumsum(widths)];
    edgesIdx(edgesIdx > Nt) = Nt;
    edgesIdx = unique(edgesIdx);
    if edgesIdx(end) < Nt
        edgesIdx = [edgesIdx Nt];
    end
    if numel(edgesIdx) < 2
        edgesIdx = [0 Nt];
    end

    nBins = numel(edgesIdx) - 1;
    csum = cumsum([0; counts]);
    cBin = zeros(nBins, 1);
    for i = 1:nBins
        a = edgesIdx(i) + 1;
        b = edgesIdx(i+1);
        if b < a
            cBin(i) = 0;
            continue;
        end
        cBin(i) = csum(b+1) - csum(a);
    end

    edgesNs = edgesIdx * dtNs;
    tStart = edgesNs(1:end-1);
    tEnd = edgesNs(2:end);
    tBinNs = sqrt(max(tStart, eps) .* max(tEnd, eps));
    if tStart(1) == 0
        tBinNs(1) = tEnd(1) / 2;
    end
end

function tau0 = parseTau0(str)
    tau0 = [];
    if isempty(str)
        return;
    end
    if isnumeric(str)
        tau0 = str(:).';
        return;
    end
    try
        tau0 = str2num(str); %#ok<ST2NM>
    catch
        tau0 = [];
    end
    if isempty(tau0) || any(~isfinite(tau0))
        tau0 = [];
    end
end

function [tauFit, coeff, fitCounts] = fitRoiTCSPC(edgesNs, counts, tau0, includeBG, optimizeTau)
    tau0 = double(tau0(:));
    counts = double(counts(:));

    if optimizeTau
        p0 = tau0;
        xmin = max(0.05, p0/10);
        xmax = p0 * 10;
        tol = 1e-5;
        steps = 200 * numel(p0);
        [pfit, ~] = Simplex(@roiTcspcErr, p0, xmin, xmax, tol, steps, edgesNs, counts, includeBG);
    else
        pfit = tau0;
    end

    [~, coeff, fitCounts] = roiTcspcErr(pfit, edgesNs, counts, includeBG);
    tauFit = pfit(:).';
end

function [err, coeff, fitCounts] = roiTcspcErr(tau, edgesNs, counts, includeBG)
    tau = max(double(tau(:)), 1e-6);
    t0 = edgesNs(1:end-1);
    t1 = edgesNs(2:end);
    nExp = numel(tau);
    nBins = numel(t0);

    M = zeros(nBins, nExp + includeBG);
    for k = 1:nExp
        tk = tau(k);
        M(:, includeBG + k) = tk * (exp(-t0./tk) - exp(-t1./tk));
    end
    if includeBG
        M(:,1) = (t1 - t0);
    end

    coeff = lsqnonneg(M, counts);
    fitCounts = M * coeff;
    err = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
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
