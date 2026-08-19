function Luminosa_GUI_globalIRF_autofit_full
% Luminosa_GUI_globalIRF_autofit_full
% Standalone GUI for FLIM / ISM-FLIM / FCS analysis.
% GUI Tau0 is used only for TCSPC/pattern fitting.
% IRF is always estimated from the whole-file TCSPC using internal default lifetime guesses.
% Main changes:
% - left controls resized so buttons are fully visible
% - all action-button rows have uniform height
% - ROI TCSPC is shown as log-binned COUNT DENSITY, not raw counts per variable-width bin
% - IRF estimation caches results and reduces the TCSPC cube to a 1D global decay first
% - GPU is used where it helps most here: cube reduction, not tiny-parameter optimization

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
    app.tcspcGlobal = [];
    app.irfGlobal = [];
    app.irfGlobalModel = '';
    app.irfGlobalMeta = struct();
    app.irfGlobalDtNs = [];
    app.globalFit = [];
    app.tcspcDisplayMode = 'none';
    app.gamma = 1;
    app.displayMode = 'intensity';
    app.cbTau = [];
    app.irfCache = struct('key', '', 'irf', [], 'meta', struct());
    app.txtFitResults = [];
    app.txtLegend = [];
    app.pnlTCSPCSide = [];
    app.channelGallery = [];

    app.fig = uifigure('Name','Luminosa FLIM / ISM-FLIM / FCS', 'Position',[80 60 1360 840]);
    app.grid = uigridlayout(app.fig, [1 1]);
    app.grid.Padding = [4 4 4 4];
    app.tg = uitabgroup(app.grid);

    tabFLIM = uitab(app.tg, 'Title', 'FLIM / ISM-FLIM');
    tabFCS  = uitab(app.tg, 'Title', 'FCS');

    buildFLIMTab(tabFLIM);
    buildFCSTab(tabFCS);

    function buildFLIMTab(parent)
        g = uigridlayout(parent, [1 2]);
        g.ColumnWidth = {440, '1x'};
        g.RowHeight = {'1x'};
        g.Padding = [4 4 4 4];
        g.ColumnSpacing = 6;

        ctl = uipanel(g, 'Title', 'Controls');
        axp = uipanel(g, 'Title', 'Display');

        try
            ctl.Scrollable = 'on';
        catch
        end

        ctlInner = uipanel(ctl, 'BorderType', 'none');
        ctlInner.Units = 'pixels';

        % Keep all action-button rows at one consistent height.
        rowHeights = [42, 28, 34, 36, 36, 30, 30, 30, 30, 30, 30, 30, 36, 36, 72];
        ctlGrid = uigridlayout(ctlInner, [numel(rowHeights) 1]);
        ctlGrid.RowHeight = num2cell(rowHeights);
        ctlGrid.RowSpacing = 6;
        ctlGrid.Padding = [8 8 8 8];

        innerH = sum(rowHeights) + ctlGrid.RowSpacing*(numel(rowHeights)-1) + ctlGrid.Padding(2) + ctlGrid.Padding(4) + 8;
        ctlInner.Position = [0 0 415 innerH];

        % Row 1
        app.btnLoad = uibutton(ctlGrid, 'Text', 'Load PTU', 'ButtonPushedFcn', @onLoadPTU, 'FontSize', 12);

        % Row 2
        rowChecks = uigridlayout(ctlGrid, [1 2]);
        rowChecks.ColumnWidth = {'1x','1x'};
        rowChecks.Padding = [0 0 0 0];
        app.chkStoreTCSPC = uicheckbox(rowChecks, 'Text', 'Store TCSPC', 'Value', true);
        app.chkUseGPU = uicheckbox(rowChecks, 'Text', 'Use GPU', 'Value', false);

        % Row 3
        rowGamma = uigridlayout(ctlGrid, [1 3]);
        rowGamma.ColumnWidth = {70,'1x',55};
        rowGamma.Padding = [0 0 0 0];
        uilabel(rowGamma, 'Text', 'Gamma');
        app.sliderGamma = uislider(rowGamma, 'Limits', [0.3 1.5], 'Value', 1, ...
            'MajorTicks', [0.3 1 1.5], 'ValueChangedFcn', @onGammaChanged);
        app.lblGamma = uilabel(rowGamma, 'Text', '1.00');

        % Row 4
        rowRun = uigridlayout(ctlGrid, [1 4]);
        rowRun.ColumnWidth = {'1x','1x','1x','1x'};
        rowRun.Padding = [0 0 0 0];
        app.btnQuickFLIM = uibutton(rowRun, 'Text', 'FLIM', 'ButtonPushedFcn', @onQuickFLIM, 'FontSize', 11);
        app.btnISMFLIM = uibutton(rowRun, 'Text', 'ISM', 'ButtonPushedFcn', @onISMFLIM, 'FontSize', 11);
        app.btnShowIntensity = uibutton(rowRun, 'Text', 'Intensity', 'ButtonPushedFcn', @onShowIntensity, 'FontSize', 11);
        app.btnShowTau = uibutton(rowRun, 'Text', 'Tau', 'ButtonPushedFcn', @onShowTauMean, 'FontSize', 11);

        % Row 5
        rowROI = uigridlayout(ctlGrid, [1 3]);
        rowROI.ColumnWidth = {'1x','1x','1x'};
        rowROI.Padding = [0 0 0 0];
        app.btnSelectROI = uibutton(rowROI, 'Text', 'ROI', 'ButtonPushedFcn', @onSelectROI, 'FontSize', 11);
        app.btnShowTCSPC = uibutton(rowROI, 'Text', 'Show TCSPC', 'ButtonPushedFcn', @onShowTCSPC, 'FontSize', 11);
        app.btnFitTCSPC = uibutton(rowROI, 'Text', 'Fit TCSPC', 'ButtonPushedFcn', @onFitTCSPC, 'FontSize', 11);

        % Row 6
        rowNsub = uigridlayout(ctlGrid, [1 2]);
        rowNsub.ColumnWidth = {120,'1x'};
        rowNsub.Padding = [0 0 0 0];
        uilabel(rowNsub, 'Text', 'Log bins / octave');
        app.txtBins = uieditfield(rowNsub, 'numeric', 'Value', 8, 'Limits', [2 32], ...
            'LowerLimitInclusive', 'on', 'UpperLimitInclusive', 'on');

        % Row 7
        rowTau0 = uigridlayout(ctlGrid, [1 2]);
        rowTau0.ColumnWidth = {85,'1x'};
        rowTau0.Padding = [0 0 0 0];
        uilabel(rowTau0, 'Text', 'Tau0 (ns)');
        app.editTau0 = uieditfield(rowTau0, 'text', 'Value', '0.35 1.5 5');

        % Row 8
        rowFit = uigridlayout(ctlGrid, [1 2]);
        rowFit.ColumnWidth = {'1x','1x'};
        rowFit.Padding = [0 0 0 0];
        app.chkOptimizeTau = uicheckbox(rowFit, 'Text', 'Optimize tau', 'Value', true);
        app.chkIncludeBG = uicheckbox(rowFit, 'Text', 'Include BG', 'Value', true);

        % Row 9
        rowRange = uigridlayout(ctlGrid, [1 5]);
        rowRange.ColumnWidth = {80,60,'1x','1x',70};
        rowRange.Padding = [0 0 0 0];
        uilabel(rowRange, 'Text', 'Tau range');
        app.chkAutoTauRange = uicheckbox(rowRange, 'Text', 'Auto', 'Value', true);
        app.editTauMin = uieditfield(rowRange, 'numeric', 'Value', 1.0, 'Limits', [0 1e3]);
        app.editTauMax = uieditfield(rowRange, 'numeric', 'Value', 7.0, 'Limits', [0 1e3]);
        app.editTauMin.Tooltip = 'Min tau (ns)';
        app.editTauMax.Tooltip = 'Max tau (ns)';
        app.btnApplyTauRange = uibutton(rowRange, 'Text', 'Apply', 'ButtonPushedFcn', @onApplyTauRange);

        % Row 10
        rowPattern = uigridlayout(ctlGrid, [1 2]);
        rowPattern.ColumnWidth = {110,'1x'};
        rowPattern.Padding = [0 0 0 0];
        uilabel(rowPattern, 'Text', 'Pattern source');
        app.dropPatternSource = uidropdown(rowPattern, 'Items', {'Confocal (ptuOut)','ISM reassigned (flim)'}, ...
            'Value', 'Confocal (ptuOut)');

        % Row 11
        rowIRF = uigridlayout(ctlGrid, [1 2]);
        rowIRF.ColumnWidth = {85,'1x'};
        rowIRF.Padding = [0 0 0 0];
        uilabel(rowIRF, 'Text', 'IRF model');
        app.dropIRF = uidropdown(rowIRF, 'Items', {'GammaShifted','ExGauss','Simple'}, ...
            'Value', 'GammaShifted', 'ValueChangedFcn', @onIRFModelChanged);

        % Row 12
        rowChunk = uigridlayout(ctlGrid, [1 2]);
        rowChunk.ColumnWidth = {120,'1x'};
        rowChunk.Padding = [0 0 0 0];
        uilabel(rowChunk, 'Text', 'Read chunk (M)');
        app.editChunkM = uieditfield(rowChunk, 'numeric', 'Value', 1, 'Limits', [0.1 50]);

        % Row 13
        app.btnPatternMatch = uibutton(ctlGrid, 'Text', 'Global Fit + Pattern Match', ...
            'ButtonPushedFcn', @onPatternMatch, 'FontSize', 12);

        % Row 14
        rowSave = uigridlayout(ctlGrid, [1 2]);
        rowSave.ColumnWidth = {'1x','1x'};
        rowSave.Padding = [0 0 0 0];
        app.btnSaveMAT = uibutton(rowSave, 'Text', 'Save MAT', 'ButtonPushedFcn', @onSaveMAT, 'FontSize', 11);
        app.btnSavePNG = uibutton(rowSave, 'Text', 'Save PNG', 'ButtonPushedFcn', @onSavePNG, 'FontSize', 11);

        % Row 15
        app.status = uilabel(ctlGrid, 'Text', 'Ready.');
        app.status.WordWrap = 'on';
        app.status.HorizontalAlignment = 'left';

        axGrid = uigridlayout(axp, [2 1]);
        axGrid.RowHeight = {'2.2x','1x'};
        axGrid.Padding = [6 6 6 6];

        app.axImage = uiaxes(axGrid);
        title(app.axImage, 'Image');

        tcspcGrid = uigridlayout(axGrid, [1 2]);
        tcspcGrid.ColumnWidth = {'3.8x','1.5x'};
        tcspcGrid.ColumnSpacing = 8;
        tcspcGrid.Padding = [0 0 0 0];

        app.axTCSPC = uiaxes(tcspcGrid);
        title(app.axTCSPC, 'TCSPC');
        app.axTCSPC.Box = 'on';

        app.pnlTCSPCSide = uipanel(tcspcGrid, 'Title', 'Fit results / legend');
        sideGrid = uigridlayout(app.pnlTCSPCSide, [4 1]);
        sideGrid.RowHeight = {18, '1x', 18, '0.8x'};
        sideGrid.Padding = [4 4 4 4];
        sideGrid.RowSpacing = 4;
        uilabel(sideGrid, 'Text', 'Fit results');
        app.txtFitResults = uitextarea(sideGrid, 'Editable', 'off', 'Value', {'No fit results yet.'});
        uilabel(sideGrid, 'Text', 'Legend');
        app.txtLegend = uitextarea(sideGrid, 'Editable', 'off', 'Value', {'No legend yet.'});

        if ~gpuIsAvailable()
            app.chkUseGPU.Enable = 'off';
            app.chkUseGPU.Value = false;
        end
    end

    function buildFCSTab(parent)
        g = uigridlayout(parent, [1 2]);
        g.ColumnWidth = {440, '1x'};
        g.RowHeight = {'1x'};
        g.Padding = [4 4 4 4];
        g.ColumnSpacing = 6;

        ctl = uipanel(g, 'Title', 'Controls');
        axp = uipanel(g, 'Title', 'Display');
        ctlInner = uipanel(ctl, 'BorderType', 'none');
        ctlInner.Units = 'pixels';
        try
            ctl.Scrollable = 'on';
        catch
        end

        rowHeights = [36,36,24,28,24,28,24,28,36,36,36,72];
        ctlGrid = uigridlayout(ctlInner, [numel(rowHeights) 1]);
        ctlGrid.RowHeight = num2cell(rowHeights);
        ctlGrid.RowSpacing = 6;
        ctlGrid.Padding = [8 8 8 8];
        innerH = sum(rowHeights) + ctlGrid.RowSpacing*(numel(rowHeights)-1) + ctlGrid.Padding(2) + ctlGrid.Padding(4) + 8;
        ctlInner.Position = [0 0 415 innerH];

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
        axGrid.Padding = [6 6 6 6];
        app.axFCS = uiaxes(axGrid);
        title(app.axFCS, 'FCS');
    end


    function onLoadPTU(~, ~)
        [f, p] = uigetfile('*.ptu', 'Select PTU file');
        if isequal(f, 0)
            return;
        end
        app.lastFile = fullfile(p, f);
        addStatus('Loading PTU...');
        drawnow;

        storeTCSPC = app.chkStoreTCSPC.Value;
        photonsPerChunk = max(1e5, app.editChunkM.Value * 1e6);
        useGPU = app.chkUseGPU.Value && gpuIsAvailable() && storeTCSPC;

        try
            app.flim = [];
            app.ismRes = [];
            app.pattern = [];
            app.tcspc = [];
            app.tcspcFit = [];
            app.tcspcGlobal = [];
            app.irfCache = struct('key', '', 'irf', [], 'meta', struct());
    app.txtFitResults = [];
    app.txtLegend = [];
    app.pnlTCSPCSide = [];
            app.irfGlobal = [];
            app.irfGlobalModel = '';
            app.irfGlobalMeta = struct();
            app.irfGlobalDtNs = [];
            app.globalFit = [];
            app.tcspcDisplayMode = 'none';
            clearTCSPCSidePanel();

            if ~isempty(app.roi) && isvalid(app.roi)
                delete(app.roi);
            end
            app.roi = [];

            if exist('PTU_MultiFrameScanReadFast_nativeTCSPC', 'file') == 2
                app.ptuOut = PTU_MultiFrameScanReadFast_nativeTCSPC(app.lastFile, photonsPerChunk, storeTCSPC, useGPU, true);
            elseif exist('PTU_MultiFrameScanReadFast_multiTauTCSPC', 'file') == 2
                app.ptuOut = PTU_MultiFrameScanReadFast_multiTauTCSPC(app.lastFile, photonsPerChunk, storeTCSPC, useGPU, true, 8, true);
            elseif exist('PTU_MultiFrameScanReadFast', 'file') == 2
                app.ptuOut = PTU_MultiFrameScanReadFast(app.lastFile, photonsPerChunk, storeTCSPC, useGPU);
            else
                error('No compatible PTU reader was found on the MATLAB path.');
            end

            img = getIntensityMap();
            if isempty(img)
                addStatus('PTU loaded, but no displayable intensity map was found in ptuOut.');
                return;
            end

            showIntensityFromPTU();
            irfOK = ensureGlobalIRF(true);
            fitOK = false;
            if irfOK
                fitOK = autoFitGlobalTCSPC();
            end
            showGlobalTCSPC();
            channelMsg = getChannelDiagnosticMessage();
            if irfOK && fitOK
                statusMsg = sprintf('PTU loaded. Whole-file IRF estimated with %s model; best global fit selected automatically.', app.dropIRF.Value);
            elseif irfOK
                statusMsg = sprintf('PTU loaded. Whole-file IRF estimated with %s model, but automatic global fit failed.', app.dropIRF.Value);
            else
                statusMsg = 'PTU loaded, but global IRF estimation failed.';
            end
            if ~isempty(channelMsg)
                statusMsg = sprintf('%s %s', statusMsg, channelMsg);
            end
            addStatus(statusMsg);
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
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        hasLinear = isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix);
        hasMultiTau = isfield(app.ptuOut, 'tcspc_pix_mt') && ~isempty(app.ptuOut.tcspc_pix_mt);
        if ~hasLinear && ~hasMultiTau
            addStatus('No TCSPC data available for quick FLIM.');
            return;
        end

        addStatus('Computing quick FLIM...');
        drawnow;

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        app.flim = quickFLIMFromTCSPCFlexible(app.ptuOut, useGPU);
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
            app.flim = reassigned_flim(app.ptuOut, app.ismRes, struct());
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

    function onIRFModelChanged(~, ~)
        if isempty(app.ptuOut)
            return;
        end
        ok = ensureGlobalIRF(true);
        fitOK = false;
        if ok
            fitOK = autoFitGlobalTCSPC();
        end
        if strcmp(app.tcspcDisplayMode, 'roi_fit') && ~isempty(app.tcspc)
            onFitTCSPC();
        elseif strcmp(app.tcspcDisplayMode, 'roi') && ~isempty(app.tcspc)
            onShowTCSPC();
        else
            showGlobalTCSPC();
        end
        if ok && fitOK
            addStatus(sprintf('Whole-file IRF recomputed with %s model; automatic global fit updated.', app.dropIRF.Value));
        elseif ok
            addStatus(sprintf('Whole-file IRF recomputed with %s model.', app.dropIRF.Value));
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
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end
        if isempty(app.roi) || ~isvalid(app.roi)
            addStatus('Select an ROI first.');
            return;
        end

        [countsNative, dtNsNative, srcInfo] = tcspcFromROI_native(app.ptuOut, app.roi);
        if isempty(countsNative)
            addStatus('No ROI TCSPC available. Reload with native photon TCSPC support or Store TCSPC.');
            return;
        end

        binsPerOct = app.txtBins.Value;
        [countsShift, shiftBins, shiftNs, peakIdx, displayStartIdx] = shiftDecayForIRFDisplay(countsNative, dtNsNative);
        [tBinNs, cDensity, cInt, edgesNs, widthNs] = logBinTCSPC(countsShift, dtNsNative, binsPerOct);

        app.tcspc = struct();
        app.tcspc.tBinNs = tBinNs;
        app.tcspc.cDensity = cDensity;
        app.tcspc.cBin = cInt;
        app.tcspc.edgesNs = edgesNs;
        app.tcspc.widthNs = widthNs;
        app.tcspc.dtNs = dtNsNative;
        app.tcspc.rawCountsNative = countsNative(:);
        app.tcspc.rawCountsShift = countsShift(:);
        app.tcspc.shiftBins = shiftBins;
        app.tcspc.shiftNs = shiftNs;
        app.tcspc.peakIdx = peakIdx;
        app.tcspc.displayStartIdx = displayStartIdx;
        app.tcspc.source = srcInfo;

        cla(app.axTCSPC);
        yyaxis(app.axTCSPC, 'left');
        semilogx(app.axTCSPC, tBinNs*1e-9, cDensity, 'k.', 'MarkerSize', 10); 
        hold(app.axTCSPC, 'on');

        legendLines = {'ROI TCSPC'};
        if ensureGlobalIRF(false) && ~isempty(app.irfGlobal) && (isempty(app.irfGlobalDtNs) || abs(app.irfGlobalDtNs - dtNsNative) < 1e-6)
            irfDisp = app.irfGlobal(:);
            if shiftBins < numel(irfDisp)
                irfDisp = irfDisp(shiftBins+1:end);
            end
            irfDisp = max(irfDisp, 0);
            if ~isempty(irfDisp) && max(irfDisp) > 0
                irfScale = max(countsShift) / max(irfDisp);
                irfScaled = irfDisp * irfScale;
                [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, dtNsNative, binsPerOct);
                semilogx(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
                legendLines{end+1} = 'Whole-file IRF (scaled)';
            end
        end
        hold(app.axTCSPC, 'off');
        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Delay time (s)');
        ylabel(app.axTCSPC, 'Counts / ns');
        app.axTCSPC.YScale = 'log';
        title(app.axTCSPC, sprintf('ROI TCSPC + whole-file IRF (%s)', app.dropIRF.Value));

        yyaxis(app.axTCSPC, 'right');
        app.axTCSPC.YAxis(2).Visible = 'off';
        app.axTCSPC.YAxis(2).Ticks = [];
        yyaxis(app.axTCSPC, 'left');

        setLegendEntries(legendLines);
        setFitSummary({sprintf('ROI source: %s', srcInfo.mode), ...
                       sprintf('Native dt = %.3f ps', 1e3*dtNsNative), ...
                       sprintf('Shift = %.3f ns', shiftNs), ...
                       'Press "Fit TCSPC" to run reconvolution fit.'});
        app.tcspcDisplayMode = 'roi';

        addStatus(sprintf('ROI TCSPC built from %s. GUI Tau0 now applies only to ROI fitting.', srcInfo.mode));
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

        if ~ensureGlobalIRF(false)
            addStatus('Whole-file IRF unavailable. Load a PTU with native TCSPC and try again.');
            return;
        end
        if ~isempty(app.irfGlobalDtNs) && abs(app.tcspc.dtNs - app.irfGlobalDtNs) > 1e-6
            addStatus('Whole-file IRF resolution does not match ROI TCSPC resolution.');
            return;
        end

        irfROI = app.irfGlobal(:);
        irfMeta = app.irfGlobalMeta;

        [tauFit, coeff, fitCountsRawFull, fitCountsShift, edgesRawNs] = fitRoiTCSPCWithIRF(...
            app.tcspc.rawCountsNative, app.tcspc.shiftBins, app.tcspc.dtNs, tau0, irfROI, includeBG, optimizeTau);
        [tFitNs, fitDensity, fitCountsLog, fitEdgesNs] = logBinTCSPC(fitCountsShift, app.tcspc.dtNs, app.txtBins.Value);

        app.tcspcFit = struct('tauFit', tauFit, 'coeff', coeff, 'fitCountsRaw', fitCountsRawFull, ...
            'fitCountsShift', fitCountsShift, 'fitCountsLog', fitCountsLog, 'fitDensity', fitDensity, ...
            'edgesRawNs', edgesRawNs, 'fitEdgesNs', fitEdgesNs, 'irf', irfROI(:), 'irfMeta', irfMeta, ...
            'usedCache', true);

        dataDensity = app.tcspc.cDensity(:);
        fitDensity = fitDensity(:);
        n = min(numel(dataDensity), numel(fitDensity));
        residualDensity = dataDensity(1:n) - fitDensity(1:n);

        cla(app.axTCSPC);
        yyaxis(app.axTCSPC, 'left');
        semilogx(app.axTCSPC, app.tcspc.tBinNs*1e-9, app.tcspc.cDensity, 'k.', 'MarkerSize', 10); hold(app.axTCSPC, 'on');
        semilogx(app.axTCSPC, tFitNs*1e-9, fitDensity, 'r-', 'LineWidth', 1.3);

        legendLines = {'ROI TCSPC', 'Reconvolution fit'};
        irfDisp = irfROI(:);
        if app.tcspc.shiftBins < numel(irfDisp)
            irfDisp = irfDisp(app.tcspc.shiftBins+1:end);
        end
        irfDisp = max(irfDisp, 0);
        if ~isempty(irfDisp) && max(irfDisp) > 0
            irfScale = max(app.tcspc.rawCountsShift) / max(irfDisp);
            irfScaled = irfDisp * irfScale;
            [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, app.tcspc.dtNs, app.txtBins.Value);
            semilogx(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
            legendLines{end+1} = 'Whole-file IRF (scaled)';
        end
        hold(app.axTCSPC, 'off');
        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Delay time (s)');
        ylabel(app.axTCSPC, 'Counts / ns');
        app.axTCSPC.YScale = 'log';
        title(app.axTCSPC, sprintf('ROI TCSPC fit using whole-file %s IRF', app.dropIRF.Value));

        yyaxis(app.axTCSPC, 'right');
        semilogx(app.axTCSPC, app.tcspc.tBinNs(1:n)*1e-9, residualDensity, 'm-', 'LineWidth', 1.0);
        hold(app.axTCSPC, 'on');
        yline(app.axTCSPC, 0, 'k-');
        hold(app.axTCSPC, 'off');
        ylabel(app.axTCSPC, 'Residual (counts/ns)');
        app.axTCSPC.YAxis(2).Visible = 'on';
        app.axTCSPC.YAxis(2).Scale = 'linear';
        yyaxis(app.axTCSPC, 'left');

        nParams = numel(tauFit) + double(includeBG);
        chi2 = sum((app.tcspc.cBin(:) - fitCountsLog(:)).^2 ./ max(fitCountsLog(:), 1));
        ndf = max(numel(app.tcspc.cBin) - nParams, 1);
        chi2red = chi2 / ndf;
        amps = coeff((1+double(includeBG)):end);
        ampFrac = amps ./ max(sum(amps), eps);

        setLegendEntries(legendLines);
        setFitSummary({sprintf('ROI fit | model: %d-exp', numel(tauFit)), ...
                       sprintf('tau (ns): [%s]', strtrim(num2str(tauFit(:).', '%.4g '))), ...
                       sprintf('amp frac: [%s]', strtrim(num2str(ampFrac(:).', '%.4g '))), ...
                       sprintf('chi2/ndf = %.4g', chi2red), ...
                       sprintf('IRF model: %s', app.dropIRF.Value), ...
                       sprintf('IRF source: whole file | dt = %.3f ps', 1e3*app.tcspc.dtNs)});

        app.tcspcDisplayMode = 'roi_fit';
        addStatus(sprintf('ROI fit tau (ns): %s | whole-file IRF model: %s', num2str(tauFit(:).', '%.3g '), app.dropIRF.Value));
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
        tcspc = app.tcspc; %#ok<NASGU>
        tcspcFit = app.tcspcFit; %#ok<NASGU>
        save(outFile, 'ptuOut', 'flim', 'ismRes', 'pattern', 'tcspc', 'tcspcFit');
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
            addStatus('No tcspc_pix. Reload with "Store TCSPC".');
            return;
        end

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5 5".');
            return;
        end

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        includeBG = app.chkIncludeBG.Value;

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

        if ~ensureGlobalIRF(false)
            addStatus('Whole-file IRF unavailable. Cannot run pattern matching.');
            return;
        end
        if ~isempty(app.irfGlobalDtNs) && abs(app.irfGlobalDtNs - dtNs) > 1e-6
            addStatus('Whole-file IRF resolution does not match pattern-matching TCSPC resolution.');
            return;
        end
        irf = app.irfGlobal(:);
        irfMeta = app.irfGlobalMeta;

        addStatus('Running global fit + pattern match with whole-file IRF...');
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
            outFLIM = GlobalMultiExpPatternMatchFromTCSPC(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0, opts);
            outFLIM.irf = irf(:);
            outFLIM.irfMeta = irfMeta;
            app.pattern = outFLIM;
        catch ME
            addStatus(['Pattern match failed: ' ME.message]);
            return;
        end

        try
            cmaps = {'gray','red','green','blue','orange','cyan','magenta'};
            nMaps = size(outFLIM.Amp, 3);
            cmaps = cmaps(1:min(nMaps, numel(cmaps)));
            ShowImagesMultiCmapWithColorbars(outFLIM.Amp, cmaps, [], [], true);
        catch
        end

        showPatternOverlay();
        addStatus(sprintf('Pattern match complete using whole-file %s IRF. GUI Tau0 was used only for fit components.', app.dropIRF.Value));
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
        Ncasc = ceil(log2(maxtime / Timeunit / Nsub));

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
            timegates = eval(app.editTimegates.Value); %#ok<EVLDIR>
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

    function addStatus(msg)
        app.status.Text = msg;
    end

    function clearTCSPCSidePanel()
        setFitSummary({'No fit results yet.'});
        setLegendEntries({'No legend yet.'});
    end

    function setFitSummary(lines)
        if nargin < 1 || isempty(lines)
            lines = {'No fit results yet.'};
        end
        if ischar(lines)
            lines = cellstr(lines);
        end
        try
            app.txtFitResults.Value = lines(:);
        catch
        end
    end

    function setLegendEntries(lines)
        if nargin < 1 || isempty(lines)
            lines = {'No legend yet.'};
        end
        if ischar(lines)
            lines = cellstr(lines);
        end
        try
            app.txtLegend.Value = lines(:);
        catch
        end
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
        imgG = applyGamma(double(img), app.gamma);
        cla(app.axImage);
        imagesc(app.axImage, imgG);
        axis(app.axImage, 'image');
        app.axImage.Visible = 'on';
        colormap(app.axImage, gray);
        colorbar(app.axImage);
        title(app.axImage, 'Intensity (gamma)');
        if ~isempty(app.cbTau) && isvalid(app.cbTau)
            delete(app.cbTau);
            app.cbTau = [];
        end
        showRecordedChannelGallery('intensity');
        app.displayMode = 'intensity';
        drawnow;
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
        showRecordedChannelGallery('tau');
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

    function showGlobalTCSPC()
        if isempty(app.ptuOut)
            return;
        end

        [countsNative, dtNs, srcInfo] = tcspcFromWholeFile_native(app.ptuOut);
        if isempty(countsNative)
            return;
        end

        if isfield(srcInfo, 'binned') && srcInfo.binned
            if ~isfield(srcInfo, 'tAxisNs') || isempty(srcInfo.tAxisNs)
                return;
            end
            tBinNs = srcInfo.tAxisNs(:);
            cInt = double(countsNative(:));
            if isfield(srcInfo, 'widthNs') && ~isempty(srcInfo.widthNs)
                widthNs = double(srcInfo.widthNs(:));
            else
                widthNs = ones(size(cInt));
            end
            cDensity = cInt ./ max(widthNs, eps);
            app.tcspcGlobal = struct('tBinNs', tBinNs, 'cDensity', cDensity, 'cBin', cInt, ...
                'widthNs', widthNs, 'dtNs', [], 'rawCountsNative', countsNative(:), 'source', srcInfo);
        else
            binsPerOct = app.txtBins.Value;
            [countsShift, shiftBins, shiftNs, peakIdx, displayStartIdx] = shiftDecayForIRFDisplay(countsNative, dtNs);
            [tBinNs, cDensity, cInt, edgesNs, widthNs] = logBinTCSPC(countsShift, dtNs, binsPerOct);
            app.tcspcGlobal = struct('tBinNs', tBinNs, 'cDensity', cDensity, 'cBin', cInt, ...
                'edgesNs', edgesNs, 'widthNs', widthNs, 'dtNs', dtNs, 'rawCountsNative', countsNative(:), ...
                'rawCountsShift', countsShift(:), 'shiftBins', shiftBins, 'shiftNs', shiftNs, ...
                'peakIdx', peakIdx, 'displayStartIdx', displayStartIdx, 'source', srcInfo);
        end

        cla(app.axTCSPC);
        yyaxis(app.axTCSPC, 'left');
        semilogx(app.axTCSPC, app.tcspcGlobal.tBinNs*1e-9, app.tcspcGlobal.cDensity, 'k.', 'MarkerSize', 10); 
        hold(app.axTCSPC, 'on');
        legendLines = {'Whole-file TCSPC'};

        if ensureGlobalIRF(false) && isfield(app.tcspcGlobal, 'dtNs') && ~isempty(app.tcspcGlobal.dtNs) ...
                && abs(app.irfGlobalDtNs - app.tcspcGlobal.dtNs) < 1e-6
            irfDisp = app.irfGlobal(:);
            if isfield(app.tcspcGlobal, 'shiftBins') && app.tcspcGlobal.shiftBins < numel(irfDisp)
                irfDisp = irfDisp(app.tcspcGlobal.shiftBins+1:end);
            end
            irfDisp = max(irfDisp, 0);
            if ~isempty(irfDisp) && max(irfDisp) > 0
                irfScale = max(app.tcspcGlobal.cBin) / max(irfDisp);
                irfScaled = irfDisp * irfScale;
                [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, app.irfGlobalDtNs, app.txtBins.Value);
                semilogx(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
                legendLines{end+1} = 'Whole-file IRF (scaled)';
            end
        end

        residualDensity = [];
        if ~isempty(app.globalFit) && isfield(app.globalFit, 'fitDensity') && ~isempty(app.globalFit.fitDensity)
            semilogx(app.axTCSPC, app.globalFit.tBinNs*1e-9, app.globalFit.fitDensity, 'r-', 'LineWidth', 1.3);
            legendLines{end+1} = sprintf('Best fit (%d-exp)', app.globalFit.nExp);
            n = min(numel(app.tcspcGlobal.cDensity), numel(app.globalFit.fitDensity));
            residualDensity = app.tcspcGlobal.cDensity(1:n) - app.globalFit.fitDensity(1:n);
            setFitSummary(app.globalFit.summaryLines);
        else
            setFitSummary({sprintf('Whole-file source: %s', srcInfo.mode), ...
                           sprintf('IRF model: %s', app.dropIRF.Value), ...
                           'No automatic fit available.'});
        end
        hold(app.axTCSPC, 'off');

        grid(app.axTCSPC, 'on');
        xlabel(app.axTCSPC, 'Delay time (s)');
        ylabel(app.axTCSPC, 'Counts / ns');
        app.axTCSPC.YScale = 'log';
        title(app.axTCSPC, sprintf('Whole-file TCSPC + whole-file IRF (%s)', app.dropIRF.Value));

        yyaxis(app.axTCSPC, 'right');
        if ~isempty(residualDensity)
            semilogx(app.axTCSPC, app.tcspcGlobal.tBinNs(1:numel(residualDensity))*1e-9, residualDensity, 'm-', 'LineWidth', 1.0);
            hold(app.axTCSPC, 'on');
            yline(app.axTCSPC, 0, 'k-');
            hold(app.axTCSPC, 'off');
            ylabel(app.axTCSPC, 'Residual (counts/ns)');
            app.axTCSPC.YAxis(2).Visible = 'on';
            app.axTCSPC.YAxis(2).Scale = 'linear';
            legendLines{end+1} = 'Residual';
        else
            app.axTCSPC.YAxis(2).Visible = 'off';
            app.axTCSPC.YAxis(2).Ticks = [];
        end
        yyaxis(app.axTCSPC, 'left');

        setLegendEntries(legendLines);
        app.tcspcDisplayMode = 'global';
    end

    function ok = ensureGlobalIRF(force)
        if nargin < 1
            force = false;
        end
        ok = false;
        if isempty(app.ptuOut)
            return;
        end

        model = app.dropIRF.Value;
        if ~force && ~isempty(app.irfGlobal) && strcmp(app.irfGlobalModel, model)
            ok = true;
            return;
        end

        [countsNative, dtNs, srcInfo] = tcspcFromWholeFile_native(app.ptuOut);
        if isempty(countsNative) || isempty(dtNs)
            addStatus('Whole-file IRF unavailable (no native/coarse TCSPC).');
            app.irfGlobal = [];
            app.irfGlobalModel = model;
            return;
        end

        headIRF = app.ptuOut.head;
        headIRF.MeasDesc_Resolution = dtNs * 1e-9;
        tauGuess = getDefaultTau0ForIRF();
        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        src = sprintf('wholefile:%s', srcInfo.mode);

        try
            [irf, cacheHit, irfMeta] = estimateIRFFromDecayFast(headIRF, countsNative, ...
                tauGuess, model, useGPU, app.irfCache, src);
            app.irfCache = struct('key', irfMeta.cacheKey, 'irf', irf(:), 'meta', irfMeta);
        catch ME
            addStatus(['Whole-file IRF failed: ' ME.message]);
            app.irfGlobal = [];
            app.irfGlobalModel = model;
            return;
        end

        app.irfGlobal = irf(:);
        app.irfGlobalModel = model;
        app.irfGlobalMeta = irfMeta;
        app.irfGlobalDtNs = dtNs;
        ok = true;

        if cacheHit
            addStatus(sprintf('Whole-file IRF cache hit (%s).', model));
        else
            addStatus(sprintf('Whole-file IRF estimated from global TCSPC (%s).', model));
        end
    end

    function ok = autoFitGlobalTCSPC()
        ok = false;
        app.globalFit = [];
        if isempty(app.ptuOut) || isempty(app.irfGlobal)
            return;
        end
        [countsNative, dtNs, srcInfo] = tcspcFromWholeFile_native(app.ptuOut);
        if isempty(countsNative) || isempty(dtNs)
            return;
        end
        includeBG = app.chkIncludeBG.Value;
        try
            fitAll = fitWholeDecayModelSweep(countsNative, dtNs, app.irfGlobal, 4, includeBG);
        catch ME
            addStatus(['Automatic whole-file fit failed: ' ME.message]);
            return;
        end
        if isempty(fitAll)
            return;
        end
        [bestFit, summaryLines] = chooseBestDecayModel(fitAll);
        bestFit.summaryLines = summaryLines;
        bestFit.source = srcInfo;
        [countsShift, shiftBins, shiftNs, peakIdx, displayStartIdx] = shiftDecayForIRFDisplay(countsNative, dtNs); %#ok<ASGLU>
        [tFitNs, fitDensity] = logBinTCSPC(bestFit.fitCountsRawFull(shiftBins+1:end), dtNs, app.txtBins.Value);
        bestFit.tBinNs = tFitNs;
        bestFit.fitDensity = fitDensity;
        app.globalFit = bestFit;
        ok = true;
    end

    function addGlobalFitSummaryText()
        if isempty(app.globalFit) || ~isfield(app.globalFit, 'summaryLines') || isempty(app.globalFit.summaryLines)
            return;
        end
        text(app.axTCSPC, 1.02, 0.98, app.globalFit.summaryLines, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 9, ...
            'Interpreter', 'none', 'BackgroundColor', 'w', 'Clipping', 'off');
    end

    function img = getIntensityMap()
        img = [];

        if ~isempty(app.flim) && isfield(app.flim, 'total') && isfield(app.flim.total, 'tag') ...
                && ~isempty(app.flim.total.tag)
            img = double(app.flim.total.tag);
            return;
        end

        if isempty(app.ptuOut)
            return;
        end

        if isfield(app.ptuOut, 'tags') && ~isempty(app.ptuOut.tags)
            img = sum(double(app.ptuOut.tags), 3);
            return;
        end

        if isfield(app.ptuOut, 'tag') && ~isempty(app.ptuOut.tag)
            if ndims(app.ptuOut.tag) >= 4
                img = squeeze(sum(sum(double(app.ptuOut.tag), 3), 4));
            elseif ndims(app.ptuOut.tag) == 3
                img = sum(double(app.ptuOut.tag), 3);
            else
                img = double(app.ptuOut.tag);
            end
            return;
        end

        if isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix)
            tmp = sum(double(app.ptuOut.tcspc_pix), 3);
            if ndims(tmp) == 3
                img = sum(tmp, 3);
            else
                img = tmp;
            end
            return;
        end

        if isfield(app.ptuOut, 'tcspc_pix_mt') && ~isempty(app.ptuOut.tcspc_pix_mt)
            tmp = sum(double(app.ptuOut.tcspc_pix_mt), 3);
            if ndims(tmp) == 3
                img = sum(tmp, 3);
            else
                img = tmp;
            end
            return;
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

    function stack = getIntensityChannelStack()
        stack = [];
        nRecorded = getRecordedChannelCountFromFolder();

        if isempty(app.ptuOut)
            return;
        end

        if isfield(app.ptuOut, 'tags') && ~isempty(app.ptuOut.tags)
            stack = collapseImageChannels(double(app.ptuOut.tags), nRecorded);
            return;
        end

        if isfield(app.ptuOut, 'tag') && ~isempty(app.ptuOut.tag)
            stack = collapseImageChannels(double(app.ptuOut.tag), nRecorded);
            return;
        end

        if ~isempty(app.flim) && isfield(app.flim, 'channels') && ~isempty(app.flim.channels)
            stack = collectChannelField(app.flim.channels, 'tag');
            if ~isempty(stack)
                return;
            end
        end

        if isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix)
            stack = collapseTcspcChannels(double(app.ptuOut.tcspc_pix));
            return;
        end

        if isfield(app.ptuOut, 'tcspc_pix_mt') && ~isempty(app.ptuOut.tcspc_pix_mt)
            stack = collapseTcspcChannels(double(app.ptuOut.tcspc_pix_mt));
        end
    end

    function stack = getTauChannelStack()
        stack = [];
        if isempty(app.flim) || ~isfield(app.flim, 'channels') || isempty(app.flim.channels)
            return;
        end
        stack = collectChannelField(app.flim.channels, 'tauMean');
    end

    function showRecordedChannelGallery(mode)
        if nargin < 1 || isempty(mode)
            mode = 'intensity';
        end

        switch lower(mode)
            case 'tau'
                stack = getTauChannelStack();
                if getChannelCount(stack) <= 1
                    return;
                end
                refIntensity = collectChannelField(app.flim.channels, 'tag');
                trange = getTauRange(stack);
                cmap = jet(256);
                cmap = cmap(30:end-30,:);
                figTitle = sprintf('Per-channel tau maps (%d channels)', getChannelCount(stack));
            otherwise
                stack = getIntensityChannelStack();
                if getChannelCount(stack) <= 1
                    return;
                end
                figTitle = sprintf('Per-channel intensity maps (%d channels)', getChannelCount(stack));
        end

        nCh = getChannelCount(stack);
        labels = getChannelLabels(nCh);
        nCols = min(3, nCh);
        nRows = ceil(nCh / nCols);

        if isempty(app.channelGallery) || ~isgraphics(app.channelGallery)
            app.channelGallery = figure('Name', 'Luminosa Channel Gallery', 'NumberTitle', 'off', 'Color', 'w');
        else
            figure(app.channelGallery);
            clf(app.channelGallery);
        end

        tl = tiledlayout(app.channelGallery, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, figTitle);

        for iCh = 1:nCh
            ax = nexttile(tl);
            switch lower(mode)
                case 'tau'
                    tau = double(stack(:,:,iCh));
                    if size(refIntensity, 3) >= iCh
                        rgb = buildTauOverlayRGB(tau, refIntensity(:,:,iCh), trange);
                        image(ax, rgb);
                        axis(ax, 'image');
                        axis(ax, 'off');
                    else
                        imagesc(ax, tau, trange);
                        axis(ax, 'image');
                        axis(ax, 'off');
                        colormap(ax, cmap);
                        colorbar(ax);
                    end
                otherwise
                    imagesc(ax, applyGamma(double(stack(:,:,iCh)), app.gamma));
                    axis(ax, 'image');
                    axis(ax, 'off');
                    colormap(ax, gray);
                    colorbar(ax);
            end
            title(ax, labels{iCh}, 'Interpreter', 'none');
        end
    end

    function msg = getChannelDiagnosticMessage()
        msg = '';
        nRecorded = getRecordedChannelCountFromFolder();
        if nRecorded <= 1
            return;
        end

        nLoaded = getChannelCount(getIntensityChannelStack());
        if nLoaded >= nRecorded
            msg = sprintf('Detected %d recorded channels; all are shown in the channel gallery.', nRecorded);
        elseif nLoaded > 1
            msg = sprintf('Detected %d recorded channels on disk, but only %d were loaded into MATLAB.', nRecorded, nLoaded);
        else
            msg = sprintf('Detected %d recorded channels on disk, but the current reader exposed only one combined image.', nRecorded);
        end
    end

    function n = getRecordedChannelCountFromFolder()
        n = 0;
        if isempty(app.lastFile)
            return;
        end
        dataDir = fileparts(app.lastFile);
        listing = dir(fullfile(dataDir, 'IntensityImage*.pqdat'));
        if isempty(listing)
            listing = dir(fullfile(dataDir, 'FastLTImage*.pqdat'));
        end
        n = numel(listing);
    end

    function labels = getChannelLabels(nCh)
        labels = cell(1, nCh);
        for iCh = 1:nCh
            labels{iCh} = sprintf('Channel %d', iCh);
        end
    end

    function rgb = buildTauOverlayRGB(tau, intensity, trange)
        cmap = jet(256);
        cmap = cmap(30:end-30,:);
        bright = applyGamma(double(intensity), app.gamma);
        rgb = flim_rgb(double(tau), bright, trange, cmap);
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
        rgb = buildTauOverlayRGB(tau, intensity, trange);
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

function fitAll = fitWholeDecayModelSweep(countsRawFull, dtNs, irf, maxExp, includeBG)
    countsRawFull = double(countsRawFull(:));
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    maxExp = max(1, round(maxExp));
    fitAll = repmat(struct(), 0, 1);
    for nExp = 1:maxExp
        seeds = getDefaultTauSeedsForGlobalFit(nExp);
        [tauFit, coeff, fitCountsRawFull, stats] = fitDecayIRFModelMultiStart(countsRawFull, dtNs, seeds, irf, includeBG, true);
        fitAll(end+1,1).nExp = nExp; %#ok<AGROW>
        fitAll(end).tauFit = tauFit(:).';
        fitAll(end).coeff = coeff(:).';
        fitAll(end).fitCountsRawFull = fitCountsRawFull(:);
        fitAll(end).stats = stats;
    end
end

function [bestFit, summaryLines] = chooseBestDecayModel(fitAll)
    if isempty(fitAll)
        bestFit = [];
        summaryLines = {};
        return;
    end
    bicVals = arrayfun(@(s) s.stats.bic, fitAll);
    [minBIC, ~] = min(bicVals);
    cand = find(bicVals <= minBIC + 2);
    if isempty(cand)
        [~, idxBest] = min(bicVals);
    else
        idxBest = cand(1);
    end
    bestFit = fitAll(idxBest);

    amps = bestFit.coeff(1 + double(bestFit.stats.includeBG):end);
    ampFrac = amps ./ max(sum(amps), eps);
    tauStr = num2str(bestFit.tauFit, '%.3g ');
    ampStr = num2str(ampFrac(:).', '%.3g ');
    summaryLines = { ...
        sprintf('Best model: %d-exp', bestFit.nExp), ...
        sprintf('tau (ns): [%s]', strtrim(tauStr)), ...
        sprintf('amp frac: [%s]', strtrim(ampStr)), ...
        sprintf('chi2/ndf = %.4g', bestFit.stats.chi2red), ...
        sprintf('BIC-like = %.4g', bestFit.stats.bic), ...
        'Model sweep:'};
    for k = 1:numel(fitAll)
        summaryLines{end+1} = sprintf('  %d-exp: chi2/ndf=%.4g  BIC=%.4g', ...
            fitAll(k).nExp, fitAll(k).stats.chi2red, fitAll(k).stats.bic); %#ok<AGROW>
    end
end

function tauSeed = getDefaultTauSeedsForGlobalFit(nExp)
    switch nExp
        case 1
            tauSeed = 1.5;
        case 2
            tauSeed = [0.4 2.5];
        case 3
            tauSeed = [0.25 1.2 4.5];
        otherwise
            tauSeed = [0.15 0.6 2.0 6.0];
            tauSeed = tauSeed(1:nExp);
    end
end

function [tauFit, coeff, fitCountsRawFull, stats] = fitDecayIRFModelMultiStart(countsRawFull, dtNs, tau0, irf, includeBG, optimizeTau)
    tau0 = double(tau0(:));
    countsRawFull = double(countsRawFull(:));
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    seedMat = buildTauSeedMatrix(tau0);
    best.err = inf;
    best.tauFit = tau0(:).';
    best.coeff = [];
    best.fitCountsRawFull = [];
    best.stats = [];
    for iseed = 1:size(seedMat,2)
        tauSeed = seedMat(:,iseed);
        if optimizeTau
            p0 = tauSeed;
            xmin = max(0.03, p0 / 10);
            xmax = max(p0 * 10, p0 + 0.05);
            tol = 1e-5;
            steps = max(250, 180 * numel(p0));
            [pfit, ~] = Simplex(@roiTcspcErrRawIRF, p0, xmin, xmax, tol, steps, countsRawFull, dtNs, irf, includeBG);
        else
            pfit = tauSeed;
        end
        [err, coeffCand, fitCand] = roiTcspcErrRawIRF(pfit, countsRawFull, dtNs, irf, includeBG);
        [tauSorted, coeffSorted] = sortLifetimesAndAmps(pfit(:).', coeffCand(:).', includeBG);
        [err, coeffCand, fitCand] = roiTcspcErrRawIRF(tauSorted(:), countsRawFull, dtNs, irf, includeBG);
        statsCand = calcDecayFitStats(countsRawFull, fitCand, numel(tauSorted), includeBG);
        if err < best.err
            best.err = err;
            best.tauFit = tauSorted;
            best.coeff = coeffCand(:).';
            best.fitCountsRawFull = fitCand(:);
            best.stats = statsCand;
        end
    end
    tauFit = best.tauFit;
    coeff = best.coeff;
    fitCountsRawFull = best.fitCountsRawFull;
    stats = best.stats;
end

function seedMat = buildTauSeedMatrix(tau0)
    tau0 = double(tau0(:));
    n = numel(tau0);
    scales = [1.0, 0.75, 1.35];
    seedMat = zeros(n, numel(scales));
    for k = 1:numel(scales)
        s = scales(k);
        seed = tau0 * s;
        seed = max(seed, 0.03);
        seedMat(:,k) = seed;
    end
end

function [tauSorted, coeffSorted] = sortLifetimesAndAmps(tauFit, coeff, includeBG)
    tauSorted = double(tauFit(:).');
    coeffSorted = double(coeff(:).');
    if isempty(tauSorted)
        return;
    end
    [tauSorted, ord] = sort(tauSorted, 'ascend');
    if includeBG
        coeffSorted = [coeffSorted(1), coeffSorted(1 + ord)];
    else
        coeffSorted = coeffSorted(ord);
    end
end

function stats = calcDecayFitStats(counts, fitCounts, nExp, includeBG)
    counts = max(double(counts(:)), 0);
    fitCounts = max(double(fitCounts(:)), eps);
    n = numel(counts);
    k = nExp + double(includeBG);
    chi2 = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
    ndf = max(n - k, 1);
    chi2red = chi2 / ndf;
    bic = chi2 + k * log(max(n,2));
    aicc = chi2 + 2*k + (2*k*(k+1)) / max(n-k-1, 1);
    stats = struct('chi2', chi2, 'chi2red', chi2red, 'bic', bic, 'aicc', aicc, 'n', n, 'k', k, 'includeBG', logical(includeBG));
end

function ok = gpuIsAvailable()
    ok = false;
    try
        ok = gpuDeviceCount > 0;
    catch
        ok = false;
    end
end

function flim = quickFLIMFromTCSPCFlexible(ptuOut, useGPU)
% Use either dense linear tcspc_pix or compressed tcspc_pix_mt.
    if isfield(ptuOut, 'tcspc_pix') && ~isempty(ptuOut.tcspc_pix)
        cube = ptuOut.tcspc_pix;
        dt_ns = ptuOut.head.MeasDesc_Resolution * 1e9;
        tAxisNs = ((1:size(cube, 3)) - 0.5) * dt_ns;
    elseif isfield(ptuOut, 'tcspc_pix_mt') && ~isempty(ptuOut.tcspc_pix_mt)
        cube = ptuOut.tcspc_pix_mt;
        if isfield(ptuOut, 'tcspc_mt_centers_ns') && ~isempty(ptuOut.tcspc_mt_centers_ns)
            tAxisNs = double(ptuOut.tcspc_mt_centers_ns(:)');
        else
            error('tcspc_pix_mt exists but tcspc_mt_centers_ns is missing.');
        end
    else
        error('No TCSPC cube available.');
    end

    if useGPU
        cubeG = gpuArray(double(cube));
        if ndims(cubeG) == 4
            globalDecay = gather(squeeze(sum(sum(sum(cubeG, 1), 2), 4)));
            tag = gather(sum(sum(cubeG, 3), 4));
        else
            globalDecay = gather(squeeze(sum(sum(cubeG, 1), 2)));
            tag = gather(sum(cubeG, 3));
        end
    else
        cubeG = double(cube);
        if ndims(cubeG) == 4
            globalDecay = squeeze(sum(sum(sum(cubeG, 1), 2), 4));
            tag = sum(sum(cubeG, 3), 4);
        else
            globalDecay = squeeze(sum(sum(cubeG, 1), 2));
            tag = sum(cubeG, 3);
        end
    end

    globalDecay = double(globalDecay(:));
    [~, t0Bin] = max(globalDecay);
    tShift = double(tAxisNs(:)) - double(tAxisNs(t0Bin));

    if useGPU
        tShift3 = gpuArray(reshape(tShift, 1, 1, []));
        if ndims(cubeG) == 4
            meanArrivalNum = sum(sum(cubeG .* tShift3, 3), 4);
            denom = max(sum(sum(cubeG, 3), 4), 1);
        else
            meanArrivalNum = sum(cubeG .* tShift3, 3);
            denom = max(sum(cubeG, 3), 1);
        end
        meanArrival = gather(meanArrivalNum ./ denom);
    else
        tShift3 = reshape(tShift, 1, 1, []);
        if ndims(cubeG) == 4
            meanArrivalNum = sum(sum(cubeG .* tShift3, 3), 4);
            denom = max(sum(sum(cubeG, 3), 4), 1);
        else
            meanArrivalNum = sum(cubeG .* tShift3, 3);
            denom = max(sum(cubeG, 3), 1);
        end
        meanArrival = meanArrivalNum ./ denom;
    end

    flim = struct();
    flim.total = struct();
    flim.total.tag = double(tag);
    flim.total.tauMean = double(meanArrival);
    flim.total.tauRMS = zeros(size(meanArrival));
    flim.total.meanArrival = double(meanArrival);
    flim.total.globalDecay = double(globalDecay(:));
    flim.total.tAxisNs = double(tAxisNs(:));
    flim.total.t0Bin = t0Bin;
    flim.channels = repmat(struct('tag', [], 'tauMean', [], 'tauRMS', [], ...
        'meanArrival', [], 'globalDecay', [], 'tAxisNs', [], 't0Bin', []), 0, 1);

    if ndims(cubeG) == 4
        nCh = size(cubeG, 4);
        flim.channels = repmat(struct('tag', [], 'tauMean', [], 'tauRMS', [], ...
            'meanArrival', [], 'globalDecay', [], 'tAxisNs', [], 't0Bin', t0Bin), nCh, 1);
        for iCh = 1:nCh
            cubeCh = cubeG(:,:,:,iCh);
            if useGPU
                tagCh = gather(sum(cubeCh, 3));
                globalDecayCh = gather(squeeze(sum(sum(cubeCh, 1), 2)));
                meanArrivalNumCh = sum(cubeCh .* tShift3, 3);
                denomCh = max(sum(cubeCh, 3), 1);
                meanArrivalCh = gather(meanArrivalNumCh ./ denomCh);
            else
                tagCh = sum(cubeCh, 3);
                globalDecayCh = squeeze(sum(sum(cubeCh, 1), 2));
                meanArrivalNumCh = sum(cubeCh .* tShift3, 3);
                denomCh = max(sum(cubeCh, 3), 1);
                meanArrivalCh = meanArrivalNumCh ./ denomCh;
            end

            flim.channels(iCh).tag = double(tagCh);
            flim.channels(iCh).tauMean = double(meanArrivalCh);
            flim.channels(iCh).tauRMS = zeros(size(meanArrivalCh));
            flim.channels(iCh).meanArrival = double(meanArrivalCh);
            flim.channels(iCh).globalDecay = double(globalDecayCh(:));
            flim.channels(iCh).tAxisNs = double(tAxisNs(:));
            flim.channels(iCh).t0Bin = t0Bin;
        end
    else
        flim.channels = flim.total;
    end
end

function stack = collapseImageChannels(data, nRecorded)
    if nargin < 2
        nRecorded = 0;
    end

    switch ndims(data)
        case 2
            stack = reshape(double(data), size(data,1), size(data,2), 1);
        case 3
            stack = double(data);
        case 4
            sz = size(data);
            if nRecorded > 0 && sz(4) == nRecorded
                stack = squeeze(sum(double(data), 3));
            elseif nRecorded > 0 && sz(3) == nRecorded
                stack = squeeze(sum(double(data), 4));
            elseif sz(4) <= sz(3)
                stack = squeeze(sum(double(data), 3));
            else
                stack = squeeze(sum(double(data), 4));
            end
            if ndims(stack) == 2
                stack = reshape(stack, size(stack,1), size(stack,2), 1);
            end
        otherwise
            stack = [];
    end
end

function stack = collapseTcspcChannels(cube)
    if ndims(cube) == 4
        stack = squeeze(sum(double(cube), 3));
    elseif ndims(cube) == 3
        img = sum(double(cube), 3);
        stack = reshape(img, size(img,1), size(img,2), 1);
    else
        stack = [];
    end
end

function stack = collectChannelField(channelStruct, fieldName)
    stack = [];
    if isempty(channelStruct)
        return;
    end

    validIdx = find(arrayfun(@(s) isfield(s, fieldName) && ~isempty(s.(fieldName)), channelStruct));
    if isempty(validIdx)
        return;
    end

    firstVal = double(channelStruct(validIdx(1)).(fieldName));
    stack = zeros([size(firstVal) numel(validIdx)]);
    for i = 1:numel(validIdx)
        stack(:,:,i) = double(channelStruct(validIdx(i)).(fieldName));
    end
end

function n = getChannelCount(stack)
    if isempty(stack)
        n = 0;
    elseif ndims(stack) < 3
        n = 1;
    else
        n = size(stack, 3);
    end
end

function [counts, dtNs, srcInfo] = tcspcFromROI_native(ptuOut, roi)
% Prefer per-photon native TCSPC when available. Fall back to tcspc_pix.
    srcInfo = struct('mode', 'unknown', 'native', false);

    pos = roi.Position;
    x1 = max(1, floor(pos(1)));
    y1 = max(1, floor(pos(2)));
    x2 = ceil(pos(1) + pos(3));
    y2 = ceil(pos(2) + pos(4));

    if isfield(ptuOut, 'im_col') && isfield(ptuOut, 'im_line')
        col = double(ptuOut.im_col(:));
        row = double(ptuOut.im_line(:));
        idx = (col >= x1) & (col <= x2) & (row >= y1) & (row <= y2);
    else
        idx = [];
    end

    % Native per-photon TCSPC preserved by the patched reader
    if ~isempty(idx) && any(idx) && isfield(ptuOut, 'im_tcspc_native') && ~isempty(ptuOut.im_tcspc_native)
        dtNs = getNativeTcspcResolutionNs(ptuOut.head);
        tNative = double(ptuOut.im_tcspc_native(idx));
        if isempty(tNative)
            counts = [];
            return;
        end
        nT = max(max(tNative), ceil((ptuOut.head.MeasDesc_GlobalResolution * 1e9) / dtNs) + 1);
        counts = accumarray(max(1, round(tNative(:))), 1, [nT, 1], @sum, 0);
        srcInfo.mode = 'native photon list';
        srcInfo.native = true;
        return;
    end

    % Fallback to native-like photon list if resolution fields are present
    if ~isempty(idx) && any(idx) && isfield(ptuOut, 'im_tcspc') && ~isempty(ptuOut.im_tcspc) && isfield(ptuOut.head, 'MeasDesc_Resolution_Original')
        dtNs = getNativeTcspcResolutionNs(ptuOut.head);
        chDiv = 1;
        if isfield(ptuOut.head, 'TCSPC_chDiv') && ~isempty(ptuOut.head.TCSPC_chDiv)
            chDiv = double(ptuOut.head.TCSPC_chDiv);
        end
        tCoarse = double(ptuOut.im_tcspc(idx));
        tNative = 1 + (tCoarse - 1) * chDiv;
        nT = max(max(tNative), ceil((ptuOut.head.MeasDesc_GlobalResolution * 1e9) / dtNs) + 1);
        counts = accumarray(max(1, round(tNative(:))), 1, [nT, 1], @sum, 0);
        srcInfo.mode = 'expanded coarse photon list';
        srcInfo.native = false;
        return;
    end

    % Final fallback: cube resolution
    if isfield(ptuOut, 'tcspc_pix') && ~isempty(ptuOut.tcspc_pix)
        cube = ptuOut.tcspc_pix;
        dtNs = ptuOut.head.MeasDesc_Resolution * 1e9;
        x2 = min(size(cube, 2), x2);
        y2 = min(size(cube, 1), y2);
        sub = cube(y1:y2, x1:x2, :, :);
        counts = squeeze(sum(sum(sum(sub, 4), 2), 1));
        counts = double(counts(:));
        srcInfo.mode = 'tcspc cube';
        srcInfo.native = false;
        return;
    end

    counts = [];
    dtNs = [];
end

function [counts, dtNs, srcInfo] = tcspcFromWholeFile_native(ptuOut)
% Whole-file TCSPC with best available resolution.
    srcInfo = struct('mode', 'unknown', 'native', false, 'binned', false, 'tAxisNs', [], 'widthNs', []);

    if isfield(ptuOut, 'im_tcspc_native') && ~isempty(ptuOut.im_tcspc_native)
        dtNs = getNativeTcspcResolutionNs(ptuOut.head);
        tNative = double(ptuOut.im_tcspc_native(:));
        if isempty(tNative)
            counts = [];
            return;
        end
        nT = max(max(tNative), ceil((ptuOut.head.MeasDesc_GlobalResolution * 1e9) / dtNs) + 1);
        counts = accumarray(max(1, round(tNative(:))), 1, [nT, 1], @sum, 0);
        srcInfo.mode = 'native photon list';
        srcInfo.native = true;
        return;
    end

    if isfield(ptuOut, 'im_tcspc') && ~isempty(ptuOut.im_tcspc) && isfield(ptuOut.head, 'MeasDesc_Resolution_Original')
        dtNs = getNativeTcspcResolutionNs(ptuOut.head);
        chDiv = 1;
        if isfield(ptuOut.head, 'TCSPC_chDiv') && ~isempty(ptuOut.head.TCSPC_chDiv)
            chDiv = double(ptuOut.head.TCSPC_chDiv);
        end
        tCoarse = double(ptuOut.im_tcspc(:));
        tNative = 1 + (tCoarse - 1) * chDiv;
        nT = max(max(tNative), ceil((ptuOut.head.MeasDesc_GlobalResolution * 1e9) / dtNs) + 1);
        counts = accumarray(max(1, round(tNative(:))), 1, [nT, 1], @sum, 0);
        srcInfo.mode = 'expanded coarse photon list';
        srcInfo.native = false;
        return;
    end

    if isfield(ptuOut, 'tcspc_pix') && ~isempty(ptuOut.tcspc_pix)
        cube = ptuOut.tcspc_pix;
        dtNs = ptuOut.head.MeasDesc_Resolution * 1e9;
        counts = squeeze(sum(sum(sum(cube, 4), 2), 1));
        counts = double(counts(:));
        srcInfo.mode = 'tcspc cube';
        srcInfo.native = false;
        return;
    end

    if isfield(ptuOut, 'tcspc_pix_mt') && ~isempty(ptuOut.tcspc_pix_mt)
        cube = ptuOut.tcspc_pix_mt;
        if ndims(cube) == 4
            counts = squeeze(sum(sum(sum(cube, 4), 2), 1));
        else
            counts = squeeze(sum(sum(cube, 2), 1));
        end
        counts = double(counts(:));
        dtNs = [];
        srcInfo.mode = 'multitau cube';
        srcInfo.native = false;
        srcInfo.binned = true;
        if isfield(ptuOut, 'tcspc_mt_centers_ns') && ~isempty(ptuOut.tcspc_mt_centers_ns)
            srcInfo.tAxisNs = double(ptuOut.tcspc_mt_centers_ns(:));
        end
        if isfield(ptuOut, 'tcspc_mt_width_ns') && ~isempty(ptuOut.tcspc_mt_width_ns)
            srcInfo.widthNs = double(ptuOut.tcspc_mt_width_ns(:));
        end
        return;
    end

    counts = [];
    dtNs = [];
end

function dtNs = getNativeTcspcResolutionNs(head)
    if isfield(head, 'MeasDesc_Resolution_Original') && ~isempty(head.MeasDesc_Resolution_Original)
        dtNs = double(head.MeasDesc_Resolution_Original) * 1e9;
    else
        dtNs = double(head.MeasDesc_Resolution) * 1e9;
    end
end

function [countsShift, shiftBins, shiftNs, peakIdx, startIdx] = shiftDecayForIRFDisplay(counts, dtNs)
% Shift the decay for display, but keep a short pre-peak region so the IRF
% leading edge remains visible on the plot. Fitting still uses the full decay.
    counts = double(counts(:));
    if isempty(counts)
        countsShift = counts;
        shiftBins = 0;
        shiftNs = 0;
        peakIdx = 1;
        startIdx = 1;
        return;
    end

    [~, peakIdx] = max(counts);
    peakIdx = max(1, peakIdx);

    prePeakNs = 0.30;
    prePeakBins = max(8, round(prePeakNs / max(dtNs, eps)));
    prePeakBins = min(prePeakBins, peakIdx - 1);

    startIdx = max(1, peakIdx - prePeakBins);
    countsShift = counts(startIdx:end);
    shiftBins = startIdx - 1;
    shiftNs = shiftBins * dtNs;
end

function [tBinNs, cDensity, cInt, edgesNs, widthNs] = logBinTCSPC(counts, dtNs, binsPerOct)
% Logarithmic binning for decay display after left-shifting.
% Display quantity is count density (counts / ns), not integrated counts.
    if nargin < 3 || isempty(binsPerOct)
        binsPerOct = 8;
    end
    binsPerOct = max(2, round(binsPerOct));

    counts = double(counts(:));
    Nt = numel(counts);
    if Nt < 2
        edgesNs = [0; dtNs];
        widthNs = diff(edgesNs);
        cInt = counts(:);
        cDensity = cInt ./ max(widthNs, eps);
        tBinNs = dtNs / 2;
        return;
    end

    nOct = max(1, log2(Nt));
    nBins = max(16, ceil(nOct * binsPerOct));
    edgesIdx = unique(round(logspace(0, log10(Nt + 1), nBins + 1)));
    edgesIdx(1) = 1;
    edgesIdx(end) = Nt + 1;
    edgesIdx = unique(edgesIdx(:));

    nBins = numel(edgesIdx) - 1;
    cInt = zeros(nBins, 1);
    widthNs = zeros(nBins, 1);
    tBinNs = zeros(nBins, 1);

    csum = cumsum([0; counts]);
    for i = 1:nBins
        a = edgesIdx(i);
        b = edgesIdx(i+1) - 1;
        cInt(i) = csum(b + 1) - csum(a);
        widthNs(i) = (b - a + 1) * dtNs;
        t0 = (a - 1) * dtNs;
        t1 = b * dtNs;
        if t0 <= 0
            tBinNs(i) = max(dtNs/2, t1/2);
        else
            tBinNs(i) = sqrt(t0 * t1);
        end
    end

    edgesNs = (edgesIdx(:) - 1) * dtNs;
    cDensity = cInt ./ max(widthNs, eps);
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

function tau0 = getDefaultTau0ForIRF()
% Internal defaults used only for whole-file IRF estimation.
    tau0 = [0.35 1.5 5];
end

function [tauFit, coeff, fitCountsRawFull, fitCountsShift, edgesRawNs] = fitRoiTCSPCWithIRF(countsRawFull, shiftBins, dtNs, tau0, irf, includeBG, optimizeTau)
% Fit full native-resolution ROI decay with IRF reconvolution.
% The displayed curve is shifted only after fitting so the IRF remains meaningful.
    tau0 = double(tau0(:));
    countsRawFull = double(countsRawFull(:));
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);
    edgesRawNs = (0:numel(countsRawFull))' * dtNs;

    if optimizeTau
        p0 = tau0;
        xmin = max(0.05, p0 / 10);
        xmax = p0 * 10;
        tol = 1e-5;
        steps = 250 * numel(p0);
        [pfit, ~] = Simplex(@roiTcspcErrRawIRF, p0, xmin, xmax, tol, steps, countsRawFull, dtNs, irf, includeBG);
    else
        pfit = tau0;
    end

    [~, coeff, fitCountsRawFull] = roiTcspcErrRawIRF(pfit, countsRawFull, dtNs, irf, includeBG);
    tauFit = pfit(:).';

    shiftBins = max(0, round(shiftBins));
    if shiftBins >= numel(fitCountsRawFull)
        fitCountsShift = fitCountsRawFull(:);
    else
        fitCountsShift = fitCountsRawFull(shiftBins+1:end);
    end
end

function [err, coeff, fitCounts] = roiTcspcErrRawIRF(tau, counts, dtNs, irf, includeBG)
    tau = max(double(tau(:)), 1e-6);
    counts = max(double(counts(:)), 0);
    irf = max(double(irf(:)), 0);
    irf = irf ./ max(sum(irf), eps);

    nExp = numel(tau);
    nBins = numel(counts);
    t0 = ((0:nBins-1)') * dtNs;
    t1 = ((1:nBins)') * dtNs;

    M = zeros(nBins, nExp + includeBG);
    if includeBG
        M(:, 1) = 1;
    end

    for k = 1:nExp
        tk = tau(k);
        decayBin = tk * (exp(-t0 ./ tk) - exp(-t1 ./ tk));
        convSig = conv(irf, decayBin, 'full');
        M(:, includeBG + k) = convSig(1:nBins);
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

function [irf, cacheHit, meta] = estimateIRFFast(head, tcspc_pix, tau0, irfModel, useGPU, irfCache, src)
    globalDecay = collapseGlobalDecayFast(tcspc_pix, useGPU);
    globalDecay = double(globalDecay(:));
    [irf, cacheHit, meta] = estimateIRFFromDecayFast(head, globalDecay, tau0, irfModel, useGPU, irfCache, src);
end

function [irf, cacheHit, meta] = estimateIRFFromDecayFast(head, decay, tau0, irfModel, useGPU, irfCache, src)
    tau0 = double(tau0(:)).';
    dtNs = head.MeasDesc_Resolution * 1e9;
    pulsePeriodNs = head.MeasDesc_GlobalResolution * 1e9;
    decay = max(double(decay(:)), 0);

    [~, peakIdx] = max(decay);
    decaySig = [numel(decay), sum(decay), peakIdx, sum((1:numel(decay))'.*decay) / max(sum(decay), eps)];
    cacheKey = sprintf('%s|%s|%s|%s|%.12g|%.12g|%d|%.12g|%.12g', ...
        irfModel, sprintf('%.6g_', tau0), src, mat2str(size(decay)), ...
        decaySig(2), decaySig(1), decaySig(3), decaySig(4), dtNs);
    cacheHit = isstruct(irfCache) && isfield(irfCache, 'key') && strcmp(irfCache.key, cacheKey) && isfield(irfCache, 'irf') && ~isempty(irfCache.irf);
    if cacheHit
        irf = irfCache.irf(:);
        meta = struct('cacheKey', cacheKey, 'model', irfModel, 'usedGPUCollapse', logical(useGPU), 'nCounts', decaySig(2), 'source', src);
        return;
    end

    switch irfModel
        case 'Simple'
            head2 = head;
            head2.Resolution = dtNs;
            if isfield(head, 'TTResult_SyncRate')
                head2.SyncRate = head.TTResult_SyncRate;
            elseif isfield(head, 'SyncRate')
                head2.SyncRate = head.SyncRate;
            end
            irf = Calc_mIRF(head2, decay');
            irf = double(irf(:));

        case 'GammaShifted'
            try
                irf = fitIRFfromGlobalDecay(decay, dtNs, pulsePeriodNs, tau0, 'GammaShifted');
            catch
                try
                    outIRF = Calc_mIRF_Global_GammaShifted_fast(head, decay, tau0, struct('nCasc',4,'nSub',6));
                catch
                    outIRF = Calc_mIRF_Global_GammaShifted(head, decay, tau0);
                end
                irf = double(outIRF.IRF(:));
            end

        case 'ExGauss'
            try
                irf = fitIRFfromGlobalDecay(decay, dtNs, pulsePeriodNs, tau0, 'ExGauss');
            catch
                try
                    outIRF = Calc_mIRF_Global_ExGauss_fast(head, decay, tau0, struct('nCasc',4,'nSub',6));
                catch
                    outIRF = Calc_mIRF_Global_ExGauss(head, decay, tau0);
                end
                irf = double(outIRF.IRF(:));
            end

        otherwise
            error('Unknown IRF model: %s', irfModel);
    end

    irf = max(irf(:), 0);
    s = sum(irf);
    if s <= 0
        error('Estimated IRF is empty or non-positive.');
    end
    irf = irf / s;
    meta = struct('cacheKey', cacheKey, 'model', irfModel, 'usedGPUCollapse', logical(useGPU), 'nCounts', decaySig(2), 'source', src);
end

function globalDecay = collapseGlobalDecayFast(tcspc_pix, useGPU)
% Sum x, y, and detector/channel dimensions down to one decay vector.
    if ndims(tcspc_pix) == 4
        if useGPU
            g = gpuArray(single(tcspc_pix));
            globalDecay = gather(squeeze(sum(sum(sum(g, 1), 2), 4)));
        else
            globalDecay = squeeze(sum(sum(sum(double(tcspc_pix), 1), 2), 4));
        end
    elseif ndims(tcspc_pix) == 3
        if useGPU
            g = gpuArray(single(tcspc_pix));
            globalDecay = gather(squeeze(sum(sum(g, 1), 2)));
        else
            globalDecay = squeeze(sum(sum(double(tcspc_pix), 1), 2));
        end
    else
        error('Unsupported tcspc_pix dimensionality.');
    end
    globalDecay = double(globalDecay(:));
end

function irf = fitIRFfromGlobalDecay(globalDecay, dtNs, pulsePeriodNs, tau0, modelName)
% Fast 1D IRF estimation from the global decay.
% This is where the speedup comes from: fitting a 1D curve is much cheaper
% than repeatedly walking the whole 4D TCSPC cube.
    y = max(double(globalDecay(:)), 0);
    n = numel(y);
    t = ((0:n-1)' + 0.5) * dtNs;

    [~, peakIdx] = max(y);
    tPeak = t(peakIdx);

    switch modelName
        case 'GammaShifted'
            p0 = [max(0, tPeak - 2*dtNs), log(2.5), log(max(dtNs, 0.12 * pulsePeriodNs / n))];
            f = @(p) irfFitErrGamma(p, y, t, tau0);
            p = fminsearch(f, p0, optimset('Display','off', 'MaxIter', 200));
            irf = makeGammaShiftedIRF(t, p(1), exp(p(2)), exp(p(3)));
        case 'ExGauss'
            p0 = [tPeak, log(max(dtNs, 2*dtNs)), log(max(dtNs, 3*dtNs))];
            f = @(p) irfFitErrExGauss(p, y, t, tau0);
            p = fminsearch(f, p0, optimset('Display','off', 'MaxIter', 250));
            irf = makeExGaussIRF(t, p(1), exp(p(2)), exp(p(3)));
        otherwise
            error('Unsupported model: %s', modelName);
    end

    irf = max(irf(:), 0);
    irf = irf ./ max(sum(irf), eps);
end

function err = irfFitErrGamma(p, y, t, tau0)
    shift = p(1);
    k = exp(p(2));
    theta = exp(p(3));
    if ~isfinite(shift) || ~isfinite(k) || ~isfinite(theta) || k <= 0 || theta <= 0
        err = realmax('double');
        return;
    end
    irf = makeGammaShiftedIRF(t, shift, k, theta);
    model = buildDecayModelFromIRF(irf, t, tau0, true);
    coeff = lsqnonneg(model, y);
    yhat = model * coeff;
    err = sum((y - yhat).^2 ./ max(yhat, 1));
end

function err = irfFitErrExGauss(p, y, t, tau0)
    mu = p(1);
    sigma = exp(p(2));
    taui = exp(p(3));
    if ~isfinite(mu) || ~isfinite(sigma) || ~isfinite(taui) || sigma <= 0 || taui <= 0
        err = realmax('double');
        return;
    end
    irf = makeExGaussIRF(t, mu, sigma, taui);
    model = buildDecayModelFromIRF(irf, t, tau0, true);
    coeff = lsqnonneg(model, y);
    yhat = model * coeff;
    err = sum((y - yhat).^2 ./ max(yhat, 1));
end

function irf = makeGammaShiftedIRF(t, shift, k, theta)
    tt = t - shift;
    irf = zeros(size(t));
    idx = tt > 0;
    irf(idx) = (tt(idx).^(k-1) .* exp(-tt(idx) ./ theta)) ./ (gamma(k) * theta.^k);
    irf = max(irf, 0);
    irf = irf ./ max(sum(irf), eps);
end

function irf = makeExGaussIRF(t, mu, sigma, taui)
    lambda = 1 / max(taui, eps);
    z = (mu + lambda * sigma.^2 - t) ./ (sqrt(2) * sigma);
    irf = (lambda / 2) .* exp((lambda / 2) .* (2*mu + lambda*sigma.^2 - 2*t)) .* erfc(z);
    irf = max(irf, 0);
    irf = irf ./ max(sum(irf), eps);
end

function M = buildDecayModelFromIRF(irf, t, tau0, includeBG)
    n = numel(t);
    nExp = numel(tau0);
    M = zeros(n, nExp + includeBG);
    if includeBG
        M(:,1) = 1;
    end
    for k = 1:nExp
        basis = exp(-t / max(tau0(k), eps));
        convSig = conv(irf(:), basis(:), 'full');
        M(:, includeBG + k) = convSig(1:n);
    end
end
