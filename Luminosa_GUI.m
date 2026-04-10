function Luminosa_GUI
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
    app.cbImage = [];
    app.irfCache = struct('key', '', 'irf', [], 'meta', struct());
    app.txtFitResults = [];
    app.txtLegend = [];
    app.pnlTCSPCSide = [];
    app.txtVerbose = [];
    app.lblBusy = [];
    app.verboseLines = {};
    app.dividerLine = [];
    app.axResidual = [];
    app.pnlTCSPCAxes = [];
    app.roiListeners = [];
    app.roiUpdateBusy = false;
    app.roiUpdateTic = tic;
    app.lastRoiUpdateSec = -inf;
    app.defaultDataPath = 'D:\Luminosa\Data\260323\t1_20260323-103604\t1_20260323-103604';
    app.smallTitleFont = 8;
    
    % Series (video) related variables
    app.seriesData = [];        % Array of ptuOut structures for each frame
    app.seriesFolderPath = '';  % Path to series folder
    app.seriesFiles = [];       % List of series files
    app.currentFrame = 1;       % Current frame index
    app.frameSlider = [];       % Frame navigation slider
    app.seriesIntensityCache = [];  % Pre-computed intensity images for fast scrolling
    app.seriesTcspcPixCache = [];   % Pre-computed tcspc_pix for fast ROI access
    app.seriesFlimCache = [];       % Pre-computed FLIM images for fast scrolling
    app.seriesFrameFileMap = [];    % Map each frame to its source PTU file
    app.ptuOutOriginal = [];        % Original full-file data for global IRF calculation

    app.fig = uifigure('Name','Luminosa FLIM / ISM-FLIM / FCS', 'Position',[80 60 1080 680]);
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
        rowLoad = uigridlayout(ctlGrid, [1 2]);
        rowLoad.ColumnWidth = {'1x', '1x'};
        rowLoad.Padding = [0 0 0 0];
        app.btnLoad = uibutton(rowLoad, 'Text', 'Load PTU', 'ButtonPushedFcn', @onLoadPTU, 'FontSize', 12);
        app.btnLoadSeries = uibutton(rowLoad, 'Text', 'Load Series', 'ButtonPushedFcn', @onLoadSeries, 'FontSize', 12);

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
        rowROI.ColumnWidth = {70, 90, '1x'};
        rowROI.Padding = [0 0 0 0];
        app.btnSelectROI = uibutton(rowROI, 'Text', 'ROI', 'ButtonPushedFcn', @onSelectROI, 'FontSize', 11);
        app.dropROIShape = uidropdown(rowROI, 'Items', {'Rectangle','Ellipse','Lasso'}, 'Value', 'Rectangle');
        app.dropROIShape.Tooltip = 'ROI shape';
        app.btnShowTCSPC = uibutton(rowROI, 'Text', 'Show Full TCSPC', 'ButtonPushedFcn', @onShowFullTCSPC, 'FontSize', 11);

        % Row 6
        rowNsub = uigridlayout(ctlGrid, [1 2]);
        rowNsub.ColumnWidth = {120,'1x'};
        rowNsub.Padding = [0 0 0 0];
        uilabel(rowNsub, 'Text', 'Log bins / octave');
        app.txtBins = uieditfield(rowNsub, 'numeric', 'Value', 8, 'Limits', [2 32], ...
            'LowerLimitInclusive', 'on', 'UpperLimitInclusive', 'on', ...
            'ValueChangedFcn', @onTcspcBinningChanged);

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
        rowVerbose = uigridlayout(ctlGrid, [1 2]);
        rowVerbose.ColumnWidth = {55, '1x'};
        rowVerbose.Padding = [0 0 0 0];
        app.lblBusy = uilabel(rowVerbose, 'Text', 'IDLE');
        app.lblBusy.HorizontalAlignment = 'left';
        app.lblBusy.FontWeight = 'bold';
        app.txtVerbose = uitextarea(rowVerbose, 'Editable', 'off', 'Value', {'Ready.'});

        axGrid = uigridlayout(axp, [2 1]);
        axGrid.RowHeight = {'2.2x','1x'};
        axGrid.Padding = [6 6 6 6];

        % Image area with frame slider on left
        imageGrid = uigridlayout(axGrid, [1 2]);
        imageGrid.ColumnWidth = {30, '1x'};  % Slider 30px wide, image takes remaining space
        imageGrid.Padding = [0 0 0 0];
        imageGrid.ColumnSpacing = 4;
        
        % Vertical frame slider (initially hidden) - now on left
        app.frameSlider = uislider(imageGrid, 'Orientation', 'vertical', ...
            'Limits', [1 2], 'Value', 1, 'MajorTicks', [], ...
            'ValueChangedFcn', @onFrameChanged, 'Visible', 'off');
        app.frameSlider.Tooltip = 'Video frame navigation';
        
        app.axImage = uiaxes(imageGrid);
        title(app.axImage, 'Image');

        tcspcGrid = uigridlayout(axGrid, [1 2]);
        tcspcGrid.ColumnWidth = {'3.8x','1.5x'};
        tcspcGrid.ColumnSpacing = 8;
        tcspcGrid.Padding = [0 0 0 0];

        app.pnlTCSPCAxes = uipanel(tcspcGrid, 'BorderType', 'none');
        app.axTCSPC = uiaxes(app.pnlTCSPCAxes);
        app.axResidual = uiaxes(app.pnlTCSPCAxes);
        app.axTCSPC.Units = 'normalized';
        app.axResidual.Units = 'normalized';
        app.axTCSPC.Box = 'on';
        app.axResidual.Box = 'on';
        layoutTcspcAxes(false);
        title(app.axTCSPC, 'TCSPC');
        setSmallTitles();

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
        if isfolder(app.defaultDataPath)
            [f, p] = uigetfile('*.ptu', 'Select PTU file', app.defaultDataPath);
        else
            [f, p] = uigetfile('*.ptu', 'Select PTU file');
        end
        if isequal(f, 0)
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false));
        app.lastFile = fullfile(p, f);
        addStatus('Loading PTU...');
        drawnow;

        storeTCSPC = app.chkStoreTCSPC.Value;
        photonsPerChunk = max(1e5, app.editChunkM.Value * 1e6);
        useGPU = app.chkUseGPU.Value && gpuIsAvailable() && storeTCSPC;

        try
            % Clear existing data
            app.flim = [];
            app.ismRes = [];
            app.pattern = [];
            app.tcspc = [];
            app.tcspcFit = [];
            app.tcspcGlobal = [];
            app.irfCache = struct('key', '', 'irf', [], 'meta', struct());
            app.irfGlobal = [];
            app.irfGlobalModel = '';
            app.irfGlobalMeta = struct();
            app.irfGlobalDtNs = [];
            app.globalFit = [];
            app.tcspcDisplayMode = 'none';
            clearTCSPCSidePanel();
            
            % Clear series data (switch to single file mode)
            app.seriesData = {};
            app.seriesIntensityCache = {};
            app.seriesTcspcPixCache = {};
            app.seriesFlimCache = {};
            app.seriesFrameFileMap = [];
            app.currentFrameIndex = 1;

            if ~isempty(app.roi) && isvalid(app.roi)
                delete(app.roi);
            end
            app.roi = [];
            clearRoiListeners();

            % Check file size to decide on loading method
            fileInfo = dir(app.lastFile);
            fileSizeMB = fileInfo.bytes / (1024 * 1024);
            
            if fileSizeMB < 300 && exist('mHist4', 'file') == 2
                % Use mHist4 for smaller files with potential multiframe support
                addStatus('Loading PTU with mHist4 for multiframe detection...');
                drawnow;
                
                % Load PTU header to check for multiframe
                if exist('PTU_Read_Head', 'file') == 2
                    header = PTU_Read_Head(app.lastFile);
                    if isfield(header, 'TTResult_NumberOfRecords') && header.TTResult_NumberOfRecords > 1000000
                        % Potentially multiframe - use PTU_FLIM_GPU with frame extraction
                        opts = struct();
                        opts.computePerFrame = true;
                        opts.photonsPerChunk = photonsPerChunk;
                        opts.storeTcspcPix = storeTCSPC;
                        opts.storePhotonLists = storeTCSPC;  % Need photon lists for per-frame extraction
                        opts.storeTimeCell = false;
                        opts.useGPU = useGPU;
                        opts.showWaitbar = false;
                        
                        addStatus('Loading with PTU_FLIM_GPU for multiframe analysis...');
                        drawnow;
                        
                        try
                            app.ptuOut = PTU_FLIM_GPU(app.lastFile, opts);
                            addStatus(sprintf('PTU_FLIM_GPU completed. Data structure fields: %s', strjoin(fieldnames(app.ptuOut), ', ')));
                        catch ME
                            addStatus(sprintf('PTU_FLIM_GPU failed: %s', ME.message));
                            app.frameSlider.Visible = 'off';
                            loadSingleFramePTU();
                            return;
                        end
                        
                        % Store original full-file data for global IRF calculation
                        app.ptuOutOriginal = app.ptuOut;
                        
                        % Check if we have multiple frames with robust error handling
                        hasNFramesField = false;
                        nFramesValue = 1;
                        hasMultiFrameTag = false;
                        actualFrames = 1;
                        
                        try
                            if isfield(app.ptuOut, 'nFrames') && ~isempty(app.ptuOut.nFrames) && isfinite(app.ptuOut.nFrames)
                                hasNFramesField = true;
                                nFramesValue = app.ptuOut.nFrames;
                            end
                        catch ME
                            addStatus(sprintf('Warning: Could not access nFrames field: %s', ME.message));
                        end
                        
                        try
                            if isfield(app.ptuOut, 'tag') && ~isempty(app.ptuOut.tag) && ndims(app.ptuOut.tag) >= 4
                                actualFrames = size(app.ptuOut.tag, 4);
                                if actualFrames > 1
                                    hasMultiFrameTag = true;
                                end
                            end
                        catch ME
                            addStatus(sprintf('Warning: Could not analyze tag dimensions: %s', ME.message));
                        end
                        
                        % Determine if this is a multiframe file
                        isMultiFrame = (hasNFramesField && nFramesValue > 1) || (hasMultiFrameTag && actualFrames > 1);
                        
                        if isMultiFrame
                            % Use the maximum of nFrames and actual tag frames for safety
                            numFrames = max(nFramesValue, actualFrames);
                            addStatus(sprintf('Detected multiframe PTU: nFrames=%d, tag frames=%d. Using %d frames.', nFramesValue, actualFrames, numFrames));
                            drawnow;
                            
                            % Extract individual frames
                            app.seriesData = cell(numFrames, 1);
                            app.seriesIntensityCache = cell(numFrames, 1);
                            app.seriesFlimCache = cell(numFrames, 1);
                            if storeTCSPC
                                app.seriesTcspcPixCache = cell(numFrames, 1);
                            end
                            
                            for frameIdx = 1:numFrames
                                addStatus(sprintf('  -> Processing frame %d/%d...', frameIdx, numFrames));
                                
                                % Create frame-specific data structure without inheriting global cubes
                                frameData = prepareExtractedFrameStruct(app.ptuOut);
                                
                                % Extract frame-specific intensity data
                                intensityImg = [];
                                if isfield(app.ptuOut, 'tag') && size(app.ptuOut.tag, 4) >= frameIdx
                                    frameData.tag = app.ptuOut.tag(:,:,:,frameIdx);
                                    frameData.tags = frameData.tag;
                                    intensityImg = sum(frameData.tag, 3);
                                    addStatus(sprintf('    Frame %d intensity: %dx%d pixels, %d total counts', frameIdx, size(intensityImg,1), size(intensityImg,2), sum(intensityImg(:))));
                                elseif isfield(app.ptuOut, 'tags') && ~isempty(app.ptuOut.tags)
                                    % Fallback to global tags if per-frame not available
                                    frameData.tags = app.ptuOut.tags;
                                    intensityImg = frameData.tags;
                                    addStatus(sprintf('    Frame %d: Using global intensity data as fallback', frameIdx));
                                end
                                
                                % Always cache intensity image (even if empty)
                                if ~isempty(intensityImg)
                                    app.seriesIntensityCache{frameIdx} = intensityImg;
                                else
                                    % Create intensity from other available data
                                    intensityImg = getIntensityMapFromPTUData(frameData);
                                    app.seriesIntensityCache{frameIdx} = intensityImg;
                                    if ~isempty(intensityImg)
                                        addStatus(sprintf('    Frame %d: Generated intensity from available data (%dx%d)', frameIdx, size(intensityImg,1), size(intensityImg,2)));
                                    else
                                        addStatus(sprintf('    Frame %d: No intensity data available', frameIdx));
                                    end
                                end
                                
                                % Extract frame-specific tau data
                                if isfield(app.ptuOut, 'tau') && size(app.ptuOut.tau, 4) >= frameIdx
                                    frameData.tau = app.ptuOut.tau(:,:,:,frameIdx);
                                    frameData.taus = frameData.tau;  % Ensure taus field exists
                                end
                                app.seriesFlimCache{frameIdx} = getTauMeanMapFromPTUData(frameData);
                                
                                % Keep frame-specific photons for per-frame TCSPC/ROI work.
                                if storeTCSPC
                                    tcspcDataAssigned = false;
                                    if isfield(app.ptuOut, 'im_frame') && ~isempty(app.ptuOut.im_frame)
                                        framePhotonMask = (app.ptuOut.im_frame == frameIdx);
                                        if any(framePhotonMask)
                                            frameData.im_sync = app.ptuOut.im_sync(framePhotonMask);
                                            frameData.im_tcspc = app.ptuOut.im_tcspc(framePhotonMask);
                                            frameData.im_chan = app.ptuOut.im_chan(framePhotonMask);
                                            frameData.im_line = app.ptuOut.im_line(framePhotonMask);
                                            frameData.im_col = app.ptuOut.im_col(framePhotonMask);
                                            frameData.im_frame = ones(sum(framePhotonMask), 1, 'like', app.ptuOut.im_frame);
                                            tcspcDataAssigned = true;
                                            addStatus(sprintf('    Frame %d: Extracted %d photons for frame-specific TCSPC', frameIdx, numel(frameData.im_tcspc)));
                                        else
                                            addStatus(sprintf('    Frame %d: No photons found for this frame', frameIdx));
                                        end
                                    end
                                    
                                    app.seriesTcspcPixCache{frameIdx} = [];
                                    if ~tcspcDataAssigned
                                        frameData.im_sync = [];
                                        frameData.im_tcspc = uint16([]);
                                        frameData.im_chan = uint8([]);
                                        frameData.im_line = uint16([]);
                                        frameData.im_col = uint16([]);
                                        frameData.im_frame = uint16([]);
                                        addStatus(sprintf('    Frame %d: No frame-specific TCSPC data available', frameIdx));
                                    end
                                end
                                [frameData.cachedGlobalDecayCounts, frameData.cachedGlobalDecayDtNs, frameData.cachedGlobalDecaySrcInfo] = ...
                                    cachedFrameWholeDecay(frameData);
                                
                                app.seriesData{frameIdx} = frameData;
                            end
                            
                            app.seriesFrameFileMap = ones(numFrames, 1); % All frames from same file
                            
                            % Calculate global IRF from original full-file data
                            addStatus('Calculating global IRF from all frames...');
                            drawnow;
                            irfOK = ensureGlobalIRF(true);
                            if irfOK
                                addStatus('Global IRF calculated successfully.');
                            end
                            
                            % Set up frame slider
                            app.frameSlider.Limits = [1, numFrames];
                            app.frameSlider.Value = 1;
                            app.frameSlider.Visible = 'on';
                            
                            % Load first frame and calculate its FLIM
                            app.currentFrameIndex = 1;
                            loadSeriesFrame(1);
                            
                            addStatus(sprintf('PTU loaded with %d frames. Global IRF calculated. Use frame slider to navigate.', numFrames));
                        else
                            % Single frame file
                            app.frameSlider.Visible = 'off';
                            loadSingleFramePTU();
                        end
                    else
                        % Single frame file
                        app.frameSlider.Visible = 'off';
                        loadSingleFramePTU();
                    end
                else
                    % Fallback to standard loading
                    app.frameSlider.Visible = 'off';
                    loadSingleFramePTU();
                end
            else
                % Standard loading for larger files
                app.frameSlider.Visible = 'off';
                loadSingleFramePTU();
            end
            
        catch ME
            addStatus(['Load failed: ' ME.message]);
        end
        
        function loadSingleFramePTU()
            % Standard PTU loading logic
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
            if irfOK && fitOK
                addStatus(sprintf('PTU loaded. Whole-file IRF estimated with %s model; best global fit selected automatically.', app.dropIRF.Value));
            elseif irfOK
                addStatus(sprintf('PTU loaded. Whole-file IRF estimated with %s model, but automatic global fit failed.', app.dropIRF.Value));
            else
                addStatus('PTU loaded, but global IRF estimation failed.');
            end
        end
    end
    
    function onLoadSeries(~, ~)
        folderPath = uigetdir(app.defaultDataPath, 'Select folder containing Series_*.ptu files');
        if isequal(folderPath, 0)
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false));
        
        addStatus('Loading series...');
        drawnow;
        
        % Find and sort series files
        pattern = 'Series_*.ptu';
        files = dir(fullfile(folderPath, pattern));
        if isempty(files)
            addStatus(sprintf('No files matching %s found in selected folder.', pattern));
            return;
        end
        
        % Sort by numeric index in Series_xx.ptu
        idx = nan(numel(files), 1);
        for k = 1:numel(files)
            tok = regexp(files(k).name, 'Series_(\d+)\.ptu', 'tokens', 'once');
            if ~isempty(tok)
                idx(k) = str2double(tok{1});
            end
        end
        if all(isfinite(idx))
            [~, order] = sort(idx);
            files = files(order);
        else
            [~, order] = sort({files.name});
            files = files(order);
        end
        
        app.seriesFolderPath = folderPath;
        app.seriesFiles = files;
        app.currentFrame = 1;
        app.ptuOutOriginal = [];
        
        % Initialize variables for frame extraction
        addStatus(sprintf('Found %d PTU files. Extracting multiframes...', numel(files)));
        drawnow;
        
        % Pre-load ALL PTU files and extract individual frames
        storeTCSPC = app.chkStoreTCSPC.Value;
        photonsPerChunk = max(1e5, app.editChunkM.Value * 1e6);
        useGPU = app.chkUseGPU.Value && gpuIsAvailable() && storeTCSPC;
        
        allFrameData = {};
        allIntensityCache = {};
        allTcspcPixCache = {};
        allFlimCache = {};
        frameFileMap = [];  % Track which file each frame came from
        totalFrameCount = 0;
        
        for fileIdx = 1:numel(files)
            filename = fullfile(app.seriesFolderPath, files(fileIdx).name);
            addStatus(sprintf('Loading PTU file %d/%d: %s', fileIdx, numel(files), files(fileIdx).name));
            drawnow;
            
            try
                % Load PTU data with per-frame computation enabled
                if exist('PTU_FLIM_GPU', 'file') == 2
                    % Use PTU_FLIM_GPU which naturally separates frames
                    opts = struct();
                    opts.computePerFrame = true;  % Enable per-frame computation
                    opts.storeTcspcPix = storeTCSPC;
                    opts.storePhotonLists = storeTCSPC;
                    opts.storeTimeCell = false;
                    opts.useGPU = useGPU;
                    ptuData = PTU_FLIM_GPU(filename, opts);
                elseif exist('PTU_MultiFrameScanReadFast_nativeTCSPC', 'file') == 2
                    ptuData = PTU_MultiFrameScanReadFast_nativeTCSPC(filename, photonsPerChunk, storeTCSPC, useGPU, true);
                elseif exist('PTU_MultiFrameScanReadFast_multiTauTCSPC', 'file') == 2
                    ptuData = PTU_MultiFrameScanReadFast_multiTauTCSPC(filename, photonsPerChunk, storeTCSPC, useGPU, true, 8, true);
                elseif exist('PTU_MultiFrameScanReadFast', 'file') == 2
                    ptuData = PTU_MultiFrameScanReadFast(filename, photonsPerChunk, storeTCSPC, useGPU);
                else
                    error('No compatible PTU reader was found on the MATLAB path.');
                end
                
                % Extract individual frames from the multiframe PTU with robust error handling
                numFramesInFile = 1;
                hasMultipleFrames = false;
                
                % Check for nFrames field safely
                try
                    if isfield(ptuData, 'nFrames') && ~isempty(ptuData.nFrames) && isfinite(ptuData.nFrames) && ptuData.nFrames > 0
                        numFramesInFile = ptuData.nFrames;
                        hasMultipleFrames = (numFramesInFile > 1);
                    end
                catch ME
                    addStatus(sprintf('  Warning: Could not access nFrames in %s: %s', files(fileIdx).name, ME.message));
                end
                
                % Also check tag dimensions as fallback
                try
                    if isfield(ptuData, 'tag') && ~isempty(ptuData.tag) && ndims(ptuData.tag) >= 4
                        tagFrames = size(ptuData.tag, 4);
                        if tagFrames > 1
                            numFramesInFile = max(numFramesInFile, tagFrames);
                            hasMultipleFrames = true;
                        end
                    end
                catch ME
                    addStatus(sprintf('  Warning: Could not analyze tag dimensions in %s: %s', files(fileIdx).name, ME.message));
                end
                
                if hasMultipleFrames
                    addStatus(sprintf('  -> Extracting %d frames from %s', numFramesInFile, files(fileIdx).name));
                    
                    for frameIdx = 1:numFramesInFile
                        totalFrameCount = totalFrameCount + 1;
                        frameFileMap(totalFrameCount) = fileIdx;
                        
                        % Create individual frame data structure without inheriting global cubes
                        frameData = prepareExtractedFrameStruct(ptuData);
                        
                        % Extract frame-specific data
                        if isfield(ptuData, 'tag') && size(ptuData.tag, 4) >= frameIdx
                            frameData.tag = ptuData.tag(:,:,:,frameIdx);
                            frameData.tags = frameData.tag;  % Single frame intensity
                        elseif isfield(ptuData, 'tags')
                            frameData.tags = ptuData.tags;  % Use accumulated if per-frame not available
                        end
                        
                        if isfield(ptuData, 'tau') && size(ptuData.tau, 4) >= frameIdx
                            frameData.tau = ptuData.tau(:,:,:,frameIdx);
                            frameData.taus = frameData.tau;  % Single frame tau
                        elseif isfield(ptuData, 'taus')
                            frameData.taus = ptuData.taus;  % Use accumulated if per-frame not available
                        end
                        
                        % Extract frame-specific photons if available
                        if isfield(ptuData, 'im_frame') && ~isempty(ptuData.im_frame)
                            framePhotonMask = (ptuData.im_frame == frameIdx);
                            if any(framePhotonMask)
                                frameData.im_sync = ptuData.im_sync(framePhotonMask);
                                frameData.im_tcspc = ptuData.im_tcspc(framePhotonMask);
                                frameData.im_chan = ptuData.im_chan(framePhotonMask);
                                frameData.im_line = ptuData.im_line(framePhotonMask);
                                frameData.im_col = ptuData.im_col(framePhotonMask);
                                frameData.im_frame = uint16(ones(sum(framePhotonMask), 1));  % All frame 1 for this extracted frame
                            else
                                % No photons for this frame
                                frameData.im_sync = [];
                                frameData.im_tcspc = uint16([]);
                                frameData.im_chan = uint8([]);
                                frameData.im_line = uint16([]);
                                frameData.im_col = uint16([]);
                                frameData.im_frame = uint16([]);
                            end
                        end
                        [frameData.cachedGlobalDecayCounts, frameData.cachedGlobalDecayDtNs, frameData.cachedGlobalDecaySrcInfo] = ...
                            cachedFrameWholeDecay(frameData);
                        
                        % Store the extracted frame
                        allFrameData{totalFrameCount} = frameData;
                        
                        % Cache intensity image for this frame
                        intensityImg = getIntensityMapFromPTUData(frameData);
                        allIntensityCache{totalFrameCount} = intensityImg;
                        
                        allTcspcPixCache{totalFrameCount} = [];
                        allFlimCache{totalFrameCount} = getTauMeanMapFromPTUData(frameData);
                    end
                else
                    % Single frame file - treat as one frame
                    totalFrameCount = totalFrameCount + 1;
                    frameFileMap(totalFrameCount) = fileIdx;
                    allFrameData{totalFrameCount} = ptuData;
                    allIntensityCache{totalFrameCount} = getIntensityMapFromPTUData(ptuData);
                    
                    if isfield(ptuData, 'tcspc_pix') && ~isempty(ptuData.tcspc_pix)
                        allTcspcPixCache{totalFrameCount} = ptuData.tcspc_pix;
                    elseif isfield(ptuData, 'tcspc_pix_mt') && ~isempty(ptuData.tcspc_pix_mt)
                        allTcspcPixCache{totalFrameCount} = ptuData.tcspc_pix_mt;
                    end
                    
                    if ~isempty(allTcspcPixCache{totalFrameCount})
                        try
                            flimData = quickFLIMFromTCSPCFlexible(ptuData, useGPU);
                            if isfield(flimData, 'total') && isfield(flimData.total, 'tauMean')
                                allFlimCache{totalFrameCount} = flimData.total.tauMean;
                            end
                        catch
                            allFlimCache{totalFrameCount} = [];
                        end
                    else
                        allFlimCache{totalFrameCount} = [];
                    end
                end
                
            catch ME
                addStatus(['Failed to load PTU file ' num2str(fileIdx) ': ' ME.message]);
                continue;  % Skip this file and continue with others
            end
        end
        
        % Store all extracted frames
        app.seriesData = allFrameData;
        app.seriesIntensityCache = allIntensityCache;
        app.seriesTcspcPixCache = allTcspcPixCache;
        app.seriesFlimCache = allFlimCache;
        app.seriesFrameFileMap = frameFileMap;  % Track which file each frame came from
        
        % Set current frame to first frame and update display
        app.currentFrame = 1;
        if totalFrameCount > 0
            app.ptuOut = app.seriesData{1};
        end
        
        % Setup frame slider
        if totalFrameCount > 1
            app.frameSlider.Limits = [1 totalFrameCount];
            app.frameSlider.Value = 1;
            app.frameSlider.MajorTicks = 1:max(1, floor(totalFrameCount/10)):totalFrameCount;
            app.frameSlider.Visible = 'on';
        else
            app.frameSlider.Visible = 'off';
        end
        
        % Initialize display with first frame
        if totalFrameCount > 0
            showIntensityFromCache(1);
            
            % Setup global IRF and fitting for first frame
            irfOK = ensureGlobalIRF(true);
            fitOK = false;
            if irfOK
                fitOK = autoFitGlobalTCSPC();
            end
            showGlobalTCSPC();
            
            % Update title to show frame info
            if ~isempty(frameFileMap) && frameFileMap(1) <= numel(files)
                sourceFile = files(frameFileMap(1)).name;
                title(app.axImage, sprintf('Frame 1/%d from %s', totalFrameCount, sourceFile));
            else
                title(app.axImage, sprintf('Frame 1/%d', totalFrameCount));
            end
        end
        
        addStatus(sprintf('Series loaded: %d total frames extracted from %d PTU files. Use slider for instant navigation.', totalFrameCount, numel(files)));
    end
    
    function onFrameChanged(~, ~)
        if isempty(app.seriesData)
            return;
        end
        frameNum = round(app.frameSlider.Value);
        if frameNum ~= app.currentFrame
            loadSeriesFrame(frameNum);
        end
    end
    
    function loadSeriesFrame(frameNum)
        if frameNum < 1 || frameNum > numel(app.seriesData)
            return;
        end
        
        % Fast frame switching using cached data
        if frameNum <= numel(app.seriesData) && ~isempty(app.seriesData{frameNum})
            requestedDisplayMode = app.displayMode;
            
            % Set current frame data - no loading needed since it's cached
            app.currentFrame = frameNum;
            app.currentFrameIndex = frameNum;  % Keep both variables in sync
            app.ptuOut = app.seriesData{frameNum};
            app.frameSlider.Value = frameNum;
            
            % Ensure frame has proper TCSPC data structure
            if ~isfield(app.ptuOut, 'tcspc_pix') && frameNum <= numel(app.seriesTcspcPixCache) && ...
               ~isempty(app.seriesTcspcPixCache{frameNum})
                app.ptuOut.tcspc_pix = app.seriesTcspcPixCache{frameNum};
            end
            
            % Clear previous analysis results for new frame (but preserve global IRF)
            app.flim = [];
            app.ismRes = [];
            app.pattern = [];
            app.tcspc = [];
            app.tcspcFit = [];
            app.tcspcGlobal = [];
            app.globalFit = [];
            app.tcspcDisplayMode = 'none';
            clearTCSPCSidePanel();
            % Note: Keep app.irfGlobal, app.irfGlobalModel, app.irfGlobalMeta, app.irfGlobalDtNs
            % for multiframe files since IRF should be calculated from full file data

            if ~isempty(app.roi) && isvalid(app.roi)
                delete(app.roi);
            end
            app.roi = [];
            clearRoiListeners();
            
            % Show cached intensity for instant switching
            showIntensityFromCache(frameNum);
            
            % Setup IRF using original full-file data for proper fitting
            if ~isempty(app.ptuOutOriginal)
                % Temporarily switch to original data for IRF calculation
                tempPtuOut = app.ptuOut;
                app.ptuOut = app.ptuOutOriginal;
                irfOK = ensureGlobalIRF(false);  % Don't force recalculation if already done
                app.ptuOut = tempPtuOut;  % Switch back to frame data
            else
                irfOK = ensureGlobalIRF(true);
            end
            
            % Calculate and display global TCSPC for current frame
            fitOK = false;
            if irfOK
                fitOK = autoFitGlobalTCSPC();
            end
            showGlobalTCSPC();
            
            % Prefer the per-frame FLIM extracted with the frame data.
            cachedTau = [];
            if ~isempty(app.seriesFlimCache) && frameNum <= numel(app.seriesFlimCache)
                cachedTau = app.seriesFlimCache{frameNum};
            end
            if ~isempty(cachedTau)
                app.flim = flimStructFromTauMap(cachedTau, getIntensityMapFromPTUData(app.ptuOut));
                addStatus(sprintf('Using cached FLIM for frame %d.', frameNum));
            elseif hasQuickFLIMData(app.ptuOut)
                addStatus(sprintf('Computing FLIM for frame %d...', frameNum));
                drawnow;
                
                try
                    useGPU = app.chkUseGPU.Value && gpuIsAvailable();
                    flimData = quickFLIMFromTCSPCFlexible(app.ptuOut, useGPU);
                    
                    if numel(app.seriesFlimCache) < frameNum
                        app.seriesFlimCache{frameNum} = [];
                    end
                    if isfield(flimData, 'total') && isfield(flimData.total, 'tauMean')
                        app.seriesFlimCache{frameNum} = flimData.total.tauMean;
                    end
                    
                    app.flim = flimData;
                    addStatus(sprintf('FLIM computed for frame %d.', frameNum));
                catch ME
                    addStatus(sprintf('FLIM computation failed for frame %d: %s', frameNum, ME.message));
                    app.flim = [];
                end
            else
                app.flim = [];
                addStatus(sprintf('No frame-specific FLIM data available for frame %d.', frameNum));
            end
            
            % Update title to show frame info with source file
            totalFrames = numel(app.seriesData);
            if ~isempty(app.seriesFrameFileMap) && frameNum <= length(app.seriesFrameFileMap) && ...
               app.seriesFrameFileMap(frameNum) <= numel(app.seriesFiles)
                sourceFile = app.seriesFiles(app.seriesFrameFileMap(frameNum)).name;
                title(app.axImage, sprintf('Frame %d/%d from %s', frameNum, totalFrames, sourceFile));
            else
                title(app.axImage, sprintf('Frame %d/%d', frameNum, totalFrames));
            end

            if strcmp(requestedDisplayMode, 'tau') && ~isempty(app.flim)
                showTauMean();
            end
            
            if irfOK && fitOK
                addStatus(sprintf('Switched to frame %d/%d - TCSPC and IRF ready', frameNum, totalFrames));
            else
                addStatus(sprintf('Switched to frame %d/%d (instant from cache)', frameNum, totalFrames));
            end
        else
            addStatus(sprintf('Frame %d not cached. Use Load Series to pre-load all frames.', frameNum));
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

        if ~hasQuickFLIMData(app.ptuOut)
            app.flim = flimStructFromTauMap(getTauMeanMapFromPTUData(app.ptuOut), getIntensityMapFromPTUData(app.ptuOut));
            if isempty(app.flim)
                addStatus('No frame-specific FLIM data available.');
                return;
            end
            addStatus('Using cached per-frame FLIM map.');
            showTauMean();
            return;
        end

        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
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

        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
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

    function onShowFullTCSPC(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
        showGlobalTCSPC();
        addStatus('Showing whole-file TCSPC.');
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
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
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

    function onTcspcBinningChanged(~, ~)
        if isempty(app.ptuOut)
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        try
            switch app.tcspcDisplayMode
                case 'roi_fit'
                    if ~isempty(app.roi) && isvalid(app.roi)
                        onShowTCSPC();
                        onFitTCSPC();
                    else
                        showGlobalTCSPC();
                    end
                case 'roi'
                    if ~isempty(app.roi) && isvalid(app.roi)
                        onShowTCSPC();
                    else
                        showGlobalTCSPC();
                    end
                case 'global'
                    showGlobalTCSPC();
                otherwise
                    if ~isempty(app.tcspcGlobal) || ~isempty(app.ptuOut)
                        showGlobalTCSPC();
                    end
            end
            addStatus(sprintf('TCSPC display updated to %.0f log bins/octave.', app.txtBins.Value));
        catch ME
            addStatus(['Failed to update TCSPC log binning: ' ME.message]);
        end
    end

    function onSelectROI(~, ~)
        if isempty(app.axImage) || isempty(app.axImage.Children)
            addStatus('Display an image first.');
            return;
        end
        clearRoiListeners();
        if ~isempty(app.roi) && isvalid(app.roi)
            delete(app.roi);
        end
        switch app.dropROIShape.Value
            case 'Ellipse'
                app.roi = drawellipse(app.axImage, 'Color', 'y');
            case 'Lasso'
                if exist('drawfreehand', 'file') == 2
                    app.roi = drawfreehand(app.axImage, 'Color', 'y');
                else
                    app.roi = drawpolygon(app.axImage, 'Color', 'y');
                end
            otherwise
                app.roi = drawrectangle(app.axImage, 'Color', 'y');
        end
        addStatus(sprintf('ROI selected (%s).', app.dropROIShape.Value));
        attachRoiListeners(app.roi);
        updateRoiTcspcLive(false);
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
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        [countsNative, dtNsNative, srcInfo] = tcspcFromROI_native(app.ptuOut, app.roi);
        if isempty(countsNative)
            addStatus('No ROI TCSPC available. Reload with native photon TCSPC support or Store TCSPC.');
            return;
        end

        binsPerOct = app.txtBins.Value;
        [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(countsNative, dtNsNative);
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
        app.tcspc.riseIdx = riseIdx;
        app.tcspc.displayStartIdx = riseIdx;
        app.tcspc.source = srcInfo;

        layoutTcspcAxes(false);
        cla(app.axTCSPC);
        if ~isempty(app.axResidual) && isvalid(app.axResidual)
            cla(app.axResidual);
        end
        plot(app.axTCSPC, tBinNs*1e-9, cDensity, 'k.', 'MarkerSize', 10);
        hold(app.axTCSPC, 'on');

        legendLines = {'ROI TCSPC'};
        [irfDispBase, ~, ~] = getWholeFileIRFForTargetGrid(dtNsNative, numel(countsNative));
        if ~isempty(irfDispBase)
            irfDisp = applyDisplayShiftToDecay(irfDispBase(:), shiftBins);
            irfDisp = max(irfDisp, 0);
            if ~isempty(irfDisp) && max(irfDisp) > 0
                irfScale = max(countsShift) / max(irfDisp);
                irfScaled = irfDisp * irfScale;
                [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, dtNsNative, binsPerOct);
                irfDensity = clipIrfDisplay(irfDensity, 1e-1);
                plot(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
                legendLines{end+1} = 'Whole-file IRF (scaled)';
            end
        end
        hold(app.axTCSPC, 'off');
        grid(app.axTCSPC, 'on');
        app.axTCSPC.XScale = 'linear';
        xlabel(app.axTCSPC, 'Delay time (s)');
        ylabel(app.axTCSPC, 'Counts / ns');
        applyLogYGrid(app.axTCSPC);
        title(app.axTCSPC, sprintf('ROI TCSPC + whole-file IRF (%s)', app.dropIRF.Value));
        setSmallTitles();
        setTcspcXLim(tBinNs*1e-9);

        drawResidualDivider(false);
        setLegendEntries(legendLines);
        setFitSummary({sprintf('ROI source: %s', srcInfo.mode), ...
                       sprintf('Native dt = %.3f ps', 1e3*dtNsNative), ...
                       sprintf('Shift = %.3f ns', shiftNs), ...
                       'Reconvolution fit will run automatically when ROI changes.'});
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
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5 5".');
            return;
        end

        includeBG = app.chkIncludeBG.Value;
        optimizeTau = app.chkOptimizeTau.Value;

        [irfROI, irfMeta, ~] = getWholeFileIRFForTargetGrid(app.tcspc.dtNs, numel(app.tcspc.rawCountsNative));
        if isempty(irfROI)
            addStatus('Whole-file IRF unavailable for ROI TCSPC resolution.');
            return;
        end

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

        layoutTcspcAxes(true);
        cla(app.axTCSPC);
        if ~isempty(app.axResidual) && isvalid(app.axResidual)
            cla(app.axResidual);
        end
        plot(app.axTCSPC, app.tcspc.tBinNs*1e-9, app.tcspc.cDensity, 'k.', 'MarkerSize', 10);
        hold(app.axTCSPC, 'on');
        plot(app.axTCSPC, tFitNs*1e-9, fitDensity, 'r-', 'LineWidth', 1.3);

        legendLines = {'ROI TCSPC', 'Reconvolution fit'};
        irfDisp = applyDisplayShiftToDecay(irfROI(:), app.tcspc.shiftBins);
        irfDisp = max(irfDisp, 0);
        if ~isempty(irfDisp) && max(irfDisp) > 0
            irfScale = max(app.tcspc.rawCountsShift) / max(irfDisp);
            irfScaled = irfDisp * irfScale;
            [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, app.tcspc.dtNs, app.txtBins.Value);
            irfDensity = clipIrfDisplay(irfDensity, 1e-1);
            plot(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
            legendLines{end+1} = 'Whole-file IRF (scaled)';
        end
        hold(app.axTCSPC, 'off');
        grid(app.axTCSPC, 'on');
        app.axTCSPC.XScale = 'linear';
        ylabel(app.axTCSPC, 'Counts / ns');
        applyLogYGrid(app.axTCSPC);
        title(app.axTCSPC, sprintf('ROI TCSPC fit using whole-file %s IRF', app.dropIRF.Value));
        setSmallTitles();
        setTcspcXLim(app.tcspc.tBinNs*1e-9);
        if ~isempty(app.axResidual) && isvalid(app.axResidual)
            plot(app.axResidual, app.tcspc.tBinNs(1:n)*1e-9, residualDensity, 'm-', 'LineWidth', 1.0);
            hold(app.axResidual, 'on');
            yline(app.axResidual, 0, 'k-');
            hold(app.axResidual, 'off');
            app.axResidual.XScale = 'linear';
            app.axResidual.YScale = 'linear';
            xlabel(app.axResidual, 'Delay time (s)');
            ylabel(app.axResidual, 'Residual (counts/ns)');
            grid(app.axResidual, 'on');
            rmax = max(abs(residualDensity));
            if isfinite(rmax) && rmax > 0
                app.axResidual.YLim = [-rmax rmax];
            end
        end
        drawResidualDivider(true);

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
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

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
            tcspcSrc = 'ISM reassigned TCSPC';
            dtNs = app.ptuOut.head.MeasDesc_Resolution * 1e9;
        else
            [tcspc_pix, dtNs, tcspcSrc, cubeErr] = resolvePatternMatchTcspcCube();
            if isempty(tcspc_pix)
                addStatus(cubeErr);
                return;
            end
        end

        head = app.ptuOut.head;
        pulsePeriodNs = head.MeasDesc_GlobalResolution * 1e9;

        [irf, irfMeta, irfNote] = getWholeFileIRFForTargetGrid(dtNs, size(tcspc_pix, 3));
        if isempty(irf)
            addStatus('Whole-file IRF unavailable for pattern-matching TCSPC resolution.');
            return;
        end
        if ~isempty(irfNote) && ~strcmp(irfNote, 'native')
            addStatus(sprintf('Pattern matching using %s whole-file IRF.', irfNote));
        end

        addStatus(sprintf('Running global fit + pattern match with whole-file IRF using %s...', tcspcSrc));
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
        addStatus(sprintf('Pattern match complete using whole-file %s IRF (%s). GUI Tau0 was used only for fit components.', ...
            app.dropIRF.Value, tcspcSrc));
    end
    function onPointFCS(~, ~)
        if isempty(app.ptuOut)
            app.statusFCS.Text = 'Load a PTU first.';
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
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
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
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
        appendVerbose(msg);
    end

    function clearRoiListeners()
        if isempty(app.roiListeners)
            return;
        end
        try
            for k = 1:numel(app.roiListeners)
                if isvalid(app.roiListeners(k))
                    delete(app.roiListeners(k));
                end
            end
        catch
        end
        app.roiListeners = [];
    end

    function attachRoiListeners(roi)
        clearRoiListeners();
        if isempty(roi) || ~isvalid(roi)
            return;
        end
        listeners = [];
        try
            listeners = addlistener(roi, 'ROIMoved', @(~,~) updateRoiTcspcLive(false));
        catch
        end
        try
            listeners(end+1) = addlistener(roi, 'MovingROI', @(~,~) updateRoiTcspcLive(true));
        catch
        end
        app.roiListeners = listeners;
    end

    function updateRoiTcspcLive(isMoving)
        if nargin < 1
            isMoving = false;
        end
        if isempty(app.ptuOut)
            return;
        end
        if app.roiUpdateBusy
            return;
        end
        if isMoving
            tNow = toc(app.roiUpdateTic);
            if (tNow - app.lastRoiUpdateSec) < 0.3
                return;
            end
            app.lastRoiUpdateSec = tNow;
        end
        app.roiUpdateBusy = true;
        cleanupFlag = onCleanup(@() setfieldFlag(false)); %#ok<NASGU>
        try
            onShowTCSPC();
        catch
        end
        try
            onFitTCSPC();
        catch
        end
    end

    function setfieldFlag(val)
        app.roiUpdateBusy = val;
    end

    function clearTCSPCSidePanel()
        setFitSummary({'No fit results yet.'});
        setLegendEntries({'No legend yet.'});
    end

    function setBusy(flag)
        if isempty(app.lblBusy) || ~isvalid(app.lblBusy)
            return;
        end
        if flag
            app.lblBusy.Text = 'BUSY';
            app.lblBusy.FontColor = [0.85 0.2 0.2];
        else
            app.lblBusy.Text = 'IDLE';
            app.lblBusy.FontColor = [0 0 0];
        end
        drawnow limitrate;
    end

    function appendVerbose(msg)
        if isempty(app.txtVerbose) || ~isvalid(app.txtVerbose)
            return;
        end
        if nargin < 1
            msg = '';
        end
        if ischar(msg)
            msg = {msg};
        end
        if isstring(msg)
            msg = cellstr(msg);
        end
        app.verboseLines = [app.verboseLines; msg(:)];
        maxLines = 12;
        if numel(app.verboseLines) > maxLines
            app.verboseLines = app.verboseLines(end-maxLines+1:end);
        end
        app.txtVerbose.Value = app.verboseLines;
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
        appendVerbose(lines(:));
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
        appendVerbose(lines(:));
    end

    function setSmallTitles()
        try
            if ~isempty(app.axImage) && isvalid(app.axImage)
                app.axImage.Title.FontSize = app.smallTitleFont;
            end
        catch
        end
        try
            if ~isempty(app.axTCSPC) && isvalid(app.axTCSPC)
                app.axTCSPC.Title.FontSize = app.smallTitleFont;
            end
        catch
        end
    end

    function layoutTcspcAxes(showResidual)
        if nargin < 1
            showResidual = true;
        end
        if isempty(app.axTCSPC) || ~isvalid(app.axTCSPC) || isempty(app.axResidual) || ~isvalid(app.axResidual)
            return;
        end
        left   = 0.09;
        right  = 0.03;
        top    = 0.05;
        bottom = 0.08;

        gap  = 0.01;     % Reduced gap to join the plots better
        resH = 0.35;     % Increased residual height from 0.24 to 0.35

        if showResidual
            mainH = 1 - top - bottom - gap - resH;
            mainH = max(mainH, 0.30);  % Reduced minimum main height to allow more space for residuals

            app.axTCSPC.Position    = [left, bottom + resH + gap, 1 - left - right, mainH];
            app.axResidual.Position = [left, bottom,              1 - left - right, resH];

            app.axResidual.Visible = 'on';
            app.axTCSPC.XTickLabel = {};
            app.axTCSPC.XLabel.String = '';
            app.axTCSPC.Box = 'on';
            app.axResidual.Box = 'on';
            app.axResidual.Color = [1 1 1];  % Make residual background white to match main plot
        else
            mainH = 1 - top - bottom;
            app.axTCSPC.Position    = [left, bottom, 1 - left - right, mainH];
            app.axResidual.Position = [left, bottom, 1 - left - right, 0.001];
            app.axResidual.Visible  = 'off';
            app.axTCSPC.XTickLabelMode = 'auto';
            app.axTCSPC.XLabel.String = 'Delay time (s)';
        end
    end

    function setTcspcXLim(t)
        if isempty(app.axTCSPC) || ~isvalid(app.axTCSPC)
            return;
        end
        t = double(t(:));
        t = t(isfinite(t));
        if isempty(t)
            return;
        end
        left = max(0, min(t));
        right = max(t);
        if ~isfinite(left) || ~isfinite(right)
            return;
        end
        if left == right
            right = left + 1e-9;
        end
        app.axTCSPC.XLim = [left right];
        if ~isempty(app.axResidual) && isvalid(app.axResidual)
            app.axResidual.XLim = [left right];
        end
    end

    function hideRightAxisTicks()
        if isempty(app.axTCSPC) || ~isvalid(app.axTCSPC)
            return;
        end
        try
            if numel(app.axTCSPC.YAxis) >= 2
                if isprop(app.axTCSPC.YAxis(2), 'TickValues')
                    app.axTCSPC.YAxis(2).TickValues = [];
                elseif isprop(app.axTCSPC.YAxis(2), 'Ticks')
                    app.axTCSPC.YAxis(2).Ticks = [];
                end
            end
        catch
        end
    end

    function drawResidualDivider(show)
        % Remove divider line to make plots appear joined
        if ~isempty(app.dividerLine) && isvalid(app.dividerLine)
            delete(app.dividerLine);
        end
        app.dividerLine = [];
    end

    function applyLogYGrid(ax)
        if isempty(ax) || ~isvalid(ax)
            return;
        end
        try
            drawnow limitrate;
        catch
        end
        ax.YScale = 'log';
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        yLim = ax.YLim;
        if numel(yLim) < 2 || any(~isfinite(yLim)) || yLim(1) <= 0
            return;
        end
        p1 = floor(log10(yLim(1)));
        p2 = ceil(log10(yLim(2)));
        if p2 < p1
            p2 = p1;
        end
        ax.YTick = 10.^(p1:p2);
    end

    function y = clipIrfDisplay(y, floorVal)
        if nargin < 2 || isempty(floorVal)
            floorVal = 1e-1;
        end
        y = double(y);
        y(y < floorVal) = NaN;
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
        clearImageColorbars();
        cla(app.axImage);
        app.axImage.CLimMode = 'auto';
        imagesc(app.axImage, imgG);
        axis(app.axImage, 'on');
        axis(app.axImage, 'image');
        app.axImage.Visible = 'on';
        colormap(app.axImage, gray);
        app.cbImage = colorbar(app.axImage);
        title(app.axImage, 'Intensity (gamma)');
        setSmallTitles();
        app.displayMode = 'intensity';
        drawnow;
    end

    function showIntensityFromCache(frameNum)
        % Show intensity image from pre-computed cache for fast navigation
        addStatus(sprintf('Showing intensity for frame %d (cache size: %d)', frameNum, numel(app.seriesIntensityCache)));
        
        % Check if cache exists for this frame
        if frameNum > numel(app.seriesIntensityCache)
            addStatus(sprintf('Frame %d exceeds cache size (%d). Using fallback.', frameNum, numel(app.seriesIntensityCache)));
            showIntensityFromPTU();
            return;
        end
        
        % Get cached intensity image
        img = app.seriesIntensityCache{frameNum};
        if isempty(img)
            addStatus(sprintf('Frame %d intensity cache is empty. Trying to generate from current data.', frameNum));
            % Try to generate intensity from current frame data
            if ~isempty(app.ptuOut)
                img = getIntensityMapFromPTUData(app.ptuOut);
                if ~isempty(img)
                    % Update cache with generated intensity
                    app.seriesIntensityCache{frameNum} = img;
                    addStatus(sprintf('Generated and cached intensity for frame %d (%dx%d)', frameNum, size(img,1), size(img,2)));
                else
                    showIntensityFromPTU();
                    return;
                end
            else
                showIntensityFromPTU();
                return;
            end
        end
        
        % Display the intensity image
        imgG = applyGamma(double(img), app.gamma);
        clearImageColorbars();
        cla(app.axImage);
        app.axImage.CLimMode = 'auto';
        imagesc(app.axImage, imgG);
        axis(app.axImage, 'on');
        axis(app.axImage, 'image');
        app.axImage.Visible = 'on';
        colormap(app.axImage, gray);
        app.cbImage = colorbar(app.axImage);
        title(app.axImage, sprintf('Frame %d Intensity (cached)', frameNum));
        setSmallTitles();
        app.displayMode = 'intensity';
        drawnow;
        
        addStatus(sprintf('Successfully displayed intensity for frame %d', frameNum));
    end

    function img = getIntensityMapFromPTUData(ptuData)
        % Extract intensity map from a specific PTU data structure
        img = [];
        
        if isfield(ptuData, 'tags') && ~isempty(ptuData.tags)
            img = sum(double(ptuData.tags), 3);
            return;
        end
        
        if isfield(ptuData, 'tag') && ~isempty(ptuData.tag)
            if ndims(ptuData.tag) >= 4
                img = squeeze(sum(sum(double(ptuData.tag), 3), 4));
            elseif ndims(ptuData.tag) == 3
                img = sum(double(ptuData.tag), 3);
            else
                img = double(ptuData.tag);
            end
            return;
        end
        
        if isfield(ptuData, 'tcspc_pix') && ~isempty(ptuData.tcspc_pix)
            tmp = sum(double(ptuData.tcspc_pix), 3);
            if ndims(tmp) == 3
                img = sum(tmp, 3);
            else
                img = tmp;
            end
            return;
        end
        
        if isfield(ptuData, 'tcspc_pix_mt') && ~isempty(ptuData.tcspc_pix_mt)
            tmp = sum(double(ptuData.tcspc_pix_mt), 3);
            if ndims(tmp) == 3
                img = sum(tmp, 3);
            else
                img = tmp;
            end
            return;
        end
    end

    function frameData = prepareExtractedFrameStruct(ptuData)
        frameData = struct();
        metaFields = {'head', 'dind', 'Ngate', 'Resolution_ns', 'options', ...
            'tcspc_mt_centers_ns', 'tcspc_mt_width_ns'};
        for ii = 1:numel(metaFields)
            if isfield(ptuData, metaFields{ii})
                frameData.(metaFields{ii}) = ptuData.(metaFields{ii});
            end
        end
        frameData.nFrames = 1;
    end

    function tauMean = getTauMeanMapFromPTUData(ptuData)
        tauMean = [];
        tauCube = [];
        tagCube = [];

        if isfield(ptuData, 'taus') && ~isempty(ptuData.taus)
            tauCube = double(ptuData.taus);
        elseif isfield(ptuData, 'tau') && ~isempty(ptuData.tau)
            tauCube = double(ptuData.tau);
        end

        if isempty(tauCube)
            return;
        end

        if ndims(tauCube) >= 4 && size(tauCube, 4) == 1
            tauCube = tauCube(:,:,:,1);
        end

        if isfield(ptuData, 'tags') && ~isempty(ptuData.tags)
            tagCube = double(ptuData.tags);
        elseif isfield(ptuData, 'tag') && ~isempty(ptuData.tag)
            tagCube = double(ptuData.tag);
        end

        if ndims(tagCube) >= 4 && size(tagCube, 4) == 1
            tagCube = tagCube(:,:,:,1);
        end

        if ndims(tauCube) <= 2
            tauMean = double(tauCube);
            return;
        end

        if ~isempty(tagCube) && size(tagCube, 1) == size(tauCube, 1) && ...
                size(tagCube, 2) == size(tauCube, 2) && size(tagCube, 3) == size(tauCube, 3)
            denom = sum(tagCube, 3);
            tauMean = sum(tauCube .* tagCube, 3) ./ max(denom, 1);
            zeroMask = denom <= 0;
            if any(zeroMask(:))
                tauFallback = mean(tauCube, 3, 'omitnan');
                tauMean(zeroMask) = tauFallback(zeroMask);
            end
        else
            tauMean = mean(tauCube, 3, 'omitnan');
        end
    end

    function flim = flimStructFromTauMap(tauMean, intensityImg)
        flim = [];
        if isempty(tauMean)
            return;
        end

        flim = struct();
        flim.total = struct();
        flim.total.tauMean = double(tauMean);
        flim.total.meanArrival = double(tauMean);
        flim.total.tauRMS = zeros(size(tauMean), 'double');
        flim.total.globalDecay = [];
        flim.total.tAxisNs = [];
        flim.total.t0Bin = [];
        if nargin >= 2 && ~isempty(intensityImg)
            flim.total.tag = double(intensityImg);
        else
            flim.total.tag = [];
        end
    end

    function tf = hasQuickFLIMData(ptuData)
        tf = false;
        if isempty(ptuData)
            return;
        end

        tf = (isfield(ptuData, 'tcspc_pix') && ~isempty(ptuData.tcspc_pix)) || ...
             (isfield(ptuData, 'tcspc_pix_mt') && ~isempty(ptuData.tcspc_pix_mt)) || ...
             (isfield(ptuData, 'im_tcspc') && ~isempty(ptuData.im_tcspc) && ...
              isfield(ptuData, 'im_col') && ~isempty(ptuData.im_col) && ...
              isfield(ptuData, 'im_line') && ~isempty(ptuData.im_line) && ...
              isfield(ptuData, 'im_chan') && ~isempty(ptuData.im_chan));
    end

    function [tcspcPix, dtNs, srcLabel, errMsg] = resolvePatternMatchTcspcCube()
        tcspcPix = [];
        dtNs = [];
        srcLabel = '';
        errMsg = '';

        if isempty(app.ptuOut)
            errMsg = 'Load a PTU first.';
            return;
        end

        if isfield(app.ptuOut, 'tcspc_pix') && ~isempty(app.ptuOut.tcspc_pix)
            tcspcPix = app.ptuOut.tcspc_pix;
            dtNs = app.ptuOut.head.MeasDesc_Resolution * 1e9;
            srcLabel = 'stored tcspc_pix';
            return;
        end

        if iscell(app.seriesTcspcPixCache) && ~isempty(app.currentFrame) && app.currentFrame >= 1 && ...
                app.currentFrame <= numel(app.seriesTcspcPixCache) && ...
                ~isempty(app.seriesTcspcPixCache{app.currentFrame})
            tcspcPix = app.seriesTcspcPixCache{app.currentFrame};
            app.ptuOut.tcspc_pix = tcspcPix;
            dtNs = app.ptuOut.head.MeasDesc_Resolution * 1e9;
            srcLabel = sprintf('cached frame %d tcspc_pix', app.currentFrame);
            return;
        end

        if isfield(app.ptuOut, 'im_tcspc') && ~isempty(app.ptuOut.im_tcspc) && ...
                isfield(app.ptuOut, 'im_col') && ~isempty(app.ptuOut.im_col) && ...
                isfield(app.ptuOut, 'im_line') && ~isempty(app.ptuOut.im_line) && ...
                isfield(app.ptuOut, 'im_chan') && ~isempty(app.ptuOut.im_chan)
            addStatus('Reconstructing linear tcspc_pix from photon lists for pattern matching...');
            drawnow;
            try
                [tcspcPix, dtNs] = buildLinearTcspcCubeFromPhotonLists(app.ptuOut);
                app.ptuOut.tcspc_pix = tcspcPix;
                srcLabel = 'photon-list reconstructed tcspc_pix';
                return;
            catch ME
                errMsg = ['Pattern-match TCSPC reconstruction failed: ' ME.message];
                return;
            end
        end

        if isfield(app.ptuOut, 'tcspc_pix_mt') && ~isempty(app.ptuOut.tcspc_pix_mt)
            errMsg = 'Pattern matching requires linear tcspc_pix. Current data only has tcspc_pix_mt; reload with "Store TCSPC".';
            return;
        end

        errMsg = 'No linear tcspc_pix is available. Reload with "Store TCSPC".';
    end

    function [counts, dtNs, srcInfo] = cachedFrameWholeDecay(ptuData)
        counts = [];
        dtNs = [];
        srcInfo = struct('mode', 'unknown', 'native', false, 'binned', false, 'tAxisNs', [], 'widthNs', []);
        try
            [counts, dtNs, srcInfo] = tcspcFromWholeFile_native(ptuData);
        catch
            counts = [];
            dtNs = [];
        end
    end

    function lbl = currentTcspcLabel()
        if ~isempty(app.seriesData) && numel(app.seriesData) > 1 && ...
                ~isempty(app.currentFrame) && app.currentFrame >= 1 && app.currentFrame <= numel(app.seriesData)
            lbl = sprintf('Frame %d', app.currentFrame);
        else
            lbl = 'Whole-file';
        end
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
            cmap = jet(256);
            cmap = cmap(30:end-30,:);
            trange = getTauRange(img);
            clearImageColorbars();
            cla(app.axImage);
            imagesc(app.axImage, img);
            axis(app.axImage, 'on');
            axis(app.axImage, 'image');
            colormap(app.axImage, cmap);
            caxis(app.axImage, trange);
            app.cbTau = colorbar(app.axImage);
            ylabel(app.cbTau, 'Tau (ns)');
            title(app.axImage, 'Tau mean (ns)');
            setSmallTitles();
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
            [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(countsNative, dtNs);
            [tBinNs, cDensity, cInt, edgesNs, widthNs] = logBinTCSPC(countsShift, dtNs, binsPerOct);
            app.tcspcGlobal = struct('tBinNs', tBinNs, 'cDensity', cDensity, 'cBin', cInt, ...
                'edgesNs', edgesNs, 'widthNs', widthNs, 'dtNs', dtNs, 'rawCountsNative', countsNative(:), ...
                'rawCountsShift', countsShift(:), 'shiftBins', shiftBins, 'shiftNs', shiftNs, ...
                'peakIdx', peakIdx, 'riseIdx', riseIdx, 'displayStartIdx', riseIdx, 'source', srcInfo);
        end

        residualDensity = [];
        fitDensityPlot = [];
        if ~isempty(app.globalFit)
            app.globalFit = refreshGlobalFitDisplayCache(app.globalFit, app.tcspcGlobal, app.txtBins.Value);
        end
        if ~isempty(app.globalFit) && isfield(app.globalFit, 'fitDensity') && ~isempty(app.globalFit.fitDensity)
            fitDensityPlot = app.globalFit.fitDensity(:);
            n = min(numel(app.tcspcGlobal.cDensity), numel(fitDensityPlot));
            residualDensity = app.tcspcGlobal.cDensity(1:n) - fitDensityPlot(1:n);
            summaryLines = app.globalFit.summaryLines;
            if isfield(app.globalFit, 'displayChi2red') && ~isempty(app.globalFit.displayChi2red) && isfinite(app.globalFit.displayChi2red)
                summaryLines = [{sprintf('Display chi2/ndf = %.4g', app.globalFit.displayChi2red)}, summaryLines];
            end
            summaryLines{end+1} = sprintf('Display bins/octave = %.0f', app.txtBins.Value);
            setFitSummary(summaryLines);
        else
            setFitSummary({sprintf('%s source: %s', currentTcspcLabel(), srcInfo.mode), ...
                           sprintf('IRF model: %s', app.dropIRF.Value), ...
                           sprintf('Display bins/octave = %.0f', app.txtBins.Value), ...
                           'No automatic fit available.'});
        end

        showResidual = ~isempty(residualDensity);
        layoutTcspcAxes(showResidual);
        cla(app.axTCSPC);
        if ~isempty(app.axResidual) && isvalid(app.axResidual)
            cla(app.axResidual);
        end
        plot(app.axTCSPC, app.tcspcGlobal.tBinNs*1e-9, app.tcspcGlobal.cDensity, 'k.', 'MarkerSize', 10);
        hold(app.axTCSPC, 'on');
        legendLines = {'Whole-file TCSPC'};

        [irfDispBase, ~, ~] = getWholeFileIRFForTargetGrid(app.tcspcGlobal.dtNs, numel(app.tcspcGlobal.rawCountsNative));
        if ~isempty(irfDispBase)
            irfDisp = applyDisplayShiftToDecay(irfDispBase(:), app.tcspcGlobal.shiftBins);
            irfDisp = max(irfDisp, 0);
            if ~isempty(irfDisp) && max(irfDisp) > 0
                irfScale = max(app.tcspcGlobal.cBin) / max(irfDisp);
                irfScaled = irfDisp * irfScale;
                [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, app.tcspcGlobal.dtNs, app.txtBins.Value);
                irfDensity = clipIrfDisplay(irfDensity, 1e-1);
                plot(app.axTCSPC, tIRFNs*1e-9, irfDensity, 'b--', 'LineWidth', 1.1);
                legendLines{end+1} = 'Whole-file IRF (scaled)';
            end
        end

        if ~isempty(fitDensityPlot)
            plot(app.axTCSPC, app.globalFit.tBinNs*1e-9, fitDensityPlot, 'r-', 'LineWidth', 1.3);
            legendLines{end+1} = sprintf('Best fit (%d-exp)', app.globalFit.nExp);
        end
        hold(app.axTCSPC, 'off');

        grid(app.axTCSPC, 'on');
        app.axTCSPC.XScale = 'linear';
        ylabel(app.axTCSPC, 'Counts / ns');
        applyLogYGrid(app.axTCSPC);
        title(app.axTCSPC, sprintf('%s TCSPC + whole-file IRF (%s)', currentTcspcLabel(), app.dropIRF.Value));
        setSmallTitles();
        setTcspcXLim(app.tcspcGlobal.tBinNs*1e-9);

        if showResidual && ~isempty(app.axResidual) && isvalid(app.axResidual)
            plot(app.axResidual, app.tcspcGlobal.tBinNs(1:numel(residualDensity))*1e-9, residualDensity, 'm-', 'LineWidth', 1.0);
            hold(app.axResidual, 'on');
            yline(app.axResidual, 0, 'k-');
            hold(app.axResidual, 'off');
            app.axResidual.XScale = 'linear';
            app.axResidual.YScale = 'linear';
            xlabel(app.axResidual, 'Delay time (s)');
            ylabel(app.axResidual, 'Residual (counts/ns)');
            grid(app.axResidual, 'on');
            rmax = max(abs(residualDensity));
            if isfinite(rmax) && rmax > 0
                app.axResidual.YLim = [-rmax rmax];
            end
            legendLines{end+1} = 'Residual';
            drawResidualDivider(true);
        else
            xlabel(app.axTCSPC, 'Delay time (s)');
            drawResidualDivider(false);
        end

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

    function [irfUse, irfMetaUse, note] = getWholeFileIRFForTargetGrid(targetDtNs, targetNBins)
        irfUse = [];
        irfMetaUse = struct();
        note = '';

        if nargin < 2 || isempty(targetNBins)
            targetNBins = [];
        end
        if isempty(targetDtNs) || ~isfinite(targetDtNs) || targetDtNs <= 0
            return;
        end
        if ~ensureGlobalIRF(false) || isempty(app.irfGlobal)
            return;
        end

        if isempty(targetNBins)
            targetNBins = numel(app.irfGlobal);
        else
            targetNBins = max(1, round(double(targetNBins)));
        end

        irfMetaUse = app.irfGlobalMeta;
        sameDt = isempty(app.irfGlobalDtNs) || abs(app.irfGlobalDtNs - targetDtNs) < 1e-6;
        if sameDt && numel(app.irfGlobal) == targetNBins
            irfUse = app.irfGlobal(:);
            note = 'native';
            return;
        end

        [irfUse, okRebin] = rebinIRFToTargetGrid(app.irfGlobal(:), app.irfGlobalDtNs, targetDtNs, targetNBins);
        if ~okRebin
            irfUse = [];
            irfMetaUse = struct();
            note = '';
            return;
        end

        if isstruct(irfMetaUse)
            irfMetaUse.rebinnedFromDtNs = app.irfGlobalDtNs;
            irfMetaUse.rebinnedToDtNs = targetDtNs;
            irfMetaUse.rebinnedToNBins = targetNBins;
        end
        note = sprintf('rebinned from %.3f ps to %.3f ps', 1e3*app.irfGlobalDtNs, 1e3*targetDtNs);
    end

    function [irfDst, ok] = rebinIRFToTargetGrid(irfSrc, dtSrcNs, dtDstNs, nDst)
        irfDst = [];
        ok = false;

        irfSrc = max(double(irfSrc(:)), 0);
        if isempty(irfSrc) || ~isfinite(dtSrcNs) || ~isfinite(dtDstNs) || dtSrcNs <= 0 || dtDstNs <= 0 || ...
                isempty(nDst) || ~isfinite(nDst) || nDst < 1
            return;
        end

        nDst = max(1, round(double(nDst)));
        srcEdges = (0:numel(irfSrc))' * dtSrcNs;
        srcCum = [0; cumsum(irfSrc)];
        dstEdges = (0:nDst)' * dtDstNs;

        dstCum = zeros(size(dstEdges));
        inRange = dstEdges > srcEdges(1) & dstEdges < srcEdges(end);
        if any(inRange)
            dstCum(inRange) = interp1(srcEdges, srcCum, dstEdges(inRange), 'linear');
        end
        dstCum(dstEdges >= srcEdges(end)) = srcCum(end);

        irfDst = diff(dstCum);
        irfDst = max(irfDst(:), 0);
        s = sum(irfDst);
        if ~(isfinite(s) && s > 0)
            irfDst = [];
            return;
        end
        irfDst = irfDst ./ s;
        ok = true;
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
        [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(countsNative, dtNs); %#ok<ASGLU>
        tcspcGlobalTmp = struct('dtNs', dtNs, 'rawCountsNative', countsNative(:), 'rawCountsShift', countsShift(:), ...
            'shiftBins', shiftBins, 'shiftNs', shiftNs, 'peakIdx', peakIdx, 'riseIdx', riseIdx, 'cBin', [], 'cDensity', []);
        bestFit = refreshGlobalFitDisplayCache(bestFit, tcspcGlobalTmp, app.txtBins.Value);
        app.globalFit = bestFit;
        ok = true;
    end

    function fitOut = refreshGlobalFitDisplayCache(fitIn, tcspcGlobalData, binsPerOct)
        fitOut = fitIn;
        if isempty(fitIn) || ~isfield(fitIn, 'fitCountsRawFull') || isempty(fitIn.fitCountsRawFull) || ...
                isempty(tcspcGlobalData) || ~isfield(tcspcGlobalData, 'dtNs') || isempty(tcspcGlobalData.dtNs) || ...
                ~isfield(tcspcGlobalData, 'shiftBins') || isempty(tcspcGlobalData.shiftBins)
            return;
        end

        fitCountsShift = applyDisplayShiftToDecay(fitIn.fitCountsRawFull, tcspcGlobalData.shiftBins);
        [tFitNs, fitDensity, fitCountsLog, fitEdgesNs] = logBinTCSPC(fitCountsShift, tcspcGlobalData.dtNs, binsPerOct);
        fitOut.fitCountsShift = fitCountsShift(:);
        fitOut.tBinNs = tFitNs(:);
        fitOut.fitDensity = fitDensity(:);
        fitOut.fitCountsLog = fitCountsLog(:);
        fitOut.fitEdgesNs = fitEdgesNs(:);
        fitOut.displayBinsPerOct = binsPerOct;
        fitOut.displayChi2red = NaN;

        if isfield(tcspcGlobalData, 'cBin') && ~isempty(tcspcGlobalData.cBin) && isfield(fitOut, 'coeff')
            n = min(numel(tcspcGlobalData.cBin), numel(fitCountsLog));
            if n > 0
                nParams = numel(fitOut.coeff);
                chi2 = sum((double(tcspcGlobalData.cBin(1:n)) - double(fitCountsLog(1:n))).^2 ./ max(double(fitCountsLog(1:n)), 1));
                ndf = max(n - nParams, 1);
                fitOut.displayChi2red = chi2 / ndf;
            end
        end
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
        clearImageColorbars();
        cla(app.axImage);
        image(app.axImage, rgb);
        axis(app.axImage, 'image');
        axis(app.axImage, 'off');
        title(app.axImage, titleStr);
        setSmallTitles();
        colormap(app.axImage, cmap);
        caxis(app.axImage, trange);
        app.cbTau = colorbar(app.axImage);
        ylabel(app.cbTau, 'Tau (ns)');
    end

    function clearImageColorbars()
        if ~isempty(app.cbTau) && isvalid(app.cbTau)
            delete(app.cbTau);
        end
        app.cbTau = [];
        if ~isempty(app.cbImage) && isvalid(app.cbImage)
            delete(app.cbImage);
        end
        app.cbImage = [];
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
    elseif isfield(ptuOut, 'im_tcspc') && ~isempty(ptuOut.im_tcspc)
        [cube, dt_ns] = buildLinearTcspcCubeFromPhotonLists(ptuOut);
        tAxisNs = ((1:size(cube, 3)) - 0.5) * dt_ns;
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
end

function [cube, dt_ns] = buildLinearTcspcCubeFromPhotonLists(ptuOut)
    if ~isfield(ptuOut, 'im_col') || ~isfield(ptuOut, 'im_line') || ~isfield(ptuOut, 'im_chan') || ...
            ~isfield(ptuOut, 'im_tcspc') || isempty(ptuOut.im_tcspc)
        error('Per-photon TCSPC data is missing.');
    end

    nx = double(ptuOut.head.ImgHdr_PixX);
    ny = double(ptuOut.head.ImgHdr_PixY);
    nPix = nx * ny;
    dt_ns = ptuOut.head.MeasDesc_Resolution * 1e9;

    if isfield(ptuOut, 'dind') && ~isempty(ptuOut.dind)
        dind = double(ptuOut.dind(:));
    else
        dind = unique(double(ptuOut.im_chan(:)));
    end
    nCh = numel(dind);

    if isfield(ptuOut, 'Ngate') && ~isempty(ptuOut.Ngate)
        Ngate = double(ptuOut.Ngate);
    else
        Ngate = double(max(ptuOut.im_tcspc(:)));
    end
    Ngate = max(1, Ngate);

    chMapLUT = zeros(1, 256, 'double');
    chMapLUT(double(dind) + 1) = 1:nCh;

    x = double(ptuOut.im_col(:));
    y = double(ptuOut.im_line(:));
    g = double(ptuOut.im_tcspc(:));
    chLocal = chMapLUT(double(ptuOut.im_chan(:)) + 1);
    chLocal = chLocal(:);

    x = x(:);
    y = y(:);
    g = g(:);

    valid = (chLocal >= 1) & (chLocal <= nCh) & ...
            (x >= 1) & (x <= nx) & ...
            (y >= 1) & (y <= ny) & ...
            (g >= 1) & (g <= Ngate);

    x = x(valid);
    y = y(valid);
    g = g(valid);
    chLocal = chLocal(valid);

    grp3 = x + (y - 1) * nx + (g - 1) * nPix + (chLocal - 1) * (nPix * Ngate);
    counts = accumarray(grp3, 1, [nPix * Ngate * nCh, 1], @sum, 0);
    cube = reshape(counts, [nx, ny, Ngate, nCh]);
end

function [counts, dtNs, srcInfo] = tcspcFromROI_native(ptuOut, roi)
% Prefer per-photon native TCSPC when available. Fall back to tcspc_pix.
    srcInfo = struct('mode', 'unknown', 'native', false);

    if isfield(ptuOut, 'im_col') && isfield(ptuOut, 'im_line')
        idx = roiPhotonMask(ptuOut, roi);
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
        [maskSub, x1, x2, y1, y2] = roiPixelMaskSub(roi, size(cube,1), size(cube,2));
        if isempty(maskSub)
            y1 = 1; y2 = size(cube,1);
            x1 = 1; x2 = size(cube,2);
            maskSub = true(y2-y1+1, x2-x1+1);
        end
        sub = double(cube(y1:y2, x1:x2, :, :));
        maskSub = double(maskSub);
        counts = squeeze(sum(sum(sum(sub .* reshape(maskSub, [size(maskSub,1), size(maskSub,2), 1, 1]), 4), 2), 1));
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

    if isfield(ptuOut, 'cachedGlobalDecayCounts') && ~isempty(ptuOut.cachedGlobalDecayCounts)
        counts = double(ptuOut.cachedGlobalDecayCounts(:));
        if isfield(ptuOut, 'cachedGlobalDecayDtNs') && ~isempty(ptuOut.cachedGlobalDecayDtNs)
            dtNs = ptuOut.cachedGlobalDecayDtNs;
        else
            dtNs = [];
        end
        if isfield(ptuOut, 'cachedGlobalDecaySrcInfo') && ~isempty(ptuOut.cachedGlobalDecaySrcInfo)
            srcInfo = ptuOut.cachedGlobalDecaySrcInfo;
        else
            srcInfo.mode = 'cached frame decay';
        end
        return;
    end

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

function [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(counts, dtNs)
% Circular-shift the decay for display while preserving the full TCSPC.
% The rise onset is moved left so no more than ~0.5 ns precedes it.
    counts = double(counts(:));
    if isempty(counts) || isempty(dtNs) || ~isfinite(dtNs) || dtNs <= 0
        countsShift = counts;
        shiftBins = 0;
        shiftNs = 0;
        peakIdx = 1;
        riseIdx = 1;
        return;
    end

    [~, peakIdx] = max(counts);
    peakIdx = max(1, peakIdx);
    riseIdx = estimateDecayRiseStartIdx(counts, dtNs, peakIdx);

    preRiseNs = 0.50;
    preRiseBins = max(0, round(preRiseNs / max(dtNs, eps)));
    targetRiseIdx = min(numel(counts), preRiseBins + 1);
    shiftBins = max(0, riseIdx - targetRiseIdx);
    shiftNs = shiftBins * dtNs;
    countsShift = applyDisplayShiftToDecay(counts, shiftBins);
end

function riseIdx = estimateDecayRiseStartIdx(counts, dtNs, peakIdx)
    y = max(double(counts(:)), 0);
    n = numel(y);
    if n < 2
        riseIdx = 1;
        return;
    end

    peakIdx = max(1, min(n, round(peakIdx)));
    if peakIdx <= 1
        riseIdx = 1;
        return;
    end

    smoothBins = max(1, min(n, ceil(0.08 / max(dtNs, eps))));
    if smoothBins > 1 && mod(smoothBins, 2) == 0
        smoothBins = smoothBins - 1;
    end
    if smoothBins > 1
        kernel = ones(smoothBins, 1) / smoothBins;
        ySmooth = conv(y, kernel, 'same');
    else
        ySmooth = y;
    end

    ySorted = sort(ySmooth);
    nBase = max(1, min(numel(ySorted), round(0.25 * numel(ySorted))));
    baselineSlice = ySorted(1:nBase);
    baseline = mean(baselineSlice);
    noise = std(baselineSlice);
    amp = max(ySmooth(peakIdx) - baseline, 0);
    if amp <= 0
        riseIdx = peakIdx;
        return;
    end

    runBins = max(2, min([12, peakIdx, ceil(0.08 / max(dtNs, eps))]));
    riseIdx = sustainedThresholdCrossing(ySmooth, peakIdx, baseline + max(4 * noise, 0.02 * amp), runBins);
    if isempty(riseIdx)
        riseIdx = sustainedThresholdCrossing(ySmooth, peakIdx, baseline + max(2 * noise, 0.01 * amp), runBins);
    end
    if isempty(riseIdx)
        riseIdx = find(ySmooth(1:peakIdx) >= (baseline + 0.01 * amp), 1, 'first');
    end
    if isempty(riseIdx)
        riseIdx = peakIdx;
    end
    riseIdx = max(1, min(peakIdx, riseIdx));
end

function idx = sustainedThresholdCrossing(y, peakIdx, threshold, runBins)
    idx = [];
    if peakIdx < 1
        return;
    end

    above = y(1:peakIdx) >= threshold;
    if isempty(above)
        return;
    end

    runBins = max(1, min(numel(above), round(runBins)));
    if runBins == 1
        idx = find(above, 1, 'first');
        return;
    end

    runHits = conv(double(above), ones(runBins, 1), 'valid');
    idx = find(runHits >= runBins, 1, 'first');
end

function [tBinNs, cDensity, cInt, edgesNs, widthNs] = logBinTCSPC(counts, dtNs, binsPerOct)
% Mixed linear/log binning for decay display.
% The first 0.5 ns remain at native resolution; later bins expand logarithmically.
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

    linearNs = 0.50;
    nLinear = min(Nt, max(1, ceil(linearNs / max(dtNs, eps))));
    linearEdgesIdx = (1:(nLinear + 1)).';

    if nLinear < Nt
        remNt = Nt - nLinear;
        nOct = max(1, log2(remNt));
        nBinsLog = max(16, ceil(nOct * binsPerOct));
        logEdgesRel = unique(round(logspace(0, log10(remNt + 1), nBinsLog + 1)));
        logEdgesRel(1) = 1;
        logEdgesRel(end) = remNt + 1;
        logEdgesIdx = nLinear + logEdgesRel(:);
        edgesIdx = unique([linearEdgesIdx; logEdgesIdx(2:end)]);
    else
        edgesIdx = linearEdgesIdx;
    end

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
        if widthNs(i) <= (1 + 1e-12) * dtNs || t0 < linearNs
            tBinNs(i) = 0.5 * (t0 + t1);
        elseif t0 <= 0
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
    fitCountsShift = applyDisplayShiftToDecay(fitCountsRawFull, shiftBins);
end

function countsShift = applyDisplayShiftToDecay(counts, shiftBins)
    counts = double(counts(:));
    if isempty(counts)
        countsShift = counts;
        return;
    end
    shiftBins = max(0, round(shiftBins));
    if shiftBins == 0
        countsShift = counts;
        return;
    end
    countsShift = circshift(counts, -shiftBins);
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
    col = double(ptuOut.im_col(:));
    row = double(ptuOut.im_line(:));
    idx = false(size(col));

    if isempty(roi) || ~isvalid(roi)
        return;
    end

    if isprop(roi, 'Center') && isprop(roi, 'SemiAxes')
        c = roi.Center;
        a = roi.SemiAxes(1);
        b = roi.SemiAxes(2);
        if a <= 0 || b <= 0
            return;
        end
        theta = 0;
        if isprop(roi, 'RotationAngle')
            theta = roi.RotationAngle * (pi / 180);
        end
        x = col - c(1);
        y = row - c(2);
        xr =  cos(theta) * x + sin(theta) * y;
        yr = -sin(theta) * x + cos(theta) * y;
        idx = (xr ./ a).^2 + (yr ./ b).^2 <= 1;
        return;
    end

    if isprop(roi, 'Position')
        pos = roi.Position;
        if size(pos,2) == 2 && size(pos,1) >= 3
            idx = inpolygon(col, row, pos(:,1), pos(:,2));
            return;
        end
        if numel(pos) >= 4
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = ceil(pos(1) + pos(3));
            y2 = ceil(pos(2) + pos(4));
            idx = col >= x1 & col <= x2 & row >= y1 & row <= y2;
            return;
        end
    end
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

function [maskSub, x1, x2, y1, y2] = roiPixelMaskSub(roi, nRows, nCols)
    maskSub = [];
    x1 = 1; x2 = nCols; y1 = 1; y2 = nRows;
    if isempty(roi) || ~isvalid(roi)
        return;
    end

    if isprop(roi, 'Center') && isprop(roi, 'SemiAxes')
        c = roi.Center;
        a = roi.SemiAxes(1);
        b = roi.SemiAxes(2);
        if a <= 0 || b <= 0
            return;
        end
        theta = 0;
        if isprop(roi, 'RotationAngle')
            theta = roi.RotationAngle * (pi / 180);
        end
        x1 = max(1, floor(c(1) - a));
        x2 = min(nCols, ceil(c(1) + a));
        y1 = max(1, floor(c(2) - b));
        y2 = min(nRows, ceil(c(2) + b));
        [yg, xg] = ndgrid(y1:y2, x1:x2);
        x = xg - c(1);
        y = yg - c(2);
        xr =  cos(theta) * x + sin(theta) * y;
        yr = -sin(theta) * x + cos(theta) * y;
        maskSub = (xr ./ a).^2 + (yr ./ b).^2 <= 1;
        return;
    end

    if isprop(roi, 'Position')
        pos = roi.Position;
        if size(pos,2) == 2 && size(pos,1) >= 3
            x1 = max(1, floor(min(pos(:,1))));
            x2 = min(nCols, ceil(max(pos(:,1))));
            y1 = max(1, floor(min(pos(:,2))));
            y2 = min(nRows, ceil(max(pos(:,2))));
            [yg, xg] = ndgrid(y1:y2, x1:x2);
            maskSub = inpolygon(xg, yg, pos(:,1), pos(:,2));
            return;
        end
        if numel(pos) >= 4
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = min(nCols, ceil(pos(1) + pos(3)));
            y2 = min(nRows, ceil(pos(2) + pos(4)));
            maskSub = true(y2-y1+1, x2-x1+1);
            return;
        end
    end
end
