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
    app.distFluofit = [];
    app.flimBayes = [];
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
    app.activeFlimMode = '';
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
    app.defaultDataPath = 'D:\Luminosa\Data\260323\t1_20260323-104842\test\';
    app.smallTitleFont = 8;
    
    % Series (video) related variables
    app.seriesData = [];        % Array of ptuOut structures for each frame
    app.seriesFolderPath = '';  % Path to series folder
    app.seriesFiles = [];       % List of series files
    app.currentFrame = 1;       % Current frame index
    app.frameSlider = [];       % Frame navigation slider
    app.seriesIntensityCache = [];  % Pre-computed intensity images for fast scrolling
    app.seriesTcspcPixCache = [];   % Pre-computed tcspc_pix for fast ROI access
    app.seriesFlimCache = [];       % Pre-computed tau-mean images for fast scrolling
    app.seriesDistFitCache = [];    % Cached DistFluofit-extension results per frame
    app.seriesBayesCache = [];      % Cached Bayesian FLIM results per frame
    app.seriesFrameFileMap = [];    % Map each frame to its source PTU file
    app.ptuOutOriginal = [];        % Original full-file data for global IRF calculation
    app.mietSourceEntries = struct([]);
    app.mietResultEntries = struct([]);
    app.mietCalibration = struct('file', '', 'label', '', 'heightNm', [], 'lifetimeNs', []);
    app.mietDisplayMode = 'height';
    app.axMIET = [];
    app.cbMIET = [];
    app.dropMietDataset = [];
    app.dropMietColormap = [];
    app.txtMietInfo = [];
    app.lblMietStatus = [];
    app.editMietCalibPath = [];
    app.editMietCalibMin = [];
    app.editMietCalibMax = [];
    app.chkMietAutoRange = [];
    app.editMietMin = [];
    app.editMietMax = [];
    app.lastMietExportFile = '';
    app.mietVendorApp = [];

    app.fig = uifigure('Name','Luminosa FLIM / ISM-FLIM / FCS', 'Position',[80 60 1080 680]);
    app.grid = uigridlayout(app.fig, [1 1]);
    app.grid.Padding = [4 4 4 4];
    app.tg = uitabgroup(app.grid);

    tabFLIM = uitab(app.tg, 'Title', 'FLIM / ISM-FLIM');
    tabFCS  = uitab(app.tg, 'Title', 'FCS');
    tabMIET = uitab(app.tg, 'Title', 'MIET');

    buildFLIMTab(tabFLIM);
    buildFCSTab(tabFCS);
    buildMIETTab(tabMIET);

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
        rowHeights = [42, 28, 34, 36, 36, 36, 30, 30, 30, 30, 30, 30, 30, 36, 36, 72];
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
        rowRunPrimary = uigridlayout(ctlGrid, [1 4]);
        rowRunPrimary.ColumnWidth = {'1x','1x','1x','1x'};
        rowRunPrimary.Padding = [0 0 0 0];
        app.btnQuickFLIM = uibutton(rowRunPrimary, 'Text', 'Tau mean', 'ButtonPushedFcn', @onQuickFLIM, 'FontSize', 11);
        app.btnQuickFLIMStd = uibutton(rowRunPrimary, 'Text', 'FLIM std', 'ButtonPushedFcn', @onQuickFLIMStd, 'FontSize', 11);
        app.btnDistFluofit = uibutton(rowRunPrimary, 'Text', 'DistFluofit', 'ButtonPushedFcn', @onDistFluofitExtension, 'FontSize', 10);
        app.btnFlimBayes = uibutton(rowRunPrimary, 'Text', 'FLIM_bayes', 'ButtonPushedFcn', @onFlimBayes, 'FontSize', 10);
        app.flimButtonInactiveColor = app.btnQuickFLIM.BackgroundColor;
        app.flimButtonActiveColor = [0.78 0.78 0.78];
        app.flimButtonInactiveFontColor = app.btnQuickFLIM.FontColor;

        % Row 5
        rowRunSecondary = uigridlayout(ctlGrid, [1 3]);
        rowRunSecondary.ColumnWidth = {'1x','1x','1x'};
        rowRunSecondary.Padding = [0 0 0 0];
        app.btnISMFLIM = uibutton(rowRunSecondary, 'Text', 'ISM', 'ButtonPushedFcn', @onISMFLIM, 'FontSize', 11);
        app.btnShowIntensity = uibutton(rowRunSecondary, 'Text', 'Intensity', 'ButtonPushedFcn', @onShowIntensity, 'FontSize', 11);
        app.btnShowTau = uibutton(rowRunSecondary, 'Text', 'Sum FLIM', 'ButtonPushedFcn', @onShowTauMean, 'FontSize', 11);

        % Row 6
        rowROI = uigridlayout(ctlGrid, [1 3]);
        rowROI.ColumnWidth = {70, 90, '1x'};
        rowROI.Padding = [0 0 0 0];
        app.btnSelectROI = uibutton(rowROI, 'Text', 'ROI', 'ButtonPushedFcn', @onSelectROI, 'FontSize', 11);
        app.dropROIShape = uidropdown(rowROI, 'Items', {'Rectangle','Ellipse','Lasso'}, 'Value', 'Rectangle');
        app.dropROIShape.Tooltip = 'ROI shape';
        app.btnShowTCSPC = uibutton(rowROI, 'Text', 'Show Full TCSPC', 'ButtonPushedFcn', @onShowFullTCSPC, 'FontSize', 11);

        % Row 7
        rowNsub = uigridlayout(ctlGrid, [1 2]);
        rowNsub.ColumnWidth = {120,'1x'};
        rowNsub.Padding = [0 0 0 0];
        uilabel(rowNsub, 'Text', 'Log bins / octave');
        app.txtBins = uieditfield(rowNsub, 'numeric', 'Value', 8, 'Limits', [2 32], ...
            'LowerLimitInclusive', 'on', 'UpperLimitInclusive', 'on', ...
            'ValueChangedFcn', @onTcspcBinningChanged);

        % Row 8
        rowTau0 = uigridlayout(ctlGrid, [1 2]);
        rowTau0.ColumnWidth = {85,'1x'};
        rowTau0.Padding = [0 0 0 0];
        uilabel(rowTau0, 'Text', 'Tau0 (ns)');
        app.editTau0 = uieditfield(rowTau0, 'text', 'Value', '0.35 1.5 5');

        % Row 9
        rowFit = uigridlayout(ctlGrid, [1 2]);
        rowFit.ColumnWidth = {'1x','1x'};
        rowFit.Padding = [0 0 0 0];
        app.chkOptimizeTau = uicheckbox(rowFit, 'Text', 'Optimize tau', 'Value', true);
        app.chkIncludeBG = uicheckbox(rowFit, 'Text', 'Include BG', 'Value', true);

        % Row 10
        rowRange = uigridlayout(ctlGrid, [1 5]);
        rowRange.ColumnWidth = {135,55,'1x','1x',60};
        rowRange.Padding = [0 0 0 0];
        uilabel(rowRange, 'Text', 'Tau range to display');
        app.chkAutoTauRange = uicheckbox(rowRange, 'Text', 'Auto', 'Value', true);
        app.editTauMin = uieditfield(rowRange, 'numeric', 'Value', 1.0, 'Limits', [0 1e3]);
        app.editTauMax = uieditfield(rowRange, 'numeric', 'Value', 7.0, 'Limits', [0 1e3]);
        app.editTauMin.Tooltip = 'Min tau (ns)';
        app.editTauMax.Tooltip = 'Max tau (ns)';
        app.btnApplyTauRange = uibutton(rowRange, 'Text', 'Apply', 'ButtonPushedFcn', @onApplyTauRange);

        % Row 11
        rowPattern = uigridlayout(ctlGrid, [1 2]);
        rowPattern.ColumnWidth = {110,'1x'};
        rowPattern.Padding = [0 0 0 0];
        uilabel(rowPattern, 'Text', 'Pattern source');
        app.dropPatternSource = uidropdown(rowPattern, 'Items', {'Confocal (ptuOut)','ISM reassigned (flim)'}, ...
            'Value', 'Confocal (ptuOut)');

        % Row 12
        rowIRF = uigridlayout(ctlGrid, [1 2]);
        rowIRF.ColumnWidth = {85,'1x'};
        rowIRF.Padding = [0 0 0 0];
        uilabel(rowIRF, 'Text', 'IRF model');
        app.dropIRF = uidropdown(rowIRF, 'Items', {'GammaShifted','ExGauss','Simple'}, ...
            'Value', 'GammaShifted', 'ValueChangedFcn', @onIRFModelChanged);

        % Row 13
        rowChunk = uigridlayout(ctlGrid, [1 2]);
        rowChunk.ColumnWidth = {120,'1x'};
        rowChunk.Padding = [0 0 0 0];
        uilabel(rowChunk, 'Text', 'Read chunk (M)');
        app.editChunkM = uieditfield(rowChunk, 'numeric', 'Value', 1, 'Limits', [0.1 50]);

        % Row 14
        rowAnalysis = uigridlayout(ctlGrid, [1 2]);
        rowAnalysis.ColumnWidth = {120,'1x'};
        rowAnalysis.Padding = [0 0 0 0];
        uibutton(rowAnalysis, 'Text', 'Fit TCSPC', ...
            'ButtonPushedFcn', @onFitDisplayedTCSPC, 'FontSize', 11, ...
            'Tooltip', 'Fit the TCSPC currently shown in the bottom display.');
        app.btnPatternMatch = uibutton(rowAnalysis, 'Text', 'Global Fit + Pattern Match', ...
            'ButtonPushedFcn', @onPatternMatch, 'FontSize', 11);

        % Row 15
        rowSave = uigridlayout(ctlGrid, [1 2]);
        rowSave.ColumnWidth = {'1x','1x'};
        rowSave.Padding = [0 0 0 0];
        app.btnSaveMAT = uibutton(rowSave, 'Text', 'Save MAT', 'ButtonPushedFcn', @onSaveMAT, 'FontSize', 11);
        app.btnSavePNG = uibutton(rowSave, 'Text', 'Save PNG', 'ButtonPushedFcn', @onSavePNG, 'FontSize', 11);

        % Row 16
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
        setActiveFlimMode('');
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

    function buildMIETTab(parent)
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

        rowHeights = [42, 42, 34, 30, 30, 30, 34, 36, 36, 36, 90];
        ctlGrid = uigridlayout(ctlInner, [numel(rowHeights) 1]);
        ctlGrid.RowHeight = num2cell(rowHeights);
        ctlGrid.RowSpacing = 6;
        ctlGrid.Padding = [8 8 8 8];

        innerH = sum(rowHeights) + ctlGrid.RowSpacing*(numel(rowHeights)-1) + ctlGrid.Padding(2) + ctlGrid.Padding(4) + 8;
        ctlInner.Position = [0 0 415 innerH];

        rowSource = uigridlayout(ctlGrid, [1 2]);
        rowSource.ColumnWidth = {'1x','1x'};
        rowSource.Padding = [0 0 0 0];
        uibutton(rowSource, 'Text', 'Use Current Exports', 'ButtonPushedFcn', @onMietUseCurrentExports, 'FontSize', 11);
        uibutton(rowSource, 'Text', 'Load Export MAT', 'ButtonPushedFcn', @onMietLoadExportMAT, 'FontSize', 11);

        rowVendor = uigridlayout(ctlGrid, [1 2]);
        rowVendor.ColumnWidth = {'1x','1x'};
        rowVendor.Padding = [0 0 0 0];
        uibutton(rowVendor, 'Text', 'Load Calibration', 'ButtonPushedFcn', @onMietLoadCalibration, 'FontSize', 11);
        uibutton(rowVendor, 'Text', 'Launch MIET-GUI', 'ButtonPushedFcn', @onMietLaunchVendor, 'FontSize', 11);

        rowDataset = uigridlayout(ctlGrid, [1 2]);
        rowDataset.ColumnWidth = {95, '1x'};
        rowDataset.Padding = [0 0 0 0];
        uilabel(rowDataset, 'Text', 'Dataset');
        app.dropMietDataset = uidropdown(rowDataset, ...
            'Items', {'No MIET source loaded'}, 'Value', 'No MIET source loaded', ...
            'ValueChangedFcn', @onMietDatasetChanged);

        rowCalibPath = uigridlayout(ctlGrid, [1 2]);
        rowCalibPath.ColumnWidth = {95,'1x'};
        rowCalibPath.Padding = [0 0 0 0];
        uilabel(rowCalibPath, 'Text', 'Calibration');
        app.editMietCalibPath = uieditfield(rowCalibPath, 'text', 'Value', '', 'Editable', 'off');

        rowCalibRange = uigridlayout(ctlGrid, [1 3]);
        rowCalibRange.ColumnWidth = {95,'1x','1x'};
        rowCalibRange.Padding = [0 0 0 0];
        uilabel(rowCalibRange, 'Text', 'Calib z (nm)');
        app.editMietCalibMin = uieditfield(rowCalibRange, 'numeric', 'Value', 0, 'Limits', [-1e6 1e6]);
        app.editMietCalibMax = uieditfield(rowCalibRange, 'numeric', 'Value', 250, 'Limits', [-1e6 1e6]);

        rowDisplayRange = uigridlayout(ctlGrid, [1 4]);
        rowDisplayRange.ColumnWidth = {80,55,'1x','1x'};
        rowDisplayRange.Padding = [0 0 0 0];
        uilabel(rowDisplayRange, 'Text', 'Height range');
        app.chkMietAutoRange = uicheckbox(rowDisplayRange, 'Text', 'Auto', 'Value', true, 'ValueChangedFcn', @onMietApplyRange);
        app.editMietMin = uieditfield(rowDisplayRange, 'numeric', 'Value', 0, 'Limits', [-1e6 1e6], 'ValueChangedFcn', @onMietApplyRange);
        app.editMietMax = uieditfield(rowDisplayRange, 'numeric', 'Value', 200, 'Limits', [-1e6 1e6], 'ValueChangedFcn', @onMietApplyRange);

        rowColormap = uigridlayout(ctlGrid, [1 3]);
        rowColormap.ColumnWidth = {80,'1x',110};
        rowColormap.Padding = [0 0 0 0];
        uilabel(rowColormap, 'Text', 'Colormap');
        app.dropMietColormap = uidropdown(rowColormap, ...
            'Items', {'MIET Bronze', 'Parula', 'Hot', 'Copper', 'Gray', 'Bone', 'Jet'}, ...
            'Value', 'MIET Bronze', ...
            'ValueChangedFcn', @onMietColormapChanged);
        uibutton(rowColormap, 'Text', 'Nice Image', 'ButtonPushedFcn', @onMietNiceImage, 'FontSize', 11);

        rowCompute = uigridlayout(ctlGrid, [1 2]);
        rowCompute.ColumnWidth = {'1x','1x'};
        rowCompute.Padding = [0 0 0 0];
        uibutton(rowCompute, 'Text', 'Compute Height', 'ButtonPushedFcn', @onMietCompute, 'FontSize', 11);
        uibutton(rowCompute, 'Text', 'Save MIET MAT', 'ButtonPushedFcn', @onMietSaveMAT, 'FontSize', 11);

        rowDisplayMode = uigridlayout(ctlGrid, [1 2]);
        rowDisplayMode.ColumnWidth = {'1x','1x'};
        rowDisplayMode.Padding = [0 0 0 0];
        uibutton(rowDisplayMode, 'Text', 'Show Lifetime', 'ButtonPushedFcn', @onMietShowLifetime, 'FontSize', 11);
        uibutton(rowDisplayMode, 'Text', 'Show Height', 'ButtonPushedFcn', @onMietShowHeight, 'FontSize', 11);

        app.lblMietStatus = uilabel(ctlGrid, 'Text', 'Ready.');
        app.lblMietStatus.WordWrap = 'on';

        app.txtMietInfo = uitextarea(ctlGrid, 'Editable', 'off', ...
            'Value', {'No MIET source loaded.'});

        axGrid = uigridlayout(axp, [1 1]);
        axGrid.Padding = [6 6 6 6];
        app.axMIET = uiaxes(axGrid);
        title(app.axMIET, 'MIET');
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
            app.distFluofit = [];
            app.flimBayes = [];
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
            app.seriesFolderPath = '';
            app.seriesFiles = [];
            app.seriesIntensityCache = {};
            app.seriesTcspcPixCache = {};
            app.seriesFlimCache = {};
            app.seriesDistFitCache = {};
            app.seriesBayesCache = {};
            app.seriesFrameFileMap = [];
            app.currentFrame = 1;
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
                            app.seriesDistFitCache = cell(numFrames, 1);
                            app.seriesBayesCache = cell(numFrames, 1);
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
                                app.seriesFlimCache{frameIdx} = [];
                                app.seriesDistFitCache{frameIdx} = [];
                                app.seriesBayesCache{frameIdx} = [];
                                
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
        app.distFluofit = [];
        app.flimBayes = [];
        app.seriesDistFitCache = {};
        app.seriesBayesCache = {};
        
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
        allDistFitCache = {};
        allBayesCache = {};
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
                        allFlimCache{totalFrameCount} = [];
                        allDistFitCache{totalFrameCount} = [];
                        allBayesCache{totalFrameCount} = [];
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
                    allDistFitCache{totalFrameCount} = [];
                    allBayesCache{totalFrameCount} = [];
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
        app.seriesDistFitCache = allDistFitCache;
        app.seriesBayesCache = allBayesCache;
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
            app.distFluofit = [];
            app.flimBayes = [];
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
            cachedTauStd = getTauStdMapFromPTUData(app.ptuOut);
            if ~isempty(app.seriesFlimCache) && frameNum <= numel(app.seriesFlimCache)
                cachedTau = app.seriesFlimCache{frameNum};
            end
            if ~isempty(cachedTau)
                app.flim = flimStructFromTauMap(cachedTau, getIntensityMapFromPTUData(app.ptuOut), cachedTauStd);
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

            if ~isempty(app.seriesDistFitCache) && frameNum <= numel(app.seriesDistFitCache)
                app.distFluofit = app.seriesDistFitCache{frameNum};
            end
            if ~isempty(app.seriesBayesCache) && frameNum <= numel(app.seriesBayesCache)
                app.flimBayes = app.seriesBayesCache{frameNum};
            end

            if strcmp(requestedDisplayMode, 'file_tau')
                showFileSummaryOverlay();
            elseif strcmp(requestedDisplayMode, 'tau') && ~isempty(app.flim)
                showTauMean();
            elseif strcmp(requestedDisplayMode, 'tau_std')
                showTauStd();
            elseif strcmp(requestedDisplayMode, 'distfluofit')
                if ~isempty(app.distFluofit)
                    showDistFluofitOverlay();
                else
                    onDistFluofitExtension();
                end
            elseif strcmp(requestedDisplayMode, 'bayes')
                if ~isempty(app.flimBayes)
                    showFlimBayesOverlay();
                else
                    onFlimBayes();
                end
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

        if ~ensureFlimMetric('tauMean')
            return;
        end

        showTauMean();
    end

    function onQuickFLIMStd(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        if ~ensureFlimMetric('tauRMS')
            return;
        end

        showTauStd();
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

    function onDistFluofitExtension(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        cacheIdx = [];
        if ~isempty(app.seriesData) && ~isempty(app.currentFrame)
            cacheIdx = app.currentFrame;
            if cacheIdx <= numel(app.seriesDistFitCache) && ~isempty(app.seriesDistFitCache{cacheIdx})
                app.distFluofit = app.seriesDistFitCache{cacheIdx};
                addStatus(sprintf('Using cached DistFluofit result for frame %d.', cacheIdx));
                showDistFluofitOverlay();
                return;
            end
        end

        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        [tcspc_pix, dtNs, tcspcSrc, cubeErr] = resolvePatternMatchTcspcCube();
        if isempty(tcspc_pix)
            addStatus(cubeErr);
            return;
        end

        pulsePeriodNs = app.ptuOut.head.MeasDesc_GlobalResolution * 1e9;
        [irf, irfMeta, irfNote] = getWholeFileIRFForTargetGrid(dtNs, size(tcspc_pix, 3));
        if isempty(irf)
            addStatus('Whole-file IRF unavailable for DistFluofit extension.');
            return;
        end
        if ~isempty(irfNote) && ~strcmp(irfNote, 'native')
            addStatus(sprintf('DistFluofit using %s whole-file IRF.', irfNote));
        end

        addStatus(sprintf('Running DistFluofit extension with whole-file IRF using %s...', tcspcSrc));
        drawnow;

        opts = struct();
        opts.useGPU = useGPU;
        opts.mode = 'PIRLS';
        opts.nTau = 100;
        opts.shiftBounds = [-10 10];
        opts.batchSize = 1024;

        try
            outDist = distfluofit_extension(tcspc_pix, irf, pulsePeriodNs, dtNs, opts);
            outDist.irf = irf(:);
            outDist.irfMeta = irfMeta;
            app.distFluofit = outDist;
            if ~isempty(cacheIdx)
                app.seriesDistFitCache{cacheIdx} = outDist;
            end
        catch ME
            addStatus(['DistFluofit extension failed: ' ME.message]);
            return;
        end

        showDistFluofitOverlay();
        addStatus(sprintf('DistFluofit extension complete (%s, %s).', tcspcSrc, outDist.patternMatchInfo.method));
    end

    function onFlimBayes(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5" for FLIM_bayes.');
            return;
        end

        tau0Requested = sort(double(tau0(:)).', 'ascend');
        if numel(tau0Requested) > 2
            addStatus(sprintf('FLIM_bayes uses a two-state low-photon model; using shortest/longest Tau0 values [%.4g %.4g] ns.', ...
                tau0Requested(1), tau0Requested(end)));
            tau0Use = [tau0Requested(1), tau0Requested(end)];
        else
            tau0Use = tau0Requested;
        end

        cacheIdx = [];
        if ~isempty(app.seriesData) && ~isempty(app.currentFrame)
            cacheIdx = app.currentFrame;
            if cacheIdx <= numel(app.seriesBayesCache) && ~isempty(app.seriesBayesCache{cacheIdx})
                app.flimBayes = app.seriesBayesCache{cacheIdx};
                addStatus(sprintf('Using cached FLIM_bayes result for frame %d.', cacheIdx));
                showFlimBayesOverlay();
                return;
            end
        end

        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        useGPU = app.chkUseGPU.Value && gpuIsAvailable();
        [tcspc_pix, dtNs, tcspcSrc, cubeErr] = resolvePatternMatchTcspcCube();
        if isempty(tcspc_pix)
            addStatus(cubeErr);
            return;
        end

        pulsePeriodNs = app.ptuOut.head.MeasDesc_GlobalResolution * 1e9;
        [irf, irfMeta, irfNote] = getWholeFileIRFForTargetGrid(dtNs, size(tcspc_pix, 3));
        if isempty(irf)
            addStatus('Whole-file IRF unavailable for FLIM_bayes.');
            return;
        end
        if ~isempty(irfNote) && ~strcmp(irfNote, 'native')
            addStatus(sprintf('FLIM_bayes using %s whole-file IRF.', irfNote));
        end

        addStatus(sprintf('Running FLIM_bayes with whole-file IRF using %s...', tcspcSrc));
        drawnow;

        opts = struct();
        opts.useGPU = useGPU;
        opts.batchSize = 2048;
        opts.includeBackground = app.chkIncludeBG.Value;
        opts.optimizeTau = app.chkOptimizeTau.Value;
        opts.signalGrid = linspace(0.0, 1.0, 26);
        opts.fractionGrid = linspace(0.0, 1.0, 41);
        opts.shiftBounds = [-5 5];
        opts.singleExpTauGrid = [];

        try
            outBayes = flim_bayes_lowphoton(tcspc_pix, irf, pulsePeriodNs, dtNs, tau0Use, opts);
            outBayes.irf = irf(:);
            outBayes.irfMeta = irfMeta;
            outBayes.tau0Requested = tau0Requested(:).';
            app.flimBayes = outBayes;
            if ~isempty(cacheIdx)
                app.seriesBayesCache{cacheIdx} = outBayes;
            end
        catch ME
            addStatus(['FLIM_bayes failed: ' ME.message]);
            return;
        end

        showFlimBayesOverlay();
        addStatus(sprintf('FLIM_bayes complete (%s, %s).', tcspcSrc, outBayes.posteriorInfo.method));
    end

    function onShowIntensity(~, ~)
        showIntensityFromPTU();
    end

    function onShowTauMean(~, ~)
        showFileSummaryOverlay();
    end

    function onGammaChanged(~, ~)
        app.gamma = app.sliderGamma.Value;
        app.lblGamma.Text = sprintf('%.2f', app.gamma);
        switch app.displayMode
            case 'tau'
                showTauMean();
            case 'tau_std'
                showTauStd();
            case 'file_tau'
                showFileSummaryOverlay();
            case 'pattern'
                showPatternOverlay();
            case 'distfluofit'
                showDistFluofitOverlay();
            case 'bayes'
                showFlimBayesOverlay();
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
        elseif strcmp(app.displayMode, 'tau_std')
            showTauStd();
        elseif strcmp(app.displayMode, 'file_tau')
            showFileSummaryOverlay();
        elseif strcmp(app.displayMode, 'pattern')
            showPatternOverlay();
        elseif strcmp(app.displayMode, 'distfluofit')
            showDistFluofitOverlay();
        elseif strcmp(app.displayMode, 'bayes')
            showFlimBayesOverlay();
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
        tcspcDensityDisp = clipIrfDisplay(cDensity, 1e-1);
        maskTcspc = isfinite(tBinNs(:)) & isfinite(tcspcDensityDisp(:));
        plot(app.axTCSPC, tBinNs(maskTcspc)*1e-9, tcspcDensityDisp(maskTcspc), 'k.', 'MarkerSize', 10);
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
                if plotClippedIrfDensity(app.axTCSPC, tIRFNs, irfDensity, 1e-1, 'b--', 1.1)
                    legendLines{end+1} = 'Whole-file IRF (scaled)';
                end
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
                       sprintf('Display shift = %.3f ns (full measured curve shown)', shiftNs), ...
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
        tcspcDensityDisp = clipIrfDisplay(app.tcspc.cDensity, 1e-1);
        maskTcspc = isfinite(app.tcspc.tBinNs(:)) & isfinite(tcspcDensityDisp(:));
        plot(app.axTCSPC, app.tcspc.tBinNs(maskTcspc)*1e-9, tcspcDensityDisp(maskTcspc), 'k.', 'MarkerSize', 10);
        hold(app.axTCSPC, 'on');
        fitDensityDisp = clipIrfDisplay(fitDensity, 1e-1);
        maskFit = isfinite(tFitNs(:)) & isfinite(fitDensityDisp(:));
        plot(app.axTCSPC, tFitNs(maskFit)*1e-9, fitDensityDisp(maskFit), 'r-', 'LineWidth', 1.3);

        legendLines = {'ROI TCSPC', 'Reconvolution fit'};
        irfDisp = applyDisplayShiftToDecay(irfROI(:), app.tcspc.shiftBins);
        irfDisp = max(irfDisp, 0);
        if ~isempty(irfDisp) && max(irfDisp) > 0
            irfScale = max(app.tcspc.rawCountsShift) / max(irfDisp);
            irfScaled = irfDisp * irfScale;
            [tIRFNs, irfDensity] = logBinTCSPC(irfScaled, app.tcspc.dtNs, app.txtBins.Value);
            if plotClippedIrfDensity(app.axTCSPC, tIRFNs, irfDensity, 1e-1, 'b--', 1.1)
                legendLines{end+1} = 'Whole-file IRF (scaled)';
            end
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
        summaryLines = {sprintf('ROI fit | model: %d-exp', numel(tauFit)), ...
                        sprintf('tau (ns): [%s]', strtrim(num2str(tauFit(:).', '%.4g '))), ...
                        sprintf('amp frac: [%s]', strtrim(num2str(ampFrac(:).', '%.4g ')))};
        if includeBG
            bgOffset = fittedBackgroundOffsetPerBin(coeff, numel(app.tcspc.rawCountsNative), includeBG);
            summaryLines{end+1} = sprintf('bg offset ~= %.4g counts/bin', bgOffset);
        end
        summaryLines = [summaryLines, {sprintf('chi2/ndf = %.4g', chi2red), ...
                        sprintf('IRF model: %s', app.dropIRF.Value), ...
                        sprintf('IRF source: whole file | dt = %.3f ps', 1e3*app.tcspc.dtNs)}];
        setFitSummary(summaryLines);

        app.tcspcDisplayMode = 'roi_fit';
        addStatus(sprintf('ROI fit tau (ns): %s | whole-file IRF model: %s', num2str(tauFit(:).', '%.3g '), app.dropIRF.Value));
    end

    function onFitDisplayedTCSPC(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end

        switch app.tcspcDisplayMode
            case {'roi', 'roi_fit'}
                onFitTCSPC();
            case 'global'
                onFitGlobalTCSPC();
            otherwise
                if ~isempty(app.roi) && isvalid(app.roi)
                    onShowTCSPC();
                    onFitTCSPC();
                else
                    showGlobalTCSPC();
                    onFitGlobalTCSPC();
                end
            end
    end

    function onFitGlobalTCSPC(~, ~)
        if isempty(app.ptuOut)
            addStatus('Load a PTU first.');
            return;
        end
        setBusy(true);
        cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>

        if ~ensureGlobalIRF(false) || isempty(app.irfGlobal)
            addStatus('Whole-file IRF unavailable for fitting the displayed TCSPC.');
            return;
        end

        tau0 = parseTau0(app.editTau0.Value);
        if isempty(tau0)
            addStatus('Enter tau0 (ns), e.g. "0.35 1.5 5".');
            return;
        end

        includeBG = app.chkIncludeBG.Value;
        optimizeTau = app.chkOptimizeTau.Value;

        [countsNative, dtNs, srcInfo] = tcspcFromWholeFile_native(app.ptuOut);
        if isempty(countsNative) || isempty(dtNs)
            addStatus('Whole-file TCSPC unavailable for fitting.');
            return;
        end

        addStatus(sprintf('Fitting displayed whole-file TCSPC as %d-exp using %s IRF...', numel(tau0), app.dropIRF.Value));
        try
            [fitOut, summaryLines] = fitWholeDecayModelWithTau0(countsNative, dtNs, app.irfGlobal, tau0, includeBG, optimizeTau);
        catch ME
            addStatus(['Whole-file TCSPC fit failed: ' ME.message]);
            return;
        end

        if isempty(fitOut)
            addStatus('Whole-file TCSPC fit failed.');
            return;
        end

        fitOut.summaryLines = summaryLines;
        fitOut.source = srcInfo;
        [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(countsNative, dtNs); %#ok<ASGLU>
        tcspcGlobalTmp = struct('dtNs', dtNs, 'rawCountsNative', countsNative(:), 'rawCountsShift', countsShift(:), ...
            'shiftBins', shiftBins, 'shiftNs', shiftNs, 'peakIdx', peakIdx, 'riseIdx', riseIdx, 'cBin', [], 'cDensity', []);
        fitOut = refreshGlobalFitDisplayCache(fitOut, tcspcGlobalTmp, app.txtBins.Value);
        app.globalFit = fitOut;

        showGlobalTCSPC();
        addStatus(sprintf('Whole-file TCSPC fit updated using the %d-component model from Tau0.', app.globalFit.nExp));
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
        distFluofit = app.distFluofit; %#ok<NASGU>
        flimBayes = app.flimBayes; %#ok<NASGU>
        tcspc = app.tcspc; %#ok<NASGU>
        tcspcFit = app.tcspcFit; %#ok<NASGU>
        frameMapExports = buildFrameMapExports(); %#ok<NASGU>
        fileSummaryMapExports = buildFileSummaryMapExports(); %#ok<NASGU>
        exportLabels = struct( ... %#ok<NASGU>
            'frameMapsVariable', 'frameMapExports', ...
            'fileSummaryVariable', 'fileSummaryMapExports', ...
            'intensityField', 'intensity', ...
            'lifetimeField', 'lifetimeNs', ...
            'lifetimeUnits', 'ns');
        save(outFile, 'ptuOut', 'flim', 'ismRes', 'pattern', 'distFluofit', 'flimBayes', 'tcspc', 'tcspcFit', ...
            'frameMapExports', 'fileSummaryMapExports', 'exportLabels', '-v7.3');
        app.lastMietExportFile = outFile;
        addStatus(sprintf('Saved: %s | %d frame map(s), %d file summary map(s).', ...
            outFile, numel(frameMapExports), numel(fileSummaryMapExports)));
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

    function onMietUseCurrentExports(~, ~)
        if isempty(app.ptuOut)
            setMietStatus('Load a PTU or series first, or load a saved export MAT.');
            return;
        end

        entries = combineMietSourceEntries(buildFrameMapExports(), buildFileSummaryMapExports());
        if isempty(entries)
            setMietStatus('No lifetime/intensity exports are available for MIET.');
            return;
        end

        app.mietSourceEntries = entries;
        app.mietResultEntries = struct([]);
        refreshMietDatasetItems();
        app.mietDisplayMode = 'lifetime';
        showCurrentMietMap();
        setMietStatus(sprintf('Loaded %d MIET source map(s) from the current Luminosa session.', numel(entries)));
    end

    function onMietLoadExportMAT(~, ~)
        [f, p] = uigetfile('*.mat', 'Load Luminosa export MAT for MIET');
        if isequal(f, 0)
            return;
        end

        fileName = fullfile(p, f);
        try
            s = load(fileName, 'frameMapExports', 'fileSummaryMapExports');
        catch ME
            setMietStatus(['Failed to load MIET export MAT: ' ME.message]);
            return;
        end

        frameMapExports = [];
        fileSummaryMapExports = [];
        if isfield(s, 'frameMapExports')
            frameMapExports = s.frameMapExports;
        end
        if isfield(s, 'fileSummaryMapExports')
            fileSummaryMapExports = s.fileSummaryMapExports;
        end

        entries = combineMietSourceEntries(frameMapExports, fileSummaryMapExports);
        if isempty(entries)
            setMietStatus('Selected MAT file does not contain frameMapExports or fileSummaryMapExports.');
            return;
        end

        app.mietSourceEntries = entries;
        app.mietResultEntries = struct([]);
        app.lastMietExportFile = fileName;
        refreshMietDatasetItems();
        app.mietDisplayMode = 'lifetime';
        showCurrentMietMap();
        setMietStatus(sprintf('Loaded %d MIET source map(s) from %s.', numel(entries), f));
    end

    function onMietLoadCalibration(~, ~)
        [f, p] = uigetfile({'*.mat;*.csv;*.txt;*.dat', 'Calibration files (*.mat,*.csv,*.txt,*.dat)'}, ...
            'Load MIET calibration');
        if isequal(f, 0)
            return;
        end

        fileName = fullfile(p, f);
        [calib, errMsg] = loadMietCalibrationFile(fileName);
        if ~isempty(errMsg)
            setMietStatus(errMsg);
            return;
        end

        app.mietCalibration = calib;
        app.editMietCalibPath.Value = fileName;
        if ~isempty(calib.heightNm)
            app.editMietCalibMin.Value = min(calib.heightNm);
            app.editMietCalibMax.Value = max(calib.heightNm);
        end
        setMietStatus(sprintf('Loaded MIET calibration: %s (%d samples).', calib.label, numel(calib.heightNm)));
        setMietInfoLines(buildMietInfoLines());
    end

    function onMietLaunchVendor(~, ~)
        vendorRoot = luminosa_miet_vendor_root();
        if isempty(vendorRoot)
            setMietStatus('MIET-GUI was not found on the MATLAB path. Put it in external/miet-gui or add it to the MATLAB path first.');
            return;
        end

        candidates = {'Luminosa_MIET_gui', 'MIET_gui', 'MIETGUI', 'miet_gui', 'MIET_GUI'};
        for ii = 1:numel(candidates)
            name = candidates{ii};
            if exist(name, 'class') == 8 || exist(name, 'file') == 2
                try
                    switch name
                        case 'Luminosa_MIET_gui'
                            app.mietVendorApp = Luminosa_MIET_gui();
                        case 'MIET_gui'
                            app.mietVendorApp = MIET_gui();
                        otherwise
                            app.mietVendorApp = [];
                            feval(name);
                    end
                    setMietStatus(sprintf('Launched external MIET-GUI via %s.', name));
                    return;
                catch ME
                    setMietStatus(sprintf('Found %s, but launch failed: %s', name, ME.message));
                    return;
                end
            end
        end

        setMietStatus(sprintf('MIET-GUI launch entry point was not found after adding %s to the MATLAB path.', vendorRoot));
    end

    function onMietCompute(~, ~)
        if isempty(app.mietSourceEntries)
            setMietStatus('Load current exports or a saved export MAT before computing MIET.');
            return;
        end

        [calibHeightNm, calibLifetimeNs, calibMeta, errMsg] = getActiveMietCalibration();
        if ~isempty(errMsg)
            setMietStatus(errMsg);
            return;
        end

        results = app.mietSourceEntries;
        for idx = 1:numel(results)
            lifetimeMap = [];
            if isfield(results(idx), 'lifetimeNs')
                lifetimeMap = results(idx).lifetimeNs;
            end
            heightMap = lifetimeMapToHeightMap(lifetimeMap, calibLifetimeNs, calibHeightNm);
            results(idx).heightNm = heightMap;
            results(idx).heightLabel = 'MIET height (nm)';
            results(idx).calibrationLabel = calibMeta.label;
            results(idx).calibrationFile = calibMeta.file;
            results(idx).calibrationHeightRangeNm = calibMeta.heightRangeNm;
            results(idx).calibrationLifetimeRangeNs = calibMeta.lifetimeRangeNs;
        end

        app.mietResultEntries = results;
        app.mietDisplayMode = 'height';
        refreshMietDatasetItems();
        showCurrentMietMap();
        setMietStatus(sprintf('Computed MIET height maps for %d dataset(s) using %s.', numel(results), calibMeta.label));
    end

    function onMietSaveMAT(~, ~)
        if isempty(app.mietSourceEntries) && isempty(app.mietResultEntries)
            setMietStatus('No MIET data is available to save.');
            return;
        end

        [f, p] = uiputfile('*.mat', 'Save MIET MAT');
        if isequal(f, 0)
            return;
        end

        outFile = fullfile(p, f);
        mietSourceEntries = app.mietSourceEntries; %#ok<NASGU>
        mietResultEntries = app.mietResultEntries; %#ok<NASGU>
        mietCalibration = app.mietCalibration; %#ok<NASGU>
        save(outFile, 'mietSourceEntries', 'mietResultEntries', 'mietCalibration', '-v7.3');
        setMietStatus(['Saved MIET results: ' outFile]);
    end

    function onMietShowLifetime(~, ~)
        app.mietDisplayMode = 'lifetime';
        showCurrentMietMap();
    end

    function onMietShowHeight(~, ~)
        app.mietDisplayMode = 'height';
        showCurrentMietMap();
    end

    function onMietDatasetChanged(~, ~)
        showCurrentMietMap();
    end

    function onMietColormapChanged(~, ~)
        showCurrentMietMap();
    end

    function onMietApplyRange(~, ~)
        if strcmp(app.mietDisplayMode, 'height')
            showCurrentMietMap();
        end
    end

    function onMietNiceImage(~, ~)
        [mapData, intensityMap, titleStr, modeName, valueRange, errMsg] = getCurrentMietViewData();
        if ~isempty(errMsg)
            setMietStatus(errMsg);
            return;
        end

        cmapName = getSelectedMietColormapName();
        try
            luminosa_miet_nice_image(intensityMap, mapData, titleStr, ...
                'Mode', modeName, ...
                'ColormapName', cmapName, ...
                'ValueRange', valueRange);
            setMietStatus(sprintf('Opened %s nice image using %s.', modeName, cmapName));
        catch ME
            setMietStatus(sprintf('Failed to create nice image: %s', ME.message));
        end
    end

    function setMietStatus(msg)
        if isempty(msg)
            msg = 'Ready.';
        end
        if ~isempty(app.lblMietStatus) && isvalid(app.lblMietStatus)
            app.lblMietStatus.Text = msg;
        end
        addStatus(['MIET: ' msg]);
    end

    function setMietInfoLines(lines)
        if nargin < 1 || isempty(lines)
            lines = {'No MIET source loaded.'};
        end
        if ischar(lines)
            lines = cellstr(lines);
        end
        if ~isempty(app.txtMietInfo) && isvalid(app.txtMietInfo)
            app.txtMietInfo.Value = lines(:);
        end
    end

    function entries = combineMietSourceEntries(frameMapExports, fileSummaryMapExports)
        entries = struct([]);
        if nargin < 1 || isempty(frameMapExports)
            frameMapExports = struct([]);
        end
        if nargin < 2 || isempty(fileSummaryMapExports)
            fileSummaryMapExports = struct([]);
        end

        frameMapExports = reshape(frameMapExports, [], 1);
        fileSummaryMapExports = reshape(fileSummaryMapExports, [], 1);
        if isempty(frameMapExports)
            entries = fileSummaryMapExports;
        elseif isempty(fileSummaryMapExports)
            entries = frameMapExports;
        else
            entries = [frameMapExports; fileSummaryMapExports];
        end
    end

    function entries = getActiveMietEntries()
        if ~isempty(app.mietResultEntries)
            entries = app.mietResultEntries;
        else
            entries = app.mietSourceEntries;
        end
    end

    function refreshMietDatasetItems()
        entries = getActiveMietEntries();
        if isempty(entries)
            app.dropMietDataset.Items = {'No MIET source loaded'};
            app.dropMietDataset.Value = 'No MIET source loaded';
            setMietInfoLines({'No MIET source loaded.'});
            cla(app.axMIET);
            clearMietColorbar();
            return;
        end

        items = {entries.label};
        prevValue = app.dropMietDataset.Value;
        app.dropMietDataset.Items = items;
        if any(strcmp(items, prevValue))
            app.dropMietDataset.Value = prevValue;
        else
            app.dropMietDataset.Value = items{1};
        end
        setMietInfoLines(buildMietInfoLines());
    end

    function idx = currentMietEntryIndex()
        entries = getActiveMietEntries();
        if isempty(entries)
            idx = [];
            return;
        end
        value = app.dropMietDataset.Value;
        items = {entries.label};
        idx = find(strcmp(items, value), 1, 'first');
        if isempty(idx)
            idx = 1;
        end
    end

    function lines = buildMietInfoLines()
        lines = {'No MIET source loaded.'};
        entries = getActiveMietEntries();
        idx = currentMietEntryIndex();
        if isempty(entries) || isempty(idx)
            return;
        end

        entry = entries(idx);
        lines = { ...
            sprintf('Label: %s', getStructText(entry, 'label', '')), ...
            sprintf('Type: %s', getStructText(entry, 'mapType', '')), ...
            sprintf('Source file: %s', getStructText(entry, 'sourceFile', '')), ...
            sprintf('Frames represented: %s', scalarFieldText(entry, 'nFrames')), ...
            sprintf('Frame index in file: %s', scalarFieldText(entry, 'fileFrameIndex'))};

        if isfield(entry, 'calibrationLabel') && ~isempty(entry.calibrationLabel)
            lines{end+1} = sprintf('Calibration: %s', entry.calibrationLabel);
        elseif ~isempty(app.mietCalibration.label)
            lines{end+1} = sprintf('Calibration: %s', app.mietCalibration.label);
        else
            lines{end+1} = 'Calibration: not loaded';
        end

        if isfield(entry, 'lifetimeNs') && ~isempty(entry.lifetimeNs)
            tauVals = double(entry.lifetimeNs(isfinite(entry.lifetimeNs)));
            if ~isempty(tauVals)
                lines{end+1} = sprintf('Lifetime range (ns): %.3f .. %.3f', min(tauVals), max(tauVals));
            end
        end

        if isfield(entry, 'heightNm') && ~isempty(entry.heightNm)
            hVals = double(entry.heightNm(isfinite(entry.heightNm)));
            if ~isempty(hVals)
                lines{end+1} = sprintf('Height range (nm): %.3f .. %.3f', min(hVals), max(hVals));
            end
        end
    end

    function txt = getStructText(s, fieldName, defaultValue)
        if nargin < 3
            defaultValue = '';
        end
        if isfield(s, fieldName) && ~isempty(s.(fieldName))
            txt = s.(fieldName);
        else
            txt = defaultValue;
        end
    end

    function txt = scalarFieldText(s, fieldName)
        txt = 'n/a';
        if isfield(s, fieldName) && ~isempty(s.(fieldName))
            val = s.(fieldName);
            if isscalar(val) && isfinite(double(val))
                txt = sprintf('%.0f', double(val));
            end
        end
    end

    function showCurrentMietMap()
        [mapData, intensityMap, titleStr, modeName, ~, errMsg] = getCurrentMietViewData();
        if ~isempty(errMsg)
            return;
        end

        showMietScalarMap(mapData, intensityMap, titleStr, modeName);
        setMietInfoLines(buildMietInfoLines());
    end

    function showMietScalarMap(mapData, intensityMap, titleStr, modeName)
        if isempty(app.axMIET) || ~isvalid(app.axMIET) || isempty(mapData)
            return;
        end

        mapData = double(mapData);
        clearMietColorbar();
        cla(app.axMIET);

        valueRange = getMietDisplayRange(mapData, modeName);
        cmap = luminosa_miet_colormap(getSelectedMietColormapName(), 256, 'display');
        if strcmp(modeName, 'height')
            cbarLabel = 'Height (nm)';
        else
            cbarLabel = 'Lifetime (ns)';
        end

        if nargin >= 2 && ~isempty(intensityMap) && isequal(size(intensityMap), size(mapData))
            bright = applyGamma(double(intensityMap), app.gamma);
            rgb = flim_rgb(mapData, bright, valueRange, cmap);
            image(app.axMIET, rgb);
            axis(app.axMIET, 'image');
            axis(app.axMIET, 'off');
        else
            imagesc(app.axMIET, mapData);
            axis(app.axMIET, 'image');
            axis(app.axMIET, 'on');
        end
        colormap(app.axMIET, cmap);
        caxis(app.axMIET, valueRange);
        app.cbMIET = colorbar(app.axMIET);
        ylabel(app.cbMIET, cbarLabel);
        title(app.axMIET, titleStr, 'Interpreter', 'none');
    end

    function [mapData, intensityMap, titleStr, modeName, valueRange, errMsg] = getCurrentMietViewData()
        mapData = [];
        intensityMap = [];
        titleStr = 'MIET';
        modeName = 'height';
        valueRange = [0 1];
        errMsg = 'No MIET source loaded.';

        entries = getActiveMietEntries();
        idx = currentMietEntryIndex();
        if isempty(entries) || isempty(idx)
            return;
        end

        entry = entries(idx);
        if isfield(entry, 'intensity')
            intensityMap = entry.intensity;
        end

        if strcmp(app.mietDisplayMode, 'height') && isfield(entry, 'heightNm') && ~isempty(entry.heightNm)
            mapData = entry.heightNm;
            titleStr = sprintf('MIET height: %s', entry.label);
            modeName = 'height';
        elseif isfield(entry, 'lifetimeNs') && ~isempty(entry.lifetimeNs)
            mapData = entry.lifetimeNs;
            titleStr = sprintf('Lifetime source: %s', entry.label);
            modeName = 'lifetime';
        else
            errMsg = 'Selected dataset does not contain a displayable MIET map.';
            return;
        end

        valueRange = getMietDisplayRange(mapData, modeName);
        errMsg = '';
    end

    function cmapName = getSelectedMietColormapName()
        cmapName = 'MIET Bronze';
        if ~isempty(app.dropMietColormap) && isvalid(app.dropMietColormap) && ~isempty(app.dropMietColormap.Value)
            cmapName = app.dropMietColormap.Value;
        end
    end

    function valueRange = getMietDisplayRange(mapData, modeName)
        if strcmp(modeName, 'height')
            valueRange = getMietHeightRange(mapData);
        else
            valueRange = getTauRange(mapData);
        end
    end

    function clearMietColorbar()
        if ~isempty(app.cbMIET) && isvalid(app.cbMIET)
            delete(app.cbMIET);
        end
        app.cbMIET = [];
    end

    function trange = getMietHeightRange(heightMap)
        vals = double(heightMap(isfinite(heightMap)));
        if isempty(vals)
            trange = [0 1];
            return;
        end

        if app.chkMietAutoRange.Value
            trange = prctile(vals, [5 95]);
        else
            trange = [app.editMietMin.Value, app.editMietMax.Value];
        end

        if trange(1) == trange(2)
            trange = [min(vals) max(vals)];
        end
        if trange(2) <= trange(1)
            trange(2) = trange(1) + eps;
        end
    end

    function [calib, errMsg] = loadMietCalibrationFile(fileName)
        calib = struct('file', fileName, 'label', '', 'heightNm', [], 'lifetimeNs', []);
        errMsg = '';

        [~, baseName, ext] = fileparts(fileName);
        ext = lower(ext);

        try
            switch ext
                case '.mat'
                    s = load(fileName);
                    [heightNm, lifetimeNs, label] = parseMietCalibrationFromLoadedData(s, baseName);
                case {'.csv', '.txt', '.dat'}
                    [heightNm, lifetimeNs, label] = parseMietCalibrationFromTextFile(fileName, baseName);
                otherwise
                    errMsg = sprintf('Unsupported calibration format: %s', ext);
                    return;
            end
        catch ME
            errMsg = sprintf('Failed to load calibration: %s', ME.message);
            return;
        end

        if isempty(heightNm) || isempty(lifetimeNs)
            errMsg = 'Could not find height/lifetime calibration vectors in the selected file.';
            return;
        end

        [heightNm, order] = sort(double(heightNm(:)), 'ascend');
        lifetimeNs = double(lifetimeNs(:));
        lifetimeNs = lifetimeNs(order);
        valid = isfinite(heightNm) & isfinite(lifetimeNs);
        heightNm = heightNm(valid);
        lifetimeNs = lifetimeNs(valid);
        if numel(heightNm) < 2
            errMsg = 'Calibration file contains fewer than two valid samples.';
            return;
        end

        calib.label = label;
        calib.heightNm = heightNm;
        calib.lifetimeNs = lifetimeNs;
    end

    function [heightNm, lifetimeNs, label] = parseMietCalibrationFromLoadedData(s, defaultLabel)
        heightNm = [];
        lifetimeNs = [];
        label = defaultLabel;

        [heightNm, lifetimeNs] = parseMietCalibrationCandidate(s);
        if ~isempty(heightNm)
            return;
        end

        fields = fieldnames(s);
        for ii = 1:numel(fields)
            candidate = s.(fields{ii});
            [heightNm, lifetimeNs] = parseMietCalibrationCandidate(candidate);
            if ~isempty(heightNm)
                label = sprintf('%s:%s', defaultLabel, fields{ii});
                return;
            end
        end
    end

    function [heightNm, lifetimeNs] = parseMietCalibrationCandidate(candidate)
        heightNm = [];
        lifetimeNs = [];

        if istable(candidate)
            names = lower(candidate.Properties.VariableNames);
            heightIdx = find(contains(names, 'height') | contains(names, 'distance') | contains(names, 'znm') | strcmp(names, 'z'), 1, 'first');
            lifetimeIdx = find(contains(names, 'lifetime') | contains(names, 'tau'), 1, 'first');
            if ~isempty(heightIdx) && ~isempty(lifetimeIdx)
                heightNm = candidate{:, heightIdx};
                lifetimeNs = candidate{:, lifetimeIdx};
            end
            return;
        end

        if isstruct(candidate)
            heightFields = {'heightNm', 'height', 'distanceNm', 'distance', 'zNm', 'z'};
            lifetimeFields = {'lifetimeNs', 'lifetime', 'tauNs', 'tau'};
            h = [];
            t = [];
            for ii = 1:numel(heightFields)
                if isfield(candidate, heightFields{ii})
                    h = candidate.(heightFields{ii});
                    break;
                end
            end
            for ii = 1:numel(lifetimeFields)
                if isfield(candidate, lifetimeFields{ii})
                    t = candidate.(lifetimeFields{ii});
                    break;
                end
            end
            if ~isempty(h) && ~isempty(t)
                heightNm = h;
                lifetimeNs = t;
            end
            return;
        end

        if isnumeric(candidate) && ismatrix(candidate) && size(candidate, 2) >= 2
            heightNm = candidate(:, 1);
            lifetimeNs = candidate(:, 2);
        end
    end

    function [heightNm, lifetimeNs, label] = parseMietCalibrationFromTextFile(fileName, defaultLabel)
        heightNm = [];
        lifetimeNs = [];
        label = defaultLabel;

        try
            tbl = readtable(fileName);
            [heightNm, lifetimeNs] = parseMietCalibrationCandidate(tbl);
        catch
            heightNm = [];
            lifetimeNs = [];
        end

        if isempty(heightNm) || isempty(lifetimeNs)
            mat = readmatrix(fileName);
            if size(mat, 2) < 2
                error('Calibration text file must have at least two columns.');
            end
            heightNm = mat(:, 1);
            lifetimeNs = mat(:, 2);
            label = sprintf('%s (col1=height nm, col2=lifetime ns)', defaultLabel);
        end
    end

    function [heightNm, lifetimeNs, meta, errMsg] = getActiveMietCalibration()
        heightNm = [];
        lifetimeNs = [];
        meta = struct('label', '', 'file', '', 'heightRangeNm', [], 'lifetimeRangeNs', []);
        errMsg = '';

        if isempty(app.mietCalibration.heightNm) || isempty(app.mietCalibration.lifetimeNs)
            errMsg = 'Load a MIET calibration first.';
            return;
        end

        zMin = min(app.editMietCalibMin.Value, app.editMietCalibMax.Value);
        zMax = max(app.editMietCalibMin.Value, app.editMietCalibMax.Value);
        mask = app.mietCalibration.heightNm >= zMin & app.mietCalibration.heightNm <= zMax;
        if nnz(mask) < 2
            errMsg = 'Calibration height window contains fewer than two samples.';
            return;
        end

        heightNm = double(app.mietCalibration.heightNm(mask));
        lifetimeNs = double(app.mietCalibration.lifetimeNs(mask));

        if all(diff(lifetimeNs) >= 0)
            % already monotonic
        elseif all(diff(lifetimeNs) <= 0)
            lifetimeNs = flipud(lifetimeNs(:));
            heightNm = flipud(heightNm(:));
        else
            errMsg = 'Selected calibration height window is not monotonic in lifetime. Restrict the z-range to a monotonic branch before computing MIET.';
            return;
        end

        [lifetimeNs, ia] = unique(lifetimeNs(:), 'stable');
        heightNm = heightNm(ia);
        if numel(lifetimeNs) < 2
            errMsg = 'Calibration needs at least two unique lifetime samples after cropping.';
            return;
        end

        meta.label = app.mietCalibration.label;
        meta.file = app.mietCalibration.file;
        meta.heightRangeNm = [min(heightNm) max(heightNm)];
        meta.lifetimeRangeNs = [min(lifetimeNs) max(lifetimeNs)];
    end

    function heightMap = lifetimeMapToHeightMap(lifetimeMap, calibLifetimeNs, calibHeightNm)
        heightMap = [];
        if isempty(lifetimeMap)
            return;
        end

        lifetimeMap = double(lifetimeMap);
        heightMap = nan(size(lifetimeMap), 'double');
        valid = isfinite(lifetimeMap);
        if ~any(valid(:))
            return;
        end

        if numel(calibLifetimeNs) >= 3
            method = 'pchip';
        else
            method = 'linear';
        end

        heightMap(valid) = interp1(calibLifetimeNs, calibHeightNm, lifetimeMap(valid), method, NaN);
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
        floorVal = 1e-1;
        try
            drawnow limitrate;
        catch
        end
        ax.YScale = 'log';
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        yMax = findMaxPositiveYData(ax, floorVal);
        if ~isfinite(yMax) || yMax <= floorVal
            yLim = ax.YLim;
            if numel(yLim) >= 2 && isfinite(yLim(2)) && yLim(2) > floorVal
                yMax = yLim(2);
            else
                yMax = 10 * floorVal;
            end
        end
        ax.YLim = [floorVal yMax];
        p1 = floor(log10(floorVal));
        p2 = ceil(log10(yMax));
        if p2 < p1
            p2 = p1;
        end
        yTicks = 10.^(p1:p2);
        yTicks = yTicks(yTicks >= floorVal & yTicks <= yMax * (1 + 1e-12));
        if isempty(yTicks) || yTicks(1) > floorVal * (1 + 1e-12)
            yTicks = [floorVal, yTicks];
        end
        ax.YTick = unique(yTicks);
    end

    function yMax = findMaxPositiveYData(ax, floorVal)
        yMax = NaN;
        if nargin < 2 || isempty(floorVal)
            floorVal = 1e-1;
        end
        if isempty(ax) || ~isvalid(ax)
            return;
        end

        kids = ax.Children;
        for k = 1:numel(kids)
            child = kids(k);
            if ~isprop(child, 'YData')
                continue;
            end
            try
                y = double(child.YData(:));
            catch
                continue;
            end
            y = y(isfinite(y) & y >= floorVal);
            if isempty(y)
                continue;
            end
            ymaxChild = max(y);
            if ~isfinite(yMax) || ymaxChild > yMax
                yMax = ymaxChild;
            end
        end
    end

    function y = clipIrfDisplay(y, floorVal)
        if nargin < 2 || isempty(floorVal)
            floorVal = 1e-1;
        end
        y = double(y);
        y(~isfinite(y) | y < floorVal) = NaN;
    end

    function didPlot = plotClippedIrfDensity(ax, tNs, yDensity, floorVal, lineSpec, lineWidth)
        didPlot = false;
        if nargin < 4 || isempty(floorVal)
            floorVal = 1e-1;
        end
        if nargin < 5 || isempty(lineSpec)
            lineSpec = 'b--';
        end
        if nargin < 6 || isempty(lineWidth)
            lineWidth = 1.1;
        end
        if isempty(ax) || ~isvalid(ax)
            return;
        end

        tNs = double(tNs(:));
        yDensity = clipIrfDisplay(yDensity, floorVal);
        mask = isfinite(tNs) & isfinite(yDensity) & (yDensity > 0);
        if ~any(mask)
            return;
        end

        plot(ax, tNs(mask)*1e-9, yDensity(mask), lineSpec, 'LineWidth', lineWidth);
        didPlot = true;
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
        setActiveFlimMode('');
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
        setActiveFlimMode('');
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

        if isfield(ptuData, 'tauMean') && ~isempty(ptuData.tauMean)
            tauCube = double(ptuData.tauMean);
        elseif isfield(ptuData, 'meanArrival') && ~isempty(ptuData.meanArrival)
            tauCube = double(ptuData.meanArrival);
        elseif isfield(ptuData, 'total') && isfield(ptuData.total, 'tauMean') && ~isempty(ptuData.total.tauMean)
            tauCube = double(ptuData.total.tauMean);
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

    function tauStd = getTauStdMapFromPTUData(ptuData)
        tauStd = [];
        tauCube = [];
        tagCube = [];

        if isfield(ptuData, 'tauRMS') && ~isempty(ptuData.tauRMS)
            tauCube = double(ptuData.tauRMS);
        elseif isfield(ptuData, 'tauStd') && ~isempty(ptuData.tauStd)
            tauCube = double(ptuData.tauStd);
        elseif isfield(ptuData, 'total') && isfield(ptuData.total, 'tauRMS') && ~isempty(ptuData.total.tauRMS)
            tauCube = double(ptuData.total.tauRMS);
        elseif isfield(ptuData, 'taus') && ~isempty(ptuData.taus)
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
            tauStd = double(tauCube);
            return;
        end

        if ~isempty(tagCube) && size(tagCube, 1) == size(tauCube, 1) && ...
                size(tagCube, 2) == size(tauCube, 2) && size(tagCube, 3) == size(tauCube, 3)
            denom = sum(tagCube, 3);
            tauStd = sum(tauCube .* tagCube, 3) ./ max(denom, 1);
            zeroMask = denom <= 0;
            if any(zeroMask(:))
                tauFallback = mean(tauCube, 3, 'omitnan');
                tauStd(zeroMask) = tauFallback(zeroMask);
            end
        else
            tauStd = mean(tauCube, 3, 'omitnan');
        end
    end

    function tauMean = computeTauMeanMapOnDemand(ptuData)
        tauMean = getTauMeanMapFromPTUData(ptuData);
        if ~isempty(tauMean) || isempty(ptuData) || ~hasQuickFLIMData(ptuData)
            return;
        end

        try
            useGPU = app.chkUseGPU.Value && gpuIsAvailable();
            flimTmp = quickFLIMFromTCSPCFlexible(ptuData, useGPU);
            if isfield(flimTmp, 'total') && isfield(flimTmp.total, 'tauMean')
                tauMean = flimTmp.total.tauMean;
            end
        catch
            tauMean = [];
        end
    end

    function [tauMap, intensityMap, titleStr, statusMsg] = resolveFileSummaryOverlay()
        tauMap = [];
        intensityMap = [];
        titleStr = 'Summed intensity + FLIM overlay';
        statusMsg = '';

        if ~isempty(app.ptuOutOriginal)
            intensityMap = getIntensityMapFromPTUData(app.ptuOutOriginal);
            tauMap = computeTauMeanMapOnDemand(app.ptuOutOriginal);
            if ~isempty(tauMap) && ~isempty(intensityMap)
                sourceName = currentSourceFileName();
                titleStr = buildFileSummaryTitle(sourceName);
                statusMsg = sprintf('Showing summed intensity + FLIM overlay for all frames in %s.', sourceName);
                return;
            end
        end

        if ~isempty(app.seriesData)
            [tauMap, intensityMap, sourceName, nFramesUsed] = aggregateCurrentSourceFileOverlay();
            if ~isempty(tauMap) && ~isempty(intensityMap)
                titleStr = buildFileSummaryTitle(sourceName);
                statusMsg = sprintf('Showing summed intensity + FLIM overlay for %d frame(s) in %s.', nFramesUsed, sourceName);
                return;
            end
        end

        if ~isempty(app.ptuOut)
            intensityMap = getIntensityMapFromPTUData(app.ptuOut);
            tauMap = computeTauMeanMapOnDemand(app.ptuOut);
            if ~isempty(tauMap) && ~isempty(intensityMap)
                sourceName = currentSourceFileName();
                titleStr = buildFileSummaryTitle(sourceName);
                statusMsg = sprintf('Showing summed intensity + FLIM overlay for %s.', sourceName);
            end
        end
    end

    function [tauMap, intensityMap, sourceName, nFramesUsed] = aggregateCurrentSourceFileOverlay()
        tauMap = [];
        intensityMap = [];
        sourceName = currentSourceFileName();
        nFramesUsed = 0;

        if isempty(app.seriesData)
            return;
        end

        frameIdx = [];
        if ~isempty(app.seriesFrameFileMap) && ~isempty(app.currentFrame) && ...
                app.currentFrame >= 1 && app.currentFrame <= numel(app.seriesFrameFileMap)
            sourceIdx = app.seriesFrameFileMap(app.currentFrame);
            if isfinite(sourceIdx)
                frameIdx = find(app.seriesFrameFileMap == sourceIdx);
                if ~isempty(app.seriesFiles) && sourceIdx >= 1 && sourceIdx <= numel(app.seriesFiles)
                    sourceName = app.seriesFiles(sourceIdx).name;
                end
            end
        end

        if isempty(frameIdx)
            frameIdx = 1:numel(app.seriesData);
        end

        nFramesUsed = numel(frameIdx);
        [tauMap, intensityMap] = aggregateFrameOverlayMaps(frameIdx);
    end

    function [tauMap, intensityMap] = aggregateFrameOverlayMaps(frameIdx)
        tauMap = [];
        intensityMap = [];
        weightedTau = [];
        tauWeight = [];
        tauFallback = [];
        tauCount = [];

        for ii = 1:numel(frameIdx)
            idx = frameIdx(ii);
            if idx < 1 || idx > numel(app.seriesData)
                continue;
            end

            if idx <= numel(app.seriesIntensityCache)
                frameIntensity = app.seriesIntensityCache{idx};
            else
                frameIntensity = [];
            end
            if isempty(frameIntensity) && ~isempty(app.seriesData{idx})
                frameIntensity = getIntensityMapFromPTUData(app.seriesData{idx});
            end
            if isempty(frameIntensity)
                continue;
            end

            frameIntensity = double(frameIntensity);
            if isempty(intensityMap)
                intensityMap = zeros(size(frameIntensity), 'double');
                weightedTau = zeros(size(frameIntensity), 'double');
                tauWeight = zeros(size(frameIntensity), 'double');
                tauFallback = zeros(size(frameIntensity), 'double');
                tauCount = zeros(size(frameIntensity), 'double');
            elseif ~isequal(size(frameIntensity), size(intensityMap))
                continue;
            end

            intensityMap = intensityMap + frameIntensity;

            if idx <= numel(app.seriesFlimCache)
                frameTau = app.seriesFlimCache{idx};
            else
                frameTau = [];
            end
            if isempty(frameTau) && ~isempty(app.seriesData{idx})
                frameTau = computeTauMeanMapOnDemand(app.seriesData{idx});
            end
            if isempty(frameTau) || ~isequal(size(frameTau), size(intensityMap))
                continue;
            end

            frameTau = double(frameTau);
            finiteTau = isfinite(frameTau);
            if any(finiteTau(:))
                tauFallback(finiteTau) = tauFallback(finiteTau) + frameTau(finiteTau);
                tauCount(finiteTau) = tauCount(finiteTau) + 1;
            end

            validWeighted = finiteTau & isfinite(frameIntensity) & frameIntensity > 0;
            if any(validWeighted(:))
                weightedTau(validWeighted) = weightedTau(validWeighted) + frameTau(validWeighted) .* frameIntensity(validWeighted);
                tauWeight(validWeighted) = tauWeight(validWeighted) + frameIntensity(validWeighted);
            end
        end

        if isempty(intensityMap)
            return;
        end

        tauMap = nan(size(intensityMap), 'double');
        weightedMask = tauWeight > 0;
        tauMap(weightedMask) = weightedTau(weightedMask) ./ tauWeight(weightedMask);

        fallbackMask = ~weightedMask & tauCount > 0;
        if any(fallbackMask(:))
            tauMap(fallbackMask) = tauFallback(fallbackMask) ./ tauCount(fallbackMask);
        end
    end

    function sourceName = currentSourceFileName()
        sourceName = '';
        if ~isempty(app.seriesData) && ~isempty(app.seriesFrameFileMap) && ~isempty(app.currentFrame) && ...
                app.currentFrame >= 1 && app.currentFrame <= numel(app.seriesFrameFileMap)
            sourceIdx = app.seriesFrameFileMap(app.currentFrame);
            if isfinite(sourceIdx) && ~isempty(app.seriesFiles) && sourceIdx >= 1 && sourceIdx <= numel(app.seriesFiles)
                sourceName = app.seriesFiles(sourceIdx).name;
            end
        end

        if isempty(sourceName) && ~isempty(app.lastFile)
            [~, name, ext] = fileparts(app.lastFile);
            sourceName = [name ext];
        end

        if isempty(sourceName)
            sourceName = 'the active file';
        end
    end

    function titleStr = buildFileSummaryTitle(sourceName)
        if nargin < 1 || isempty(sourceName)
            titleStr = 'Summed intensity + FLIM overlay';
        else
            titleStr = sprintf('Summed intensity + FLIM overlay (%s)', sourceName);
        end
    end

    function frameMapExports = buildFrameMapExports()
        frameMapExports = emptyMapExport();

        if ~isempty(app.seriesData)
            outIdx = 0;
            for idx = 1:numel(app.seriesData)
                frameData = app.seriesData{idx};
                if isempty(frameData)
                    continue;
                end

                intensityMap = [];
                if idx <= numel(app.seriesIntensityCache)
                    intensityMap = app.seriesIntensityCache{idx};
                end
                if isempty(intensityMap)
                    intensityMap = getIntensityMapFromPTUData(frameData);
                end

                lifetimeMap = [];
                if idx <= numel(app.seriesFlimCache)
                    lifetimeMap = app.seriesFlimCache{idx};
                end
                if isempty(lifetimeMap)
                    lifetimeMap = computeTauMeanMapOnDemand(frameData);
                end

                if isempty(intensityMap) && isempty(lifetimeMap)
                    continue;
                end

                [sourceName, sourceIdx, fileFrameIndex] = frameSourceInfo(idx);
                outIdx = outIdx + 1;
                frameMapExports(outIdx) = makeMapExportEntry( ...
                    sprintf('%s | frame %d', sourceName, fileFrameIndex), ...
                    'per_frame', sourceName, sourceIdx, idx, fileFrameIndex, 1, intensityMap, lifetimeMap);
            end
            return;
        end

        intensityMap = getIntensityMapFromPTUData(app.ptuOut);
        lifetimeMap = computeTauMeanMapOnDemand(app.ptuOut);
        if isempty(intensityMap) && isempty(lifetimeMap)
            return;
        end

        sourceName = currentSourceFileName();
        frameMapExports(1) = makeMapExportEntry( ...
            sprintf('%s | frame 1', sourceName), ...
            'per_frame', sourceName, NaN, 1, 1, 1, intensityMap, lifetimeMap);
    end

    function fileSummaryMapExports = buildFileSummaryMapExports()
        fileSummaryMapExports = emptyMapExport();

        if ~isempty(app.ptuOutOriginal)
            intensityMap = getIntensityMapFromPTUData(app.ptuOutOriginal);
            lifetimeMap = computeTauMeanMapOnDemand(app.ptuOutOriginal);
            if isempty(intensityMap) && isempty(lifetimeMap)
                return;
            end

            sourceName = currentSourceFileName();
            nFrames = max(1, numel(app.seriesData));
            fileSummaryMapExports(1) = makeMapExportEntry( ...
                sprintf('%s | total summed', sourceName), ...
                'summed_file', sourceName, 1, [], [], nFrames, intensityMap, lifetimeMap);
            return;
        end

        if ~isempty(app.seriesData)
            if ~isempty(app.seriesFrameFileMap)
                sourceIdxList = unique(app.seriesFrameFileMap(:).', 'stable');
            else
                sourceIdxList = NaN;
            end

            outIdx = 0;
            for ii = 1:numel(sourceIdxList)
                sourceIdx = sourceIdxList(ii);
                if isfinite(sourceIdx)
                    frameIdx = find(app.seriesFrameFileMap == sourceIdx);
                else
                    frameIdx = 1:numel(app.seriesData);
                end

                [lifetimeMap, intensityMap] = aggregateFrameOverlayMaps(frameIdx);
                if isempty(intensityMap) && isempty(lifetimeMap)
                    continue;
                end

                sourceName = sourceFileNameFromIndex(sourceIdx);
                outIdx = outIdx + 1;
                fileSummaryMapExports(outIdx) = makeMapExportEntry( ...
                    sprintf('%s | total summed', sourceName), ...
                    'summed_file', sourceName, sourceIdx, [], [], numel(frameIdx), intensityMap, lifetimeMap);
            end
            return;
        end

        intensityMap = getIntensityMapFromPTUData(app.ptuOut);
        lifetimeMap = computeTauMeanMapOnDemand(app.ptuOut);
        if isempty(intensityMap) && isempty(lifetimeMap)
            return;
        end

        sourceName = currentSourceFileName();
        fileSummaryMapExports(1) = makeMapExportEntry( ...
            sprintf('%s | total summed', sourceName), ...
            'summed_file', sourceName, NaN, [], [], 1, intensityMap, lifetimeMap);
    end

    function [sourceName, sourceIdx, fileFrameIndex] = frameSourceInfo(globalFrameIndex)
        sourceIdx = NaN;
        fileFrameIndex = globalFrameIndex;

        if ~isempty(app.seriesFrameFileMap) && globalFrameIndex >= 1 && globalFrameIndex <= numel(app.seriesFrameFileMap)
            sourceIdx = double(app.seriesFrameFileMap(globalFrameIndex));
            if isfinite(sourceIdx)
                sourceFrames = find(app.seriesFrameFileMap == sourceIdx);
                localIdx = find(sourceFrames == globalFrameIndex, 1, 'first');
                if ~isempty(localIdx)
                    fileFrameIndex = localIdx;
                end
            end
        end

        sourceName = sourceFileNameFromIndex(sourceIdx);
    end

    function sourceName = sourceFileNameFromIndex(sourceIdx)
        sourceName = '';

        if nargin >= 1 && isfinite(sourceIdx) && ~isempty(app.seriesFiles) && ...
                sourceIdx >= 1 && sourceIdx <= numel(app.seriesFiles)
            sourceName = app.seriesFiles(sourceIdx).name;
        end

        if isempty(sourceName) && ~isempty(app.lastFile)
            [~, name, ext] = fileparts(app.lastFile);
            sourceName = [name ext];
        end

        if isempty(sourceName)
            sourceName = 'the active file';
        end
    end

    function entry = makeMapExportEntry(label, mapType, sourceName, sourceIdx, globalFrameIndex, fileFrameIndex, nFrames, intensityMap, lifetimeMap)
        if strcmp(mapType, 'summed_file')
            intensityLabel = 'summed intensity';
            lifetimeLabel = 'summed lifetime (ns)';
        else
            intensityLabel = 'frame intensity';
            lifetimeLabel = 'frame lifetime (ns)';
        end

        entry = struct( ...
            'label', label, ...
            'mapType', mapType, ...
            'sourceFile', sourceName, ...
            'sourceFileIndex', sourceIdx, ...
            'globalFrameIndex', globalFrameIndex, ...
            'fileFrameIndex', fileFrameIndex, ...
            'nFrames', nFrames, ...
            'intensityLabel', intensityLabel, ...
            'lifetimeLabel', lifetimeLabel, ...
            'intensity', doubleOrEmpty(intensityMap), ...
            'lifetimeNs', doubleOrEmpty(lifetimeMap));
    end

    function exports = emptyMapExport()
        exports = struct( ...
            'label', {}, ...
            'mapType', {}, ...
            'sourceFile', {}, ...
            'sourceFileIndex', {}, ...
            'globalFrameIndex', {}, ...
            'fileFrameIndex', {}, ...
            'nFrames', {}, ...
            'intensityLabel', {}, ...
            'lifetimeLabel', {}, ...
            'intensity', {}, ...
            'lifetimeNs', {});
    end

    function arr = doubleOrEmpty(arrIn)
        if isempty(arrIn)
            arr = [];
        else
            arr = double(arrIn);
        end
    end

    function flim = flimStructFromTauMap(tauMean, intensityImg, tauStd)
        flim = [];
        if nargin < 3
            tauStd = [];
        end
        if isempty(tauMean) && isempty(tauStd)
            return;
        end

        if isempty(tauMean)
            tauMean = nan(size(tauStd), 'double');
        end

        flim = struct();
        flim.total = struct();
        flim.total.tauMean = double(tauMean);
        flim.total.meanArrival = double(tauMean);
        flim.total.tauRMS = double(tauStd);
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
        setActiveFlimMode('tau');
    end

    function showTauStd()
        img = [];
        if ~isempty(app.flim) && isfield(app.flim, 'total') && isfield(app.flim.total, 'tauRMS') && ...
                ~isempty(app.flim.total.tauRMS)
            img = app.flim.total.tauRMS;
        elseif ~isempty(app.flim) && isfield(app.flim, 'reassigned') && isfield(app.flim.reassigned, 'total') && ...
                isfield(app.flim.reassigned.total, 'tauRMS') && ~isempty(app.flim.reassigned.total.tauRMS)
            img = app.flim.reassigned.total.tauRMS;
        elseif ~isempty(app.ptuOut)
            img = getTauStdMapFromPTUData(app.ptuOut);
        end

        if isempty(img)
            addStatus('No FLIM std result available.');
            return;
        end

        intensity = getIntensityMap();
        if isempty(intensity)
            intensity = getIntensityMapFromPTUData(app.ptuOut);
        end

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
            ylabel(app.cbTau, 'Tau std (ns)');
            title(app.axImage, 'FLIM std (ns)');
            setSmallTitles();
        else
            showTauOverlay(img, intensity, 'FLIM std (ns)', 'Tau std (ns)');
        end
        app.displayMode = 'tau_std';
        setActiveFlimMode('tau_std');
    end

    function showFileSummaryOverlay()
        [tauMap, intensityMap, titleStr, statusMsg] = resolveFileSummaryOverlay();
        if isempty(tauMap) || isempty(intensityMap)
            addStatus('No summed file FLIM overlay is available.');
            return;
        end
        showTauOverlay(tauMap, intensityMap, titleStr);
        app.displayMode = 'file_tau';
        setActiveFlimMode('');
        if ~isempty(statusMsg)
            addStatus(statusMsg);
        end
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
        setActiveFlimMode('');
    end

    function showDistFluofitOverlay()
        if isempty(app.distFluofit) || ~isfield(app.distFluofit, 'tauMeanArithmetic') || isempty(app.distFluofit.tauMeanArithmetic)
            addStatus('No DistFluofit result available.');
            return;
        end
        tau = app.distFluofit.tauMeanArithmetic;
        intensity = app.distFluofit.intensity;
        showTauOverlay(tau, intensity, 'DistFluofit tau mean (ns)');
        app.displayMode = 'distfluofit';
        setActiveFlimMode('distfluofit');
    end

    function showFlimBayesOverlay()
        if isempty(app.flimBayes) || ~isfield(app.flimBayes, 'tauMeanArithmetic') || isempty(app.flimBayes.tauMeanArithmetic)
            addStatus('No FLIM_bayes result available.');
            return;
        end
        tau = app.flimBayes.tauMeanArithmetic;
        intensity = app.flimBayes.intensity;
        showTauOverlay(tau, intensity, 'FLIM_bayes tau mean (ns)');
        app.displayMode = 'bayes';
        setActiveFlimMode('bayes');
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
        tcspcDensityDisp = clipIrfDisplay(app.tcspcGlobal.cDensity, 1e-1);
        maskTcspc = isfinite(app.tcspcGlobal.tBinNs(:)) & isfinite(tcspcDensityDisp(:));
        plot(app.axTCSPC, app.tcspcGlobal.tBinNs(maskTcspc)*1e-9, tcspcDensityDisp(maskTcspc), 'k.', 'MarkerSize', 10);
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
                if plotClippedIrfDensity(app.axTCSPC, tIRFNs, irfDensity, 1e-1, 'b--', 1.1)
                    legendLines{end+1} = 'Whole-file IRF (scaled)';
                end
            end
        end

        if ~isempty(fitDensityPlot)
            fitDensityDisp = clipIrfDisplay(fitDensityPlot, 1e-1);
            maskFit = isfinite(app.globalFit.tBinNs(:)) & isfinite(fitDensityDisp(:));
            plot(app.axTCSPC, app.globalFit.tBinNs(maskFit)*1e-9, fitDensityDisp(maskFit), 'r-', 'LineWidth', 1.3);
            fitLegend = sprintf('Fit (%d-exp)', app.globalFit.nExp);
            if isfield(app.globalFit, 'selectionMode')
                selectionMode = lower(char(string(app.globalFit.selectionMode)));
                switch selectionMode
                    case 'bic'
                        fitLegend = sprintf('Best fit (%d-exp)', app.globalFit.nExp);
                    case 'tau0'
                        fitLegend = sprintf('Tau0-selected fit (%d-exp)', app.globalFit.nExp);
                end
            end
            legendLines{end+1} = fitLegend;
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
        bestFit.selectionMode = 'bic';
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

    function showTauOverlay(tau, intensity, titleStr, colorbarLabel)
        cmap = jet(256);
        cmap = cmap(30:end-30,:);

        if nargin < 4 || isempty(colorbarLabel)
            colorbarLabel = 'Tau (ns)';
        end

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
        ylabel(app.cbTau, colorbarLabel);
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

    function setActiveFlimMode(modeName)
        if nargin < 1 || isempty(modeName)
            modeName = '';
        end
        app.activeFlimMode = char(modeName);

        btnList = {
            app.btnQuickFLIM, 'tau';
            app.btnQuickFLIMStd, 'tau_std';
            app.btnDistFluofit, 'distfluofit';
            app.btnFlimBayes, 'bayes'};

        for ii = 1:size(btnList, 1)
            btn = btnList{ii, 1};
            btnMode = btnList{ii, 2};
            if isempty(btn) || ~isvalid(btn)
                continue;
            end

            isActive = strcmp(app.activeFlimMode, btnMode);
            if isActive
                btn.BackgroundColor = app.flimButtonActiveColor;
                btn.FontWeight = 'bold';
                btn.FontColor = [0.15 0.15 0.15];
            else
                btn.BackgroundColor = app.flimButtonInactiveColor;
                btn.FontWeight = 'normal';
                btn.FontColor = app.flimButtonInactiveFontColor;
            end
        end
    end

    function ok = hasFlimMetric(flim, metricName)
        ok = false;
        if isempty(flim) || ~isstruct(flim)
            return;
        end

        if isfield(flim, 'total') && isfield(flim.total, metricName) && ~isempty(flim.total.(metricName))
            ok = any(isfinite(double(flim.total.(metricName)(:))));
            return;
        end

        if isfield(flim, 'reassigned') && isfield(flim.reassigned, 'total') && ...
                isfield(flim.reassigned.total, metricName) && ~isempty(flim.reassigned.total.(metricName))
            ok = any(isfinite(double(flim.reassigned.total.(metricName)(:))));
        end
    end

    function ok = ensureFlimMetric(metricName)
        ok = false;
        if isempty(metricName)
            return;
        end

        switch metricName
            case 'tauMean'
                statusLabel = 'tau mean';
            case 'tauRMS'
                statusLabel = 'FLIM std';
            otherwise
                statusLabel = metricName;
        end

        if hasFlimMetric(app.flim, metricName)
            ok = true;
            return;
        end

        if strcmp(metricName, 'tauRMS')
            tauStd = getTauStdMapFromPTUData(app.ptuOut);
            if ~isempty(tauStd)
                tauMean = [];
                if hasFlimMetric(app.flim, 'tauMean') && isfield(app.flim, 'total')
                    tauMean = app.flim.total.tauMean;
                end
                app.flim = flimStructFromTauMap(tauMean, getIntensityMapFromPTUData(app.ptuOut), tauStd);
                ok = true;
                return;
            end
        end

        if hasQuickFLIMData(app.ptuOut)
            setBusy(true);
            cleanupBusy = onCleanup(@() setBusy(false)); %#ok<NASGU>
            addStatus(sprintf('Computing %s...', statusLabel));
            drawnow;

            useGPU = app.chkUseGPU.Value && gpuIsAvailable();
            try
                app.flim = quickFLIMFromTCSPCFlexible(app.ptuOut, useGPU);
            catch ME
                addStatus(sprintf('%s computation failed: %s', statusLabel, ME.message));
                app.flim = [];
                return;
            end

            ok = hasFlimMetric(app.flim, metricName);
            if ok
                addStatus(sprintf('%s done.', statusLabel));
                return;
            end
        end

        if strcmp(metricName, 'tauMean')
            tauMean = getTauMeanMapFromPTUData(app.ptuOut);
            if ~isempty(tauMean)
                app.flim = flimStructFromTauMap(tauMean, getIntensityMapFromPTUData(app.ptuOut), getTauStdMapFromPTUData(app.ptuOut));
                ok = true;
                addStatus('Using cached tau-mean map.');
                return;
            end
            addStatus('No frame-specific tau-mean data available.');
        else
            addStatus('No FLIM std data available.');
        end
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

function [fitOut, summaryLines] = fitWholeDecayModelWithTau0(countsRawFull, dtNs, irf, tau0, includeBG, optimizeTau)
    tau0 = double(tau0(:).');
    if isempty(tau0)
        fitOut = [];
        summaryLines = {};
        return;
    end

    [tauFit, coeff, fitCountsRawFull, stats] = fitDecayIRFModelMultiStart(countsRawFull, dtNs, tau0, irf, includeBG, optimizeTau);
    fitOut = struct();
    fitOut.nExp = numel(tau0);
    fitOut.tau0 = tau0;
    fitOut.tauFit = tauFit(:).';
    fitOut.coeff = coeff(:).';
    fitOut.fitCountsRawFull = fitCountsRawFull(:);
    fitOut.stats = stats;
    fitOut.selectionMode = 'tau0';

    amps = fitOut.coeff(1 + double(fitOut.stats.includeBG):end);
    ampFrac = amps ./ max(sum(amps), eps);
    summaryLines = { ...
        sprintf('Selected model: %d-exp (from Tau0)', fitOut.nExp), ...
        sprintf('tau0 (ns): [%s]', strtrim(num2str(fitOut.tau0, '%.4g '))), ...
        sprintf('tau fit (ns): [%s]', strtrim(num2str(fitOut.tauFit, '%.4g '))), ...
        sprintf('amp frac: [%s]', strtrim(num2str(ampFrac(:).', '%.4g '))), ...
        sprintf('chi2/ndf = %.4g', fitOut.stats.chi2red)};
    if fitOut.stats.includeBG
        summaryLines{end+1} = sprintf('bg offset ~= %.4g counts/bin', ...
            fittedBackgroundOffsetPerBin(fitOut.coeff, numel(fitOut.fitCountsRawFull), fitOut.stats.includeBG));
        summaryLines{end+1} = 'Background term included';
    else
        summaryLines{end+1} = 'Background term excluded';
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
    if bestFit.stats.includeBG
        summaryLines = [summaryLines(1:3), ...
            {sprintf('bg offset ~= %.4g counts/bin', ...
            fittedBackgroundOffsetPerBin(bestFit.coeff, numel(bestFit.fitCountsRawFull), bestFit.stats.includeBG))}, ...
            summaryLines(4:end)];
    end
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
        tShiftSq3 = gpuArray(reshape(tShift.^2, 1, 1, []));
        if ndims(cubeG) == 4
            meanArrivalNum = sum(sum(cubeG .* tShift3, 3), 4);
            meanSqNum = sum(sum(cubeG .* tShiftSq3, 3), 4);
            denom = max(sum(sum(cubeG, 3), 4), 1);
        else
            meanArrivalNum = sum(cubeG .* tShift3, 3);
            meanSqNum = sum(cubeG .* tShiftSq3, 3);
            denom = max(sum(cubeG, 3), 1);
        end
        meanArrival = gather(meanArrivalNum ./ denom);
        meanSq = gather(meanSqNum ./ denom);
    else
        tShift3 = reshape(tShift, 1, 1, []);
        tShiftSq3 = reshape(tShift.^2, 1, 1, []);
        if ndims(cubeG) == 4
            meanArrivalNum = sum(sum(cubeG .* tShift3, 3), 4);
            meanSqNum = sum(sum(cubeG .* tShiftSq3, 3), 4);
            denom = max(sum(sum(cubeG, 3), 4), 1);
        else
            meanArrivalNum = sum(cubeG .* tShift3, 3);
            meanSqNum = sum(cubeG .* tShiftSq3, 3);
            denom = max(sum(cubeG, 3), 1);
        end
        meanArrival = meanArrivalNum ./ denom;
        meanSq = meanSqNum ./ denom;
    end

    tauRMS = sqrt(max(meanSq - meanArrival.^2, 0));

    flim = struct();
    flim.total = struct();
    flim.total.tag = double(tag);
    flim.total.tauMean = double(meanArrival);
    flim.total.tauRMS = double(tauRMS);
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
            [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
        return;
    end

    if isfield(ptuOut, 'tcspc_pix') && ~isempty(ptuOut.tcspc_pix)
        cube = ptuOut.tcspc_pix;
        dtNs = ptuOut.head.MeasDesc_Resolution * 1e9;
        counts = squeeze(sum(sum(sum(cube, 4), 2), 1));
        counts = double(counts(:));
        srcInfo.mode = 'tcspc cube';
        srcInfo.native = false;
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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
        [counts, srcInfo] = trimTrailingTcspcSupport(counts, srcInfo);
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

    function [countsTrim, srcInfoTrim] = trimTrailingTcspcSupport(counts, srcInfo)
        countsTrim = double(counts(:));
        srcInfoTrim = srcInfo;
        if isempty(countsTrim)
            return;
        end

        lastIdx = find(countsTrim > 0, 1, 'last');
        if isempty(lastIdx)
            lastIdx = numel(countsTrim);
        end
        lastIdx = max(1, min(numel(countsTrim), lastIdx));
        if lastIdx >= numel(countsTrim)
            return;
        end

        countsTrim = countsTrim(1:lastIdx);
        if isstruct(srcInfoTrim)
            if isfield(srcInfoTrim, 'tAxisNs') && ~isempty(srcInfoTrim.tAxisNs)
                srcInfoTrim.tAxisNs = srcInfoTrim.tAxisNs(1:min(lastIdx, numel(srcInfoTrim.tAxisNs)));
            end
            if isfield(srcInfoTrim, 'widthNs') && ~isempty(srcInfoTrim.widthNs)
                srcInfoTrim.widthNs = srcInfoTrim.widthNs(1:min(lastIdx, numel(srcInfoTrim.widthNs)));
            end
            srcInfoTrim.trimmedTrailingBins = numel(counts) - lastIdx;
        end
    end

    function [countsShift, shiftBins, shiftNs, peakIdx, riseIdx] = shiftDecayForIRFDisplay(counts, dtNs)
% Preserve the full measured TCSPC ordering for display.
% Rise/peak indices are still reported for annotations and fitting metadata.
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
    shiftBins = 0;
    shiftNs = 0;
    countsShift = counts;
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

    [coeff, fitCounts] = solveDecayCoefficientsPIRLS(M, counts);
    err = sum((counts - fitCounts).^2 ./ max(fitCounts, 1));
end

function [coeff, fitCounts] = solveDecayCoefficientsPIRLS(M, counts)
    counts = max(double(counts(:)), 0);
    M = max(double(M), 0);
    Mscaled = normalizeDecayModelColumns(M);
    coeff = solveNonnegativePIRLS(Mscaled, counts, 10);
    coeff = max(double(coeff(:)), 0);
    fitCounts = Mscaled * coeff;
    coeff = coeff(:).';
end

function Mscaled = normalizeDecayModelColumns(M)
    M = double(M);
    colScale = sum(M, 1);
    colScale(~isfinite(colScale) | colScale <= 0) = 1;
    Mscaled = M ./ colScale;
end

function coeff = solveNonnegativePIRLS(M, y, maxNumIter)
    if nargin < 3 || isempty(maxNumIter)
        maxNumIter = 10;
    end
    try
        coeff = double(PIRLSnonneg(M, y, maxNumIter));
        coeff = coeff(:);
        return;
    catch
    end
    coeff = localProjectedGradientPIRLSVec(M, y, maxNumIter, 50);
end

function beta = localProjectedGradientPIRLSVec(M, y, maxPirlsIter, maxPgIter)
    M = double(M);
    y = double(y(:));
    [~, nBasis] = size(M);
    tiny = max(1e-6, 0.1 / max(1, size(M, 1)));
    reg = 1e-6 * eye(nBasis);
    beta = max((M' * M + reg) \ (M' * y), 0);

    for iterIdx = 1:maxPirlsIter
        mu = max(M * beta, tiny);
        w = 1 ./ mu;
        Aw = M' * (M .* w) + reg;
        bw = M' * (w .* y);
        betaNew = projectedGradientNnlsVec(Aw, bw, maxPgIter, beta);
        if sum((betaNew - beta).^2) < 1e-10
            beta = betaNew;
            break;
        end
        beta = betaNew;
    end
end

function beta = projectedGradientNnlsVec(A, b, maxIter, beta)
    if nargin < 4 || isempty(beta)
        beta = max(A \ b, 0);
    else
        beta = max(double(beta(:)), 0);
    end
    L = max(sum(abs(A), 2));
    if ~isfinite(L) || L <= 0
        L = 1;
    end
    for iterIdx = 1:maxIter
        grad = A * beta - b;
        betaNew = max(beta - grad ./ L, 0);
        if sum((betaNew - beta).^2) < 1e-10
            beta = betaNew;
            break;
        end
        beta = betaNew;
    end
end

function bgPerBin = fittedBackgroundOffsetPerBin(coeff, nBins, includeBG)
    bgPerBin = 0;
    if nargin < 3 || ~includeBG || isempty(coeff) || isempty(nBins) || nBins < 1
        return;
    end
    coeff = double(coeff(:));
    bgPerBin = max(coeff(1), 0) / max(double(nBins), 1);
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
    [~, yhat] = solveDecayCoefficientsPIRLS(model, y);
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
    [~, yhat] = solveDecayCoefficientsPIRLS(model, y);
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
