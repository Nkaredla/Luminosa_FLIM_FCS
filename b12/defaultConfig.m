function cfg = defaultConfig(cfg)
%DEFAULTCONFIG Fill missing fields in the main user configuration.
% The configuration structure is permitted only at the public entry point.

    defaults = struct();
    defaults.cnum = 3;
    defaults.mitoPulseIndex = 3;
    defaults.mitoDetectorIndex = 2;
    defaults.b12PulseIndex = 2;
    defaults.b12DetectorIndex = 1;
    defaults.tcspcBinNs = 0.128;
    defaults.maxNgate = 1024;
    defaults.photonsPerChunk = 5e6;
    defaults.showWaitbar = false;
    defaults.useGPU = false;

    defaults.gateThresholdFrac = 0.15;
    defaults.gatePreNs = 1.5;
    defaults.minGateSeparationBins = 50;
    defaults.secondPeakMinFraction = 0.05;

    defaults.cellSmoothSigma = 1.5;
    % Multiplies the B12 Otsu threshold. Higher values make the detected
    % cell boundary smaller and stricter; lower values include dimmer edges.
    defaults.cellThresholdScale = 0.85;
    defaults.cellCloseRadius = 5;
    defaults.cellDilateRadius = 2;
    defaults.minCellArea = 150;

    defaults.mitoSmoothSigma = 1.0;
    defaults.mitoBackgroundRadius = 7;
    defaults.mitoThresholdScale = 0.85;
    defaults.mitoCloseRadius = 1;
    defaults.mitoDilateRadius = 0;
    defaults.minMitoArea = 3;

    % Coarse B12 bright-spot segmentation. Strong smoothing and morphology
    % intentionally merge tiny fragments into one lysosome-scale region.
    defaults.lysoSigmaSmall = 2.0;       % Higher: smoother/coarser spots.
    defaults.lysoSigmaLarge = 6.0;       % Local-background length scale.
    defaults.lysoThresholdMAD = 2.5;     % Higher: fewer contrast spots.
    defaults.lysoMinPeakFraction = 0.08; % Higher: require stronger contrast.
    defaults.lysoBrightQuantile = 0.92;  % Higher: retain only brighter spots.
    defaults.lysoExcludeMitoRadius = 0;
    defaults.lysoCloseRadius = 2;        % Higher: merge nearby spot pixels.
    defaults.lysoDilateRadius = 2;       % Higher: enlarge each coarse region.
    defaults.minLysoArea = 12;           % Higher: reject more small regions.

    defaults.irf = [];
    defaults.irfTau0 = [0.4 2.5];
    defaults.irfNCasc = 4;
    defaults.irfNSub = 6;
    defaults.minPhotonsForIRF = 500;

    defaults.maxExp = 2;
    defaults.tauSeeds1 = 1.5;
    defaults.tauSeeds2 = [0.4 2.5];
    defaults.tauBoundsNs = [0.05 10];
    defaults.includeBackground = true;
    defaults.minPhotonsPerRegionFit = 100;
    defaults.fitMaxIter = 500;
    defaults.bicImprovementForTwoExp = 2.0;
    defaults.minStateSeparationFraction = 0.10;
    defaults.minStateAmplitudeFraction = 0.02;

    defaults.bayesBatchSize = 2048;
    defaults.bayesMinPixelPhotons = 3;
    defaults.bayesMinRegionPhotons = 20;
    defaults.bayesSignalGrid = linspace(0, 1, 26);
    defaults.bayesFractionGrid = linspace(0, 1, 41);
    defaults.bayesSingleTauGrid = [];
    defaults.bayesShiftBounds = [-5 5];

    defaults.showFigures = true;
    defaults.saveOutputs = true;
    defaults.outputDir = '';
    defaults.outputStem = 'B12_MitoLysoCytosol_BayesFLIM';
    defaults.maxFitPlots = 9;
    defaults.lifetimeHistogramBins = 20;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k}) || isempty(cfg.(names{k}))
            cfg.(names{k}) = defaults.(names{k});
        end
    end
end
