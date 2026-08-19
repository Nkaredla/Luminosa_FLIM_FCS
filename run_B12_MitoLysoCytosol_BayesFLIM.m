% Example three-laser B12 compartment analysis.
ptuFile = 'D:\B12\010_B12Rho_24h_++_P15_20260807_12\FLIM_20260807-213453\RawImage.ptu';

cfg = struct();
cfg.cnum = 3;
cfg.mitoPulseIndex = 3;       % 640 nm
cfg.mitoDetectorIndex = 2;    % Det2
cfg.b12PulseIndex = 2;        % 560 nm
cfg.b12DetectorIndex = 1;     % Det1

% Cell segmentation uses only the B12 channel.
% Higher threshold scale gives a smaller/stricter boundary; lower includes
% dimmer cell edges. The final boundary is filled without internal holes.
cfg.cellThresholdScale = 0.85;
cfg.cellSmoothSigma = 2;
cfg.cellCloseRadius = 2;
cfg.cellDilateRadius = 5;
cfg.minCellArea = 650;

% Mitochondrial segmentation from 640 nm / Det2.
cfg.mitoThresholdScale = 0.85;

% Coarse lysosome segmentation from bright B12 spots.
cfg.lysoSigmaSmall = 2.0;       % Higher -> smoother, coarser regions.
cfg.lysoSigmaLarge = 6.0;       % Background-removal length scale.
cfg.lysoThresholdMAD = 2.5;     % Higher -> fewer contrast detections.
cfg.lysoMinPeakFraction = 0.08; % Higher -> require stronger local contrast.
cfg.lysoBrightQuantile = 0.92;  % Higher -> retain only brighter B12 spots.
cfg.lysoCloseRadius = 2;        % Higher -> merge neighbouring spot pixels.
cfg.lysoDilateRadius = 2;       % Higher -> enlarge coarse regions.
cfg.minLysoArea = 12;           % Higher -> remove more small regions.
cfg.lysoExcludeMitoRadius = 0;  % Increase to widen mitochondrial exclusion.

% Number of common bins in the mitochondrial/lysosome lifetime histograms.
cfg.lifetimeHistogramBins = 20;

cfg.showFigures = true;
cfg.saveOutputs = true;
result = PTU_MitoLysoCytosol_BayesFLIM(ptuFile, cfg);
