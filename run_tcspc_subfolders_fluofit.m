rootFolder = 'E:\AAA_Chaofan_Rho';
outputFolder = []; % root-level batch summary folder is chosen automatically

opts = struct();
opts.fitMethod = 'tail'; % fast tail fitting: no IRF reconvolution
opts.tailCutAfterPeakNs = 0.3; % start fitting after the TCSPC peak
opts.tailRejectLastNPoints = 20; % remove final tail bins before Tailfit/DistTailfit
opts.tailSolver = 'pirls';
opts.irfMode = 'per_curve'; % ignored when opts.fitMethod='tail'
% opts.irfMode = 'best_per_curve'; % pick the per-file IRF with the best biexponential chi2
% opts.irfMode = 'per_curve'; % estimate the IRF separately for each selected file
% opts.irfMode = 'supplied';
% opts.irfFile = 'D:\Luminosa\Data\Natasha\IRF.ptu';
opts.perCurveIrfModel = 'spad_exgauss'; % SPAD: Gaussian prompt peak + exponential tail
opts.globalIrfMethod = 'spad_exgauss';
opts.irfClipFraction = 1e-3;
opts.tau0 = [0.5 2.0]; % initial guess for direct biexponential Tailfit
opts.init = 0; % tail mode: 0 skips DistTailfit; 1 also runs DistTailfit
opts.fluofitSolver = 'pirls'; % 'mle' use 'pirls' for PIRLSnonneg amplitude fitting
opts.bestIrfTau0 = [0.4 2.0]; % biexponential screening fit used only to choose the IRF
opts.bestIrfSolver = opts.fluofitSolver;
opts.plotFits = true;
opts.saveFitFigures = true;
opts.figureFormats = {'png', 'fig'};
opts.copyFitsToBatchFolder = true;
opts.summaryAmplitudeThreshold = 0.02; % hide components below 2% in the summary sheet
opts.searchRecursively = true;
opts.multipleFileMode = 'largest_ptu'; % choose largest PTU when several nested files are found
opts.resultFolderName = 'tcspc_fluofit_results';
opts.useTcspcCache = true;
opts.forceReadTcspc = false; % set true to rebuild tcspc_curves_cache.mat

close all
clc
batch = fit_tcspc_subfolders_with_fluofit(rootFolder, outputFolder, opts);
