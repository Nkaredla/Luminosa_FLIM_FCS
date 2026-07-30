rootFolder = 'D:\Luminosa\Data\Natasha';
outputFolder = fullfile(rootFolder, 'tcspc_batch_fluofit_results');

opts = struct();
opts.irfMode = 'global';
% opts.irfMode = 'supplied';
% opts.irfFile = 'D:\Luminosa\Data\Natasha\IRF.ptu';
opts.globalIrfMethod = 'calc_mirf';
opts.fluofitSolver = 'mle'; % use 'pirls' for PIRLSnonneg amplitude fitting
opts.plotFits = false;

batch = fit_tcspc_subfolders_with_fluofit(rootFolder, outputFolder, opts);
