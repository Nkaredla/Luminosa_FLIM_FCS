rootFolder = 'D:\Luminosa\Data\Natasha';
outputFolder = []; % root-level batch summary folder is chosen automatically

opts = struct();
opts.irfMode = 'global';
% opts.irfMode = 'supplied';
% opts.irfFile = 'D:\Luminosa\Data\Natasha\IRF.ptu';
opts.globalIrfMethod = 'calc_mirf';
opts.tau0 = []; % let DistFluofit choose the required component seeds
opts.init = 1;
opts.fluofitSolver = 'mle'; % use 'pirls' for PIRLSnonneg amplitude fitting
opts.plotFits = true;
opts.saveFitFigures = true;
opts.figureFormats = {'png', 'fig'};
opts.searchRecursively = true;
opts.multipleFileMode = 'largest_ptu'; % choose largest PTU when several nested files are found
opts.resultFolderName = 'tcspc_fluofit_results';

batch = fit_tcspc_subfolders_with_fluofit(rootFolder, outputFolder, opts);
