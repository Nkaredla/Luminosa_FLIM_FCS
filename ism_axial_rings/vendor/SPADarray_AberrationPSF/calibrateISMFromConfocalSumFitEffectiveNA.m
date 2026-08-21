function calibration = calibrateISMFromConfocalSumFitEffectiveNA(stackInput, varargin)
%CALIBRATEISMFROMCONFOCALSUMFITEFFECTIVENA Fit effective NA by profiling.
%
%   calibration = calibrateISMFromConfocalSumFitEffectiveNA(stackInput, ...)
%
%   This is a deliberately simple, robust outer loop around
%   calibrateISMFromConfocalSum.  The ordinary 'objectiveNA' argument is an
%   input to the forward model; it is not, by itself, a fitted parameter. This
%   function makes the effective NA concrete by:
%
%       1) running a reduced-complexity parametric Stage-A fit over a bounded
%          NA grid,
%       2) refining the grid around the best value,
%       3) rerunning the final aberration/collection-efficiency calibration at
%          the selected NA.
%
%   The NA choice is based on Stage-A confocal 3-D PSF residual, not on the
%   detector collection-efficiency map.  This avoids letting detector QE /
%   crosstalk absorb an optical-NA error.
%
%   Common options handled by this wrapper:
%       'objectiveNA'              nominal starting NA, default 1.2
%       'effectiveNABounds'        [min max], default objectiveNA*[0.85 1.10]
%       'effectiveNACoarseStep'    default 0.03
%       'effectiveNARefineStep'    default 0.01
%       'coarseFitModes'           reduced modes for the NA profile
%       'finalFitModes'            modes for final run; [] means pass-through
%       'writeNAProfileOutputs'    save candidate figures, default false
%       'outputDir'                final output directory; candidates go under
%                                  outputDir/effective_na_profile
%
%   All other name/value pairs are passed through to calibrateISMFromConfocalSum.

    if nargin < 1 || isempty(stackInput)
        error('calibrateISMFromConfocalSumFitEffectiveNA:MissingInput', ...
            'stackInput is required.');
    end

    [opts, passArgs] = parseEffectiveNAOptions(varargin{:});
    opts = completeEffectiveNAOptions(opts, stackInput);

    if opts.verbose
        fprintf('\nEffective-NA profile calibration\n');
        fprintf('  nominal input NA %.4g; search bounds [%.4g %.4g]\n', ...
            opts.objectiveNA, opts.effectiveNABounds(1), opts.effectiveNABounds(2));
        fprintf('  coarse step %.4g; refine step %.4g\n', ...
            opts.effectiveNACoarseStep, opts.effectiveNARefineStep);
    end

    % Load the detector stack ONCE and reuse the in-memory struct for every NA
    % candidate and the final run. Without this, each iteration re-entered
    % loadFullStackISMData and re-read every PTU frame (and the dark PTU).
    stackData = loadStackOnce(stackInput, passArgs, opts);

    coarseNA = boundedGrid(opts.effectiveNABounds(1), opts.effectiveNABounds(2), ...
        opts.effectiveNACoarseStep, opts.objectiveNA);
    coarse = runNAProfile(stackData, passArgs, opts, coarseNA, ...
        opts.coarseFitModes, 'coarse');

    bestCoarseNA = bestProfileNA(coarse);
    refineBounds = [ ...
        max(opts.effectiveNABounds(1), bestCoarseNA - opts.effectiveNACoarseStep), ...
        min(opts.effectiveNABounds(2), bestCoarseNA + opts.effectiveNACoarseStep)];
    refineNA = boundedGrid(refineBounds(1), refineBounds(2), ...
        opts.effectiveNARefineStep, bestCoarseNA);
    refineNA = setdiff(refineNA, coarse.objectiveNA(:).', 'stable');
    if isempty(refineNA)
        refined = emptyProfileTable();
    else
        refined = runNAProfile(stackData, passArgs, opts, refineNA, ...
            opts.coarseFitModes, 'refine');
    end

    profile = [coarse; refined]; %#ok<AGROW>
    profile = sortrows(profile, 'objectiveNA');
    fittedNA = bestProfileNA(profile);

    if opts.verbose
        fprintf('  selected effective NA %.5g (Stage-A rel.residual %.4g)\n', ...
            fittedNA, min(profile.relResidual));
    end

    finalArgs = passArgs;
    finalArgs = setNameValue(finalArgs, 'objectiveNA', fittedNA);
    finalArgs = setNameValue(finalArgs, 'wavefrontMethod', 'parametric');
    finalArgs = setNameValue(finalArgs, 'outputDir', opts.outputDir);
    finalArgs = setNameValue(finalArgs, 'writeOutputs', opts.writeOutputs);
    finalArgs = setNameValue(finalArgs, 'verbose', opts.verbose);
    if ~isempty(opts.finalFitModes)
        finalArgs = setNameValue(finalArgs, 'fitModes', opts.finalFitModes);
    end

    calibration = calibrateISMFromConfocalSum(stackData, finalArgs{:});
    calibration.effectiveNAFit = struct();
    calibration.effectiveNAFit.method = 'bounded residual profile';
    calibration.effectiveNAFit.nominalObjectiveNA = opts.objectiveNA;
    calibration.effectiveNAFit.fittedObjectiveNA = fittedNA;
    calibration.effectiveNAFit.bounds = opts.effectiveNABounds;
    calibration.effectiveNAFit.profile = profile;
    calibration.effectiveNAFit.selectionMetric = ...
        'Stage-A confocal 3-D PSF relative residual';
    calibration.optics.effectiveObjectiveNA = fittedNA;
    calibration.centerCoordinateFitDefaults.objectiveNA = fittedNA;

    if opts.writeOutputs
        save(fullfile(opts.outputDir, 'confocal_sum_calibration.mat'), ...
            'calibration', '-v7.3');
        writetable(profile, fullfile(opts.outputDir, ...
            'stageA_effective_na_profile.csv'));
    end
end

function profile = runNAProfile(stackData, passArgs, opts, candidates, fitModes, tag)
    profile = emptyProfileTable();
    for k = 1:numel(candidates)
        na = candidates(k);
        runArgs = passArgs;
        runArgs = setNameValue(runArgs, 'objectiveNA', na);
        runArgs = setNameValue(runArgs, 'wavefrontMethod', 'parametric');
        runArgs = setNameValue(runArgs, 'fitModes', fitModes);
        runArgs = setNameValue(runArgs, 'writeOutputs', opts.writeNAProfileOutputs);
        runArgs = setNameValue(runArgs, 'verbose', opts.verbose);
        runArgs = setNameValue(runArgs, 'outputDir', fullfile( ...
            opts.profileOutputDir, sprintf('%s_NA_%0.4f', tag, na)));

        if opts.verbose
            fprintf('  [%s] NA candidate %d/%d: %.5g\n', ...
                tag, k, numel(candidates), na);
        end

        try
            cal = calibrateISMFromConfocalSum(stackData, runArgs{:});
            relResidual = scalarFieldOrNaN(cal.wavefront, 'relResidual');
            z0Um = scalarFieldOrNaN(cal.wavefront, 'z0Um');
            converged = true;
            message = "";
        catch err
            relResidual = Inf;
            z0Um = NaN;
            converged = false;
            message = string(err.message);
            if opts.verbose
                fprintf('      failed: %s\n', err.message);
            end
        end

        row = table(na, relResidual, z0Um, converged, message, ...
            'VariableNames', {'objectiveNA','relResidual','z0Um', ...
            'converged','message'});
        profile = [profile; row]; %#ok<AGROW>
    end
end

function data = loadStackOnce(stackInput, passArgs, opts)
    % Load the detector stack a single time. If the caller already passed a
    % loaded data struct, reuse it as-is; otherwise read it once with the
    % loader-relevant options pulled from the pass-through arguments.
    if isstruct(stackInput) && isfield(stackInput, 'rawCounts')
        data = stackInput;
        return;
    end
    loaderNames = {'channelIDs','darkFile','darkScale','backgroundMode', ...
        'xyPixelSizeUm','reuseCache','channelOrder','darkPerPixel', ...
        'darkMeasurementMode','flatField','flatFieldVariable','boundaryWidthPx', ...
        'ptuReaderFolder','ptuPhotonsPerChunk','cacheFile'};
    loaderArgs = {};
    for k = 1:2:numel(passArgs)
        if any(strcmpi(passArgs{k}, loaderNames))
            loaderArgs(end+1:end+2) = passArgs(k:k+1); %#ok<AGROW>
        end
    end
    if opts.verbose
        fprintf('  loading detector stack once for the whole NA profile...\n');
    end
    data = loadFullStackISMData(stackInput, loaderArgs{:}, 'verbose', opts.verbose);
end

function [opts, passArgs] = parseEffectiveNAOptions(varargin)
    opts = struct();
    opts.objectiveNA = 1.2;
    opts.effectiveNABounds = [];
    opts.effectiveNACoarseStep = 0.03;
    opts.effectiveNARefineStep = 0.01;
    opts.coarseFitModes = {'astig_x','astig_y','coma_x','coma_y','spherical'};
    opts.finalFitModes = [];
    opts.writeNAProfileOutputs = false;
    opts.outputDir = '';
    opts.writeOutputs = true;
    opts.verbose = true;

    passArgs = {};
    if mod(numel(varargin), 2) ~= 0
        error('calibrateISMFromConfocalSumFitEffectiveNA:BadNameValue', ...
            'Optional arguments must be Name/Value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = char(varargin{k});
        value = varargin{k+1};
        switch lower(name)
            case 'objectivena'
                opts.objectiveNA = double(value);
            case 'effectivenabounds'
                opts.effectiveNABounds = double(value);
            case 'effectivenacoarsestep'
                opts.effectiveNACoarseStep = double(value);
            case 'effectivenarefinestep'
                opts.effectiveNARefineStep = double(value);
            case 'coarsefitmodes'
                opts.coarseFitModes = normalizeModeList(value);
            case {'finalfitmodes','fitmodes'}
                opts.finalFitModes = normalizeModeList(value);
            case 'writenaprofileoutputs'
                opts.writeNAProfileOutputs = logical(value);
            case 'outputdir'
                opts.outputDir = char(value);
            case 'writeoutputs'
                opts.writeOutputs = logical(value);
            case 'verbose'
                opts.verbose = logical(value);
            case 'wavefrontmethod'
                if ~strcmpi(value, 'parametric')
                    warning('calibrateISMFromConfocalSumFitEffectiveNA:ForcingParametric', ...
                        ['Effective NA fitting is only defined for the ' ...
                        'parametric air-interface model; ignoring wavefrontMethod=%s.'], ...
                        char(value));
                end
            otherwise
                passArgs(end+1:end+2) = varargin(k:k+1); %#ok<AGROW>
        end
    end
end

function opts = completeEffectiveNAOptions(opts, stackInput)
    if ~isscalar(opts.objectiveNA) || ~isfinite(opts.objectiveNA) || ...
            opts.objectiveNA <= 0
        error('calibrateISMFromConfocalSumFitEffectiveNA:BadNA', ...
            'objectiveNA must be a finite positive scalar.');
    end
    if isempty(opts.effectiveNABounds)
        opts.effectiveNABounds = opts.objectiveNA * [0.85 1.10];
    end
    opts.effectiveNABounds = double(opts.effectiveNABounds(:)).';
    if numel(opts.effectiveNABounds) ~= 2 || ...
            any(~isfinite(opts.effectiveNABounds)) || ...
            opts.effectiveNABounds(1) <= 0 || ...
            opts.effectiveNABounds(2) <= opts.effectiveNABounds(1)
        error('calibrateISMFromConfocalSumFitEffectiveNA:BadBounds', ...
            'effectiveNABounds must be [min max] with 0 < min < max.');
    end
    opts.effectiveNACoarseStep = validatePositiveScalar( ...
        opts.effectiveNACoarseStep, 'effectiveNACoarseStep');
    opts.effectiveNARefineStep = validatePositiveScalar( ...
        opts.effectiveNARefineStep, 'effectiveNARefineStep');
    if isempty(opts.outputDir)
        opts.outputDir = defaultEffectiveNAOutputDir(stackInput);
    end
    opts.profileOutputDir = fullfile(opts.outputDir, 'effective_na_profile');
    if opts.writeOutputs
        if exist(opts.outputDir, 'dir') ~= 7
            mkdir(opts.outputDir);
        end
        if exist(opts.profileOutputDir, 'dir') ~= 7
            mkdir(opts.profileOutputDir);
        end
    end
end

function modes = normalizeModeList(value)
    if isempty(value)
        modes = {};
    elseif ischar(value) || isstring(value)
        modes = cellstr(value);
    else
        modes = value;
    end
    modes = modes(:).';
end

function value = validatePositiveScalar(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
        error('calibrateISMFromConfocalSumFitEffectiveNA:BadPositiveScalar', ...
            '%s must be a finite positive scalar.', name);
    end
    value = double(value);
end

function grid = boundedGrid(lo, hi, step, includeValue)
    grid = lo:step:hi;
    grid = unique([grid hi includeValue]);
    grid = grid(isfinite(grid) & grid >= lo & grid <= hi);
    grid = round(grid, 6);
    grid = unique(grid, 'stable');
end

function na = bestProfileNA(profile)
    valid = isfinite(profile.relResidual);
    if ~any(valid)
        error('calibrateISMFromConfocalSumFitEffectiveNA:NoValidNA', ...
            'All effective-NA candidates failed.');
    end
    idx = find(valid);
    [~, j] = min(profile.relResidual(idx));
    na = profile.objectiveNA(idx(j));
end

function x = scalarFieldOrNaN(s, fieldName)
    x = NaN;
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
        if isnumeric(value) && isscalar(value) && isfinite(value)
            x = double(value);
        end
    end
end

function args = setNameValue(args, name, value)
    idx = find(strcmpi(args(1:2:end), name), 1, 'first');
    if isempty(idx)
        args(end+1:end+2) = {name, value};
    else
        args{2*idx} = value;
    end
end

function T = emptyProfileTable()
    T = table(zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), strings(0,1), ...
        'VariableNames', {'objectiveNA','relResidual','z0Um', ...
        'converged','message'});
end

function outDir = defaultEffectiveNAOutputDir(stackInput)
    root = fileparts(fileparts(mfilename('fullpath')));
    [~, stem] = fileparts(char(stackInput));
    stem = regexprep(stem, '_frame_alignment$', '');
    outDir = fullfile(root, 'output_matlab', ...
        'confocal_sum_calibration_effective_na', sanitizeName(stem));
end

function name = sanitizeName(name)
    name = regexprep(char(name), '[^A-Za-z0-9_\-]+', '_');
    if isempty(name)
        name = 'calibration';
    end
end
