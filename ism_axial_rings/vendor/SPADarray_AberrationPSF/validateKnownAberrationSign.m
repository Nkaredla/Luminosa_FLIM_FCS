function report = validateKnownAberrationSign(baselineInput, inducedInput, commandedCoeffs, varargin)
%VALIDATEKNOWNABERRATIONSIGN Compare recovered changes with known commands.
%
%   report = validateKnownAberrationSign(baselineResult, inducedResult, ...
%       struct('coma_x', 0.10))
%
%   baselineInput and inducedInput can be result structures or MAT files
%   containing a result structure from estimateCenterPointISMWavefront or
%   phaseRetrieval3DBead. commandedCoeffs contains the known applied change
%   in waves RMS using the same named Zernike convention.

    p = inputParser;
    addParameter(p, 'outputCsv', '');
    addParameter(p, 'minimumCommandWaves', 1e-4);
    parse(p, varargin{:});
    opts = p.Results;

    if ~isstruct(commandedCoeffs) || isempty(fieldnames(commandedCoeffs))
        error('validateKnownAberrationSign:BadCommand', ...
            'commandedCoeffs must be a nonempty structure of waves-RMS changes.');
    end

    baseline = loadResultInput(baselineInput);
    induced = loadResultInput(inducedInput);
    baselineCoeffs = extractCoefficientStruct(baseline);
    inducedCoeffs = extractCoefficientStruct(induced);

    modes = fieldnames(commandedCoeffs);
    nModes = numel(modes);
    commandWaves = nan(nModes,1);
    baselineWaves = nan(nModes,1);
    inducedWaves = nan(nModes,1);
    recoveredDeltaWaves = nan(nModes,1);
    signAgreement = false(nModes,1);
    recoveredPerCommand = nan(nModes,1);

    for k = 1:nModes
        mode = modes{k};
        commandWaves(k) = double(commandedCoeffs.(mode));
        baselineWaves(k) = coefficientValue(baselineCoeffs, mode);
        inducedWaves(k) = coefficientValue(inducedCoeffs, mode);
        recoveredDeltaWaves(k) = inducedWaves(k) - baselineWaves(k);
        if abs(commandWaves(k)) >= opts.minimumCommandWaves
            signAgreement(k) = sign(recoveredDeltaWaves(k)) == sign(commandWaves(k));
            recoveredPerCommand(k) = recoveredDeltaWaves(k) / commandWaves(k);
        end
    end

    mode = modes(:);
    tableResult = table(mode, commandWaves, baselineWaves, inducedWaves, ...
        recoveredDeltaWaves, recoveredPerCommand, signAgreement);
    tested = abs(commandWaves) >= opts.minimumCommandWaves;

    report = struct();
    report.table = tableResult;
    report.allTestedSignsAgree = all(signAgreement(tested));
    report.nModesTested = nnz(tested);
    report.interpretation = ['A positive recoveredPerCommand confirms the ' ...
        'absolute instrument-to-model sign. Its magnitude is the recovered gain.'];

    if ~isempty(opts.outputCsv)
        writetable(tableResult, opts.outputCsv);
    end

    fprintf('\nKnown-aberration sign calibration\n');
    disp(tableResult);
    fprintf('  all tested signs agree: %d\n\n', report.allTestedSignsAgree);
end

function result = loadResultInput(inputValue)
    if isstruct(inputValue)
        result = inputValue;
        return;
    end
    if ~(ischar(inputValue) || isstring(inputValue))
        error('validateKnownAberrationSign:BadInput', ...
            'Fit inputs must be result structures or MAT-file paths.');
    end
    loaded = load(char(inputValue));
    names = fieldnames(loaded);
    preferred = {'result','res','fitResult','fit'};
    for k = 1:numel(preferred)
        if isfield(loaded, preferred{k}) && isstruct(loaded.(preferred{k}))
            result = loaded.(preferred{k});
            return;
        end
    end
    for k = 1:numel(names)
        if isstruct(loaded.(names{k}))
            result = loaded.(names{k});
            return;
        end
    end
    error('validateKnownAberrationSign:NoResult', ...
        'No result structure was found in %s.', char(inputValue));
end

function coeffs = extractCoefficientStruct(result)
    if isfield(result, 'fit') && isstruct(result.fit)
        coeffs = extractCoefficientStruct(result.fit);
        return;
    end
    if isfield(result, 'estCoeffs') && isstruct(result.estCoeffs)
        coeffs = result.estCoeffs;
        return;
    end
    error('validateKnownAberrationSign:MissingCoefficients', ...
        'The input result does not contain fit.estCoeffs or estCoeffs.');
end

function value = coefficientValue(coeffs, mode)
    value = 0;
    if isfield(coeffs, mode)
        value = double(coeffs.(mode));
    end
end
