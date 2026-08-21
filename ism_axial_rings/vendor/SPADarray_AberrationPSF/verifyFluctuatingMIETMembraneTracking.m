function results = verifyFluctuatingMIETMembraneTracking(varargin)
%VERIFYFLUCTUATINGMIETMEMBRANETRACKING Numerical checks for the MIET model.
%
%   results = verifyFluctuatingMIETMembraneTracking()
%   results = verifyFluctuatingMIETMembraneTracking('runEndToEnd',true)
%
%   The default checks are independent of the simulator implementation:
%     1. closed-form height-graph Ito drift versus a numerical divergence;
%     2. lateral-error attenuation and its moment correction;
%     3. noncentral Gaussian mass in an anisotropic tracking ellipse;
%     4. naive versus gate-conditioned diffusion MLE.
%
%   runEndToEnd=true also runs a short photon-to-topography simulation.

    p = inputParser;
    p.FunctionName = 'verifyFluctuatingMIETMembraneTracking';
    addParameter(p,'rngSeed',11);
    addParameter(p,'nCurvatureReplicates',300);
    addParameter(p,'nGateDraws',250000);
    addParameter(p,'runEndToEnd',false);
    addParameter(p,'verbose',true);
    parse(p,varargin{:});
    options = p.Results;
    rng(options.rngSeed);

    drift = verifyDriftFormula();
    curvature = verifyMomentCorrection(options.nCurvatureReplicates);
    gate = verifyGateConditioning(options.nGateDraws);
    endToEnd = struct('wasRun',false);
    if options.runEndToEnd
        endToEnd.wasRun = true;
        endToEnd.output = simulateFluctuatingMIETMembraneTracking( ...
            'rngSeed',options.rngSeed,'nFrames',350, ...
            'maxGpObservations',180,'topographyGridSize',21, ...
            'nTopographySnapshots',2,'curvatureProfilePoints',9, ...
            'makeFigure',false,'writeOutputs',false,'verbose',false);
        endToEnd.topographyRmseUm = ...
            endToEnd.output.surfaceModel.topography.rmseUm;
        endToEnd.positionRmseUm = ...
            endToEnd.output.diffusion.truthValidation.lateralLocalizationRmseUm;
        endToEnd.curvaturePerUm = ...
            endToEnd.output.surfaceModel.curvaturePerUm;
        endToEnd.diffusionUm2PerS = ...
            endToEnd.output.diffusion.DFisherMLEUm2PerS;
    end

    results = struct();
    results.options = options;
    results.drift = drift;
    results.curvature = curvature;
    results.gate = gate;
    results.endToEnd = endToEnd;
    results.passed = drift.passed && curvature.passed && gate.passed;

    if options.verbose
        printResults(results);
    end
end

% -------------------------------------------------------------------------
function result = verifyDriftFormula()
    coefficients = randn(4,4)*0.8;
    query = randn(4,2)*0.5;
    closed = zeros(4,2);
    numerical = zeros(4,2);
    relativeError = zeros(4,1);
    for test = 1:4
        [gradient,hessian] = numericGradientHessian( ...
            query(test,:),coefficients(test,:),1e-5);
        closed(test,:) = closedDrift(gradient,hessian);
        numerical(test,:) = numericalDrift( ...
            query(test,:),coefficients(test,:),2e-4);
        relativeError(test) = norm(closed(test,:)-numerical(test,:)) / ...
            max(norm(numerical(test,:)),1e-12);
    end
    result = struct();
    result.table = table((1:4).',vecnorm(closed,2,2), ...
        vecnorm(numerical,2,2),relativeError, ...
        'VariableNames',{'test','closedNorm','numericalNorm','relativeError'});
    result.maximumRelativeError = max(relativeError);
    result.passed = result.maximumRelativeError < 2e-4;
end

% -------------------------------------------------------------------------
function height = testSurface(position,coefficient)
    x = position(1);
    y = position(2);
    height = coefficient(1)*x^2 + coefficient(2)*y^2 + ...
        coefficient(3)*x*y + ...
        coefficient(4)*sin(2.1*x)*cos(1.7*y);
end

% -------------------------------------------------------------------------
function [gradient,hessian] = numericGradientHessian(position,coefficient,h)
    gradient = zeros(1,2);
    hessian = zeros(2);
    for firstAxis = 1:2
        firstOffset = zeros(1,2);
        firstOffset(firstAxis) = h;
        gradient(firstAxis) = (testSurface(position+firstOffset,coefficient)- ...
            testSurface(position-firstOffset,coefficient))/(2*h);
        for secondAxis = 1:2
            secondOffset = zeros(1,2);
            secondOffset(secondAxis) = h;
            hessian(firstAxis,secondAxis) = ( ...
                testSurface(position+firstOffset+secondOffset,coefficient)- ...
                testSurface(position+firstOffset-secondOffset,coefficient)- ...
                testSurface(position-firstOffset+secondOffset,coefficient)+ ...
                testSurface(position-firstOffset-secondOffset,coefficient)) / ...
                (4*h^2);
        end
    end
end

% -------------------------------------------------------------------------
function drift = closedDrift(gradient,hessian)
    p = gradient(:);
    metricScale = 1+p.'*p;
    drift = (-p*(trace(hessian)*metricScale-p.'*hessian*p) / ...
        metricScale^2).';
end

% -------------------------------------------------------------------------
function drift = numericalDrift(position,coefficient,h)
    drift = zeros(1,2);
    [~,~,rootDeterminant] = metricFactor(position,coefficient);
    for outputAxis = 1:2
        accumulated = 0;
        for derivativeAxis = 1:2
            offset = zeros(1,2);
            offset(derivativeAxis) = h;
            [factorPlus,~,~] = metricFactor( ...
                position+offset,coefficient);
            [factorMinus,~,~] = metricFactor( ...
                position-offset,coefficient);
            accumulated = accumulated + ...
                (factorPlus(outputAxis,derivativeAxis)- ...
                 factorMinus(outputAxis,derivativeAxis))/(2*h);
        end
        drift(outputAxis) = accumulated/rootDeterminant;
    end
end

% -------------------------------------------------------------------------
function [factor,inverseMetric,rootDeterminant] = ...
        metricFactor(position,coefficient)
    [gradient,~] = numericGradientHessian(position,coefficient,1e-5);
    p = gradient(:);
    metricScale = 1+p.'*p;
    inverseMetric = eye(2)-(p*p.')/metricScale;
    rootDeterminant = sqrt(metricScale);
    factor = rootDeterminant*inverseMetric;
end

% -------------------------------------------------------------------------
function result = verifyMomentCorrection(nReplicates)
    trueCurvature = -0.60;
    tipHeight = 0.18;
    heightSigma = 0.012;
    observationRadius = 0.28;
    n = 300;
    lateralSigma = [0.010;0.020;0.035;0.050];
    naiveMean = zeros(size(lateralSigma));
    correctedMean = zeros(size(lateralSigma));
    correctedSd = zeros(size(lateralSigma));
    for precision = 1:numel(lateralSigma)
        sigma = lateralSigma(precision);
        naive = zeros(nReplicates,1);
        corrected = zeros(nReplicates,1);
        for replicate = 1:nReplicates
            radius = observationRadius*sqrt(rand(n,1));
            angle = 2*pi*rand(n,1);
            truePosition = [radius.*cos(angle),radius.*sin(angle)];
            height = tipHeight+0.5*trueCurvature*radius.^2 + ...
                heightSigma*randn(n,1);
            measured = truePosition+sigma*randn(n,2);
            u = sum(measured.^2,2)-2*sigma^2;
            design = [ones(n,1),0.5*u];
            beta = design\height;
            regressorVariance = 4*sigma^2*max(u,0)+4*sigma^4;
            attenuation = max(1-mean(regressorVariance)/var(u),1e-3);
            naive(replicate) = beta(2);
            corrected(replicate) = beta(2)/attenuation;
        end
        naiveMean(precision) = mean(naive);
        correctedMean(precision) = mean(corrected);
        correctedSd(precision) = std(corrected);
    end
    naiveBiasPercent = 100*(naiveMean/trueCurvature-1);
    correctedBiasPercent = 100*(correctedMean/trueCurvature-1);
    result = struct();
    result.trueCurvaturePerUm = trueCurvature;
    result.table = table(1000*lateralSigma,naiveMean,naiveBiasPercent, ...
        correctedMean,correctedBiasPercent,correctedSd, ...
        'VariableNames',{'lateralSigmaNm','naiveCurvaturePerUm', ...
        'naiveBiasPercent','correctedCurvaturePerUm', ...
        'correctedBiasPercent','correctedEmpiricalSdPerUm'});
    result.maximumCorrectedBiasPercent = max(abs(correctedBiasPercent));
    result.passed = result.maximumCorrectedBiasPercent < 5;
end

% -------------------------------------------------------------------------
function result = verifyGateConditioning(nDraws)
    trueD = 0.20;
    dt = 1e-3;
    noiseCovariance = [4e-4,1.0e-4;1.0e-4,2.5e-4];
    covariance = noiseCovariance+2*trueD*dt*eye(2);
    driftTerm = [0.025,-0.015];
    meanStep = trueD*driftTerm;
    gateCovariance = [9e-4,1.5e-4;1.5e-4,6e-4];
    gateChi2 = 4;
    nAngles = 128;
    angles = ((0:nAngles-1)+0.5)*(2*pi/nAngles);
    direction = [cos(angles);sin(angles)];
    directionalPrecision = sum(direction.*(gateCovariance\direction),1);
    radialLimit = sqrt(gateChi2./directionalPrecision);

    lower = chol(covariance,'lower');
    step = meanStep+randn(nDraws,2)*lower.';
    accepted = sum((step/gateCovariance).*step,2) <= gateChi2;
    acceptedStep = step(accepted,:);
    empiricalAcceptance = mean(accepted);
    quadratureAcceptance = gateAcceptance(trueD,noiseCovariance,dt, ...
        driftTerm,radialLimit);
    standardError = sqrt(empiricalAcceptance*(1-empiricalAcceptance)/nDraws);

    naiveObjective = @(logD) stepNll(logD,acceptedStep, ...
        noiseCovariance,dt,driftTerm,[]);
    correctedObjective = @(logD) stepNll(logD,acceptedStep, ...
        noiseCovariance,dt,driftTerm,radialLimit);
    fitOptions = optimset('Display','off','TolX',1e-9);
    naiveD = exp(fminbnd(naiveObjective,log(1e-6),log(10),fitOptions));
    correctedD = exp(fminbnd(correctedObjective,log(1e-6),log(10),fitOptions));

    result = struct();
    result.trueDiffusionUm2PerS = trueD;
    result.empiricalAcceptance = empiricalAcceptance;
    result.quadratureAcceptance = quadratureAcceptance;
    result.acceptanceStandardError = standardError;
    result.naiveDiffusionUm2PerS = naiveD;
    result.correctedDiffusionUm2PerS = correctedD;
    result.naiveBiasPercent = 100*(naiveD/trueD-1);
    result.correctedBiasPercent = 100*(correctedD/trueD-1);
    result.passed = abs(empiricalAcceptance-quadratureAcceptance) < ...
        4*standardError && abs(result.correctedBiasPercent) < 5 && ...
        abs(result.correctedBiasPercent) < abs(result.naiveBiasPercent);
end

% -------------------------------------------------------------------------
function value = stepNll(logD,step,noiseCovariance,dt,driftTerm,radialLimit)
    D = exp(logD);
    covariance = noiseCovariance+2*D*dt*eye(2);
    residual = step-D*driftTerm;
    value = 0.5*size(step,1)*log(det(covariance)) + ...
        0.5*sum((residual/covariance).*residual,'all');
    if ~isempty(radialLimit)
        acceptance = gateAcceptance(D,noiseCovariance,dt, ...
            driftTerm,radialLimit);
        value = value+size(step,1)*log(max(acceptance,realmin));
    end
end

% -------------------------------------------------------------------------
function probability = gateAcceptance( ...
        D,noiseCovariance,dt,driftTerm,radialLimit)
    covariance = noiseCovariance+2*D*dt*eye(2);
    meanStep = D*driftTerm(:);
    nAngles = numel(radialLimit);
    angles = ((0:nAngles-1)+0.5)*(2*pi/nAngles);
    direction = [cos(angles);sin(angles)];
    solvedDirection = covariance\direction;
    solvedMean = covariance\meanStep;
    a = sum(direction.*solvedDirection,1);
    b = (direction.'*solvedMean).';
    c = meanStep.'*solvedMean;
    rootTwoA = sqrt(2*a);
    erfDifference = erf((a.*radialLimit-b)./rootTwoA) + ...
        erf(b./rootTwoA);
    gaussianIntegral = exp(-0.5*max(c-b.^2./a,0)) .* ...
        sqrt(pi./(2*a)).*erfDifference;
    endpointExponent = max(a.*radialLimit.^2-2*b.*radialLimit+c,0);
    radialIntegral = (b./a).*gaussianIntegral + ...
        (exp(-0.5*c)-exp(-0.5*endpointExponent))./a;
    probability = mean(max(radialIntegral,0))/sqrt(det(covariance));
    probability = min(max(probability,realmin),1);
end

% -------------------------------------------------------------------------
function printResults(results)
    fprintf('[fluctuating MIET verification]\n');
    fprintf('  drift maximum relative error: %.3g  pass=%d\n', ...
        results.drift.maximumRelativeError,results.drift.passed);
    disp(results.curvature.table);
    fprintf(['  gate acceptance: quadrature %.6f, Monte Carlo %.6f ' ...
        '+/- %.2g\n'],results.gate.quadratureAcceptance, ...
        results.gate.empiricalAcceptance,results.gate.acceptanceStandardError);
    fprintf(['  censored D: true %.4f, naive %.4f (%+.1f%%), ' ...
        'conditioned %.4f (%+.1f%%)  pass=%d\n'], ...
        results.gate.trueDiffusionUm2PerS, ...
        results.gate.naiveDiffusionUm2PerS,results.gate.naiveBiasPercent, ...
        results.gate.correctedDiffusionUm2PerS, ...
        results.gate.correctedBiasPercent,results.gate.passed);
    if results.endToEnd.wasRun
        fprintf(['  end-to-end: topography RMSE %.3g um, position RMSE ' ...
            '%.3g um, kappa %.3g 1/um, D %.3g um^2/s\n'], ...
            results.endToEnd.topographyRmseUm, ...
            results.endToEnd.positionRmseUm, ...
            results.endToEnd.curvaturePerUm, ...
            results.endToEnd.diffusionUm2PerS);
    end
    fprintf('  overall pass=%d\n',results.passed);
end
