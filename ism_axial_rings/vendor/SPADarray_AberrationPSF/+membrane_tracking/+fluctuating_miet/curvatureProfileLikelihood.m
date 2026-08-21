function identifiability = curvatureProfileLikelihood(gpData, gpFit, opts)
    import membrane_tracking.fluctuating_miet.*

%CURVATUREPROFILELIKELIHOOD Amplitude-profile interval for kappa.
%
%   The Fisher error bar assumes the log-likelihood is quadratic in kappa
%   at fixed field parameters. When a field mode is confounded with
%   curvature that assumption fails, so the interval is taken from the
%   amplitude-profile likelihood, re-fitting the field rms at every trial
%   kappa while holding correlation scales fixed. The reported ratio
%   profileSigma/fisherSigma above about 1.5 is a sign of confounding.
    centre = gpFit.curvaturePerUm;
    halfWidth = opts.curvatureProfileHalfWidth * ...
        max(gpFit.curvatureSigmaPerUm, 0.05*abs(centre) + 1e-3);
    grid = linspace(centre - halfWidth, centre + halfWidth, ...
        opts.curvatureProfilePoints).';
    profile = nan(size(grid));
    for k = 1:numel(grid)
        objective = @(logAmplitude) firstOutput(@() ...
            gpProfileNegLogLikelihood( ...
                [logAmplitude, gpFit.logParameters(2:3)], gpData, opts, ...
                grid(k), grid(k)));
        fitOptions = optimset('Display','off','MaxIter',60, ...
            'MaxFunEvals',120,'TolX',1e-5,'TolFun',1e-5);
        [~, value, exitFlag] = fminsearch(objective, ...
            gpFit.logParameters(1), fitOptions);
        if exitFlag > 0 && isfinite(value)
            profile(k) = value;
        end
    end
    identifiability = struct();
    identifiability.curvatureGridPerUm = grid;
    identifiability.profileNegLogLikelihood = profile;
    identifiability.profileSigmaPerUm = NaN;
    identifiability.profileInterval95PerUm = [NaN NaN];
    identifiability.intervalHitsGridBoundary = false;
    finite = isfinite(profile);
    if sum(finite) >= 5
        deltaNll = profile - min(profile(finite));
        inside = finite & deltaNll <= 1.92;   % 95% for one parameter
        if any(inside)
            identifiability.profileInterval95PerUm = ...
                [min(grid(inside)), max(grid(inside))];
            identifiability.profileSigmaPerUm = ...
                diff(identifiability.profileInterval95PerUm) / 3.92;
            identifiability.intervalHitsGridBoundary = ...
                inside(1) || inside(end);
        end
    end
    identifiability.fisherSigmaPerUm = gpFit.curvatureSigmaPerUm;
    identifiability.profileOverFisherSigmaRatio = ...
        identifiability.profileSigmaPerUm / ...
        max(gpFit.curvatureSigmaPerUm, realmin);
    identifiability.fieldCorrelationLengthUm = gpFit.correlationLengthUm;
    identifiability.observationRadiusUm = ...
        max(hypot(gpData.position(:,1), gpData.position(:,2)));
    identifiability.confoundingRatio = gpFit.correlationLengthUm / ...
        max(identifiability.observationRadiusUm, realmin);
    identifiability.interpretation = ['confoundingRatio near 1 means the ' ...
        'field varies on the scale of the observed patch and is quadratic ' ...
        'across it, so kappa and the field are not separable. Much less ' ...
        'than 1: fluctuations average out. Much greater than 1: the field ' ...
        'acts as a piston and is absorbed into the tip height. The interval ' ...
        'profiles field amplitude but is conditional on correlation scales.'];
end
