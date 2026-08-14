function value = fitObjective(logLifetime, counts, irf, dtNs, ...
        includeBackground, lowerBound, upperBound)
%FITOBJECTIVE Profile Poisson objective for nonlinear lifetime optimisation.
% Linear nonnegative amplitudes are solved exactly for every lifetime trial.

    lifetime = exp(double(logLifetime(:)));
    lifetime = min(max(lifetime, lowerBound), upperBound);
    model = decayModel(irf, dtNs, lifetime, includeBackground);
    coefficients = lsqnonneg(model, counts);
    value = poissonDeviance(counts, max(model * coefficients, eps));
end
