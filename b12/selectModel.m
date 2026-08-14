function [ok, exponentialCount, tauNs, tauMeanNs, coefficients, ...
        speciesFraction, fitCounts, metrics, message] = selectModel( ...
        counts, irf, dtNs, maximumExponentials, oneStateSeed, twoStateSeed, ...
        lifetimeBounds, includeBackground, minimumPhotons, maximumIterations, ...
        requiredBicImprovement, minimumSeparationFraction, minimumAmplitudeFraction)
%SELECTMODEL Fit one/two exponentials and select the parsimonious BIC model.
% metrics = [selected deviance, selected BIC, selected reduced chi-square,
%            photons, one-exp BIC, two-exp BIC, one-exp chi2, two-exp chi2].

    counts = max(double(counts(:)), 0);
    photonCount = sum(counts);
    ok = false;
    exponentialCount = 0;
    tauNs = nan(1, 2);
    tauMeanNs = NaN;
    coefficients = zeros(1, 3);
    speciesFraction = nan(1, 2);
    fitCounts = zeros(size(counts));
    metrics = [NaN, inf, NaN, photonCount, inf, inf, NaN, NaN];
    if photonCount < minimumPhotons
        message = sprintf('Too few photons (%d < %d).', round(photonCount), minimumPhotons);
        return;
    end

    [tau1, mean1, coeff1, fraction1, fitted1, deviance1, bic1, chi1] = ...
        fitFixedOrder(counts, irf, dtNs, 1, oneStateSeed, ...
        twoStateSeed, lifetimeBounds, includeBackground, maximumIterations);
    selectedOrder = 1;
    selectedValues = {tau1, mean1, coeff1, fraction1, fitted1, deviance1, bic1, chi1};
    bic2 = inf;
    chi2 = NaN;

    if maximumExponentials >= 2
        [tau2, mean2, coeff2, fraction2, fitted2, deviance2, bic2, chi2] = ...
            fitFixedOrder(counts, irf, dtNs, 2, oneStateSeed, ...
            twoStateSeed, lifetimeBounds, includeBackground, maximumIterations);
        separation = abs(tau2(2) - tau2(1)) / max(mean(tau2(1:2)), eps);
        validTwoState = separation >= minimumSeparationFraction && ...
            min(fraction2(1:2)) >= minimumAmplitudeFraction;
        if validTwoState && bic2 < bic1 - requiredBicImprovement
            selectedOrder = 2;
            selectedValues = {tau2, mean2, coeff2, fraction2, fitted2, deviance2, bic2, chi2};
        end
    end

    ok = true;
    exponentialCount = selectedOrder;
    tauNs = selectedValues{1};
    tauMeanNs = selectedValues{2};
    coefficients = selectedValues{3};
    speciesFraction = selectedValues{4};
    fitCounts = selectedValues{5};
    metrics = [selectedValues{6}, selectedValues{7}, selectedValues{8}, ...
        photonCount, bic1, bic2, chi1, chi2];
    message = '';
end
