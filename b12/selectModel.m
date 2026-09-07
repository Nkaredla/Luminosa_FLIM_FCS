function [ok, exponentialCount, tauNs, tauMeanNs, coefficients, ...
        speciesFraction, fitCounts, metrics, message] = selectModel( ...
        counts, irf, dtNs, maximumExponentials, oneStateSeed, twoStateSeed, ...
        lifetimeBounds, includeBackground, minimumPhotons, maximumIterations, ...
        requiredBicImprovement, minimumSeparationFraction, minimumAmplitudeFraction, ...
        minimumTwoExpPhotons, tailOffsetNs)
%SELECTMODEL Fit one/two exponentials and select the parsimonious BIC model.
% Optional minimumTwoExpPhotons defaults to 500. Optional tailOffsetNs fits
% only bins at least that far after this curve's peak, without reconvolution.
% Eligibility and metrics count only photons in the fitted window. Excluded
% bins and skipped fits have NaN fitCounts. Hard floors are 10/500 photons.
% metrics = [selected deviance, selected BIC, selected reduced chi-square,
%            photons, one-exp BIC, two-exp BIC, one-exp chi2, two-exp chi2].

    if nargin < 14 || isempty(minimumTwoExpPhotons)
        minimumTwoExpPhotons = 500;
    end
    if nargin < 15
        tailOffsetNs = [];
    end
    minimumPhotons = max(10, minimumPhotons);
    minimumTwoExpPhotons = max(500, minimumTwoExpPhotons);
    counts = max(double(counts(:)), 0);
    originalBinCount = numel(counts);
    firstFitBin = 1;
    if ~isempty(tailOffsetNs)
        validateattributes(tailOffsetNs, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'nonnegative'});
        validateattributes(dtNs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
        [~, peakBin] = max(counts);
        if isempty(peakBin)
            peakBin = 1;
        end
        % Round up: the first retained bin is at least the offset after peak.
        firstFitBin = min(originalBinCount + 1, peakBin + ceil(tailOffsetNs / dtNs));
        counts = counts(firstFitBin:end);
        % A delta IRF produces unconvolved exponentials on the retained tail.
        irf = zeros(size(counts));
        if ~isempty(irf)
            irf(1) = 1;
        end
    end
    photonCount = sum(counts);
    ok = false;
    exponentialCount = 0;
    tauNs = nan(1, 2);
    tauMeanNs = NaN;
    coefficients = zeros(1, 3);
    speciesFraction = nan(1, 2);
    fitCounts = nan(originalBinCount, 1);
    metrics = [NaN, inf, NaN, photonCount, inf, inf, NaN, NaN];
    if photonCount < minimumPhotons
        message = sprintf('Too few photons (%d < %d).', round(photonCount), minimumPhotons);
        return;
    end

    if numel(counts) <= 2 + double(includeBackground)
        message = 'Too few retained TCSPC bins for a single-exponential fit.';
        return;
    end

    [tau1, mean1, coeff1, fraction1, fitted1, deviance1, bic1, chi1] = ...
        fitFixedOrder(counts, irf, dtNs, 1, oneStateSeed, ...
        twoStateSeed, lifetimeBounds, includeBackground, maximumIterations);
    selectedOrder = 1;
    selectedValues = {tau1, mean1, coeff1, fraction1, fitted1, deviance1, bic1, chi1};
    bic2 = inf;
    chi2 = NaN;

    if maximumExponentials >= 2 && photonCount >= minimumTwoExpPhotons && ...
            numel(counts) > 4 + double(includeBackground)
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
    fitCounts(firstFitBin:end) = selectedValues{5};
    metrics = [selectedValues{6}, selectedValues{7}, selectedValues{8}, ...
        photonCount, bic1, bic2, chi1, chi2];
    message = '';
end
