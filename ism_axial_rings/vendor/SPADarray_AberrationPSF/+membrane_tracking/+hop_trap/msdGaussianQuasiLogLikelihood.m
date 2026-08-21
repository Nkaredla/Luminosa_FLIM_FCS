function value = msdGaussianQuasiLogLikelihood(observed, predicted, sigma)
%MSDGAUSSIANQUASILOGLIKELIHOOD Weighted Gaussian MSD quasi-likelihood.

    residual = (observed - predicted) ./ sigma;
    value = -0.5 * sum(residual.^2 + log(2*pi*sigma.^2));
end
