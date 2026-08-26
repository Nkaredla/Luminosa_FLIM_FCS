function basis = biexp_slb_basis(irf, dtNs, periodNs, nBin)
%BIEXP_SLB_BASIS Precompute what every pattern build shares.
%
% basis = biexp_slb_basis(irf, dtNs, periodNs, nBin)
%
% biexp_slb_pattern.m does conv(irf, decay) and then folds the result back into
% nBin bins. That pair of operations IS circular convolution: conv returns
% 2*nBin-1 samples and the fold adds samples nBin+1..2*nBin-1 onto 1..nBin-1,
% which is exactly the wrap of a length-nBin cyclic convolution. So the whole
% build is one fixed matrix acting on the decay vector, and the matrix does not
% depend on the lifetime.
%
% That matters for speed in a specific way. Every pixel in a variable-projection
% fit has its OWN lifetime, so patterns cannot be precomputed per grid node the
% way the grid pipeline does. But because C is shared, P pixels' patterns are
% C * G with G the 156-by-P matrix of their decays - one matrix product for the
% whole block. This is what lets a per-pixel optimizer vectorise as well as the
% grid does.
%
% A second simplification, verified rather than assumed: the four-period sum in
% biexp_slb_pattern.m has NO effect on the returned pattern. Since
% exp(-(t + r*T)/tau) = exp(-t/tau) * exp(-r*T/tau), the four repeats contribute
% only a scalar factor, and the unit-sum normalisation divides it straight back
% out (checked over 65 lifetimes spanning 0.02-8 ns: max relative difference
% 5.9e-16). It is NOT dead code though - it changes the pre-normalisation sum,
% which is what converts an amplitude into a pre-exponential species weight - so
% the factor is carried separately in periodFactor rather than dropped.
%
% FIELDS
%   C            nBin-by-nBin circulant of the IRF
%   timeNs       nBin-by-1 bin times
%   nBin, dtNs, periodNs
%   periodFactor function handle tau -> sum_{r=0..3} exp(-r*periodNs/tau)

    irf = double(irf(:));
    irf = max(irf, 0);
    if sum(irf) > 0; irf = irf / sum(irf); end
    if numel(irf) ~= nBin
        error('biexp_slb_basis:IrfLength', ...
            'The IRF has %d bins but nBin is %d.', numel(irf), nBin);
    end

    % Circulant: column j is the IRF shifted down by j-1, wrapping.
    C = zeros(nBin, nBin);
    for j = 1:nBin
        C(:, j) = circshift(irf, j - 1);
    end

    basis = struct();
    basis.C = C;
    basis.timeNs = (0:nBin - 1)' * dtNs;
    basis.nBin = nBin;
    basis.dtNs = dtNs;
    basis.periodNs = periodNs;
    basis.periodFactor = @(tau) sum(exp(-(0:3)' * (periodNs ./ tau(:)')), 1);
end
