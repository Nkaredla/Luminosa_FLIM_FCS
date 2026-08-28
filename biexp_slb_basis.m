function basis = biexp_slb_basis(irf, dtNs, periodNs, nBin, precision)
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
% PRECISION. The optional fifth argument ('double' by default, 'single'
% accepted) sets the class of C and of the time axis, and everything downstream
% inherits it because MATLAB promotes to the wider operand. Single is ample for
% 156-bin photon counts and halves the memory traffic, which is where this fit
% is actually bound: the same code measured 4.8 ms/pixel on 19000-pixel blocks
% against 2.0 ms/pixel on 4000-pixel blocks, which is a bandwidth signature and
% not an arithmetic one.
%
% It also decides whether a GPU is worth using at all here. This machine's
% Quadro T1000 Max-Q runs FP64 at 1/32 of its FP32 rate, so in double it is
% slower than the CPU - which is the likely explanation for an earlier GPU
% attempt in this project measuring 0.41-0.64x. Single is the precondition, not
% a refinement.
%
% FIELDS
%   C            nBin-by-nBin circulant of the IRF
%   timeNs       nBin-by-1 bin times
%   nBin, dtNs, periodNs
%   periodFactor function handle tau -> sum_{r=0..3} exp(-r*periodNs/tau)

    irf = double(irf(:));
    if nargin < 5 || isempty(precision); precision = 'double'; end
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
    basis.C = cast(C, precision);
    basis.precision = precision;
    basis.timeNs = cast((0:nBin - 1)' * dtNs, precision);
    basis.nBin = nBin;
    basis.dtNs = dtNs;
    basis.periodNs = periodNs;
    basis.periodFactor = @(tau) sum(exp(-(0:3)' * (periodNs ./ tau(:)')), 1);
end
