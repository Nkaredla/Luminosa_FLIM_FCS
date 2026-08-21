function out = simulateArray93FisherInformation(varargin)
%SIMULATEARRAY93FISHERINFORMATION Poisson Fisher/CRB on the 93-pixel array vs photons and aberration.
%
%   out = simulateArray93FisherInformation()
%   out = simulateArray93FisherInformation('maxZernikeOrder', 7)
%   out = simulateArray93FisherInformation('photonCounts', 10.^(2:6))
%
%   Thin wrapper around simulateArrayFisherInformation for the 93-channel
%   honeycomb array (emitter on the coverslip glass, high-NA oil objective by
%   default). See simulateArrayFisherInformation for the full option list and
%   the Poisson Fisher-information / Cramer-Rao conventions (identical to
%   simulateRetrieveCenterMicroimageAberrations.m).
%
%   The default maxZernikeOrder here is 4 (through primary spherical) to keep
%   the finite-difference Jacobian light; pass 'maxZernikeOrder', 7 for the
%   full radial-order-7 set used in compareArrayAberrationFisher. Any option
%   passed here overrides these defaults.

    out = simulateArrayFisherInformation( ...
        'detectorLayout', 'honeycomb93', ...
        'maxZernikeOrder', 4, ...
        varargin{:});
end
