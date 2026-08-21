function out = simulateArray23FisherInformation(varargin)
%SIMULATEARRAY23FISHERINFORMATION Poisson Fisher/CRB on the 23-pixel array vs photons and aberration.
%
%   out = simulateArray23FisherInformation()
%   out = simulateArray23FisherInformation('maxZernikeOrder', 5)
%   out = simulateArray23FisherInformation('photonCounts', 10.^(2:6))
%
%   Thin wrapper around simulateArrayFisherInformation for the 23-channel
%   Luminosa honeycomb array. By default it uses the same emitter-on-glass
%   high-NA oil optics as the 93-array wrapper so the two are directly
%   comparable (compareArrayAberrationFisher); pass 'sampleGeometry',
%   'homogeneous' with 'NA',1.2,'nMedium',1.33 for the legacy water context.
%
%   The default maxZernikeOrder here is 5 (through pentafoil), matching the
%   23-array mode budget. Any option passed here overrides these defaults.

    out = simulateArrayFisherInformation( ...
        'detectorLayout', 'honeycomb23', ...
        'maxZernikeOrder', 5, ...
        varargin{:});
end
