function [hEff, hExc, eta] = spadEffectivePSFArray(sim, coeffs, normalizeMode)
%--------------------------------------------------------------------------
% spadEffectivePSFArray
%
% PURPOSE
%   Compute the 3D effective excitation-detection PSF for every detector
%   pixel in the SPAD array.
%
% INPUTS
%   sim           : simulation structure from defaultParams()
%   coeffs        : aberration coefficient struct or numeric vector
%   normalizeMode : 'none', 'each', or 'global' [default: 'none']
%
% OUTPUTS
%   hEff(y,x,z,k) : effective PSF for detector channel k
%   hExc(y,x,z)   : excitation PSF
%   eta(y,x,z,k)  : detector collection efficiency
%--------------------------------------------------------------------------

    if nargin < 3 || isempty(normalizeMode)
        normalizeMode = 'none';
    end
    normalizeMode = char(normalizeMode);

    coeffs = coeffStruct(sim, coeffs);

    hExc = psfBessel(sim, coeffs, sim.lamExc);
    hDet = psfBessel(sim, coeffs, sim.lamEm);
    eta = detectorCollectionEfficiencyExplicit(sim, hDet);

    nCh = size(sim.detXY, 1);
    hEff = zeros(size(hExc,1), size(hExc,2), size(hExc,3), nCh);

    for k = 1:nCh
        hEff(:,:,:,k) = hExc .* eta(:,:,:,k);
    end

    switch lower(normalizeMode)
        case 'none'
            return;

        case 'global'
            m = max(hEff(:));
            if m > 0
                hEff = hEff / m;
            end

        case 'each'
            for k = 1:nCh
                hk = hEff(:,:,:,k);
                m = max(hk(:));
                if m > 0
                    hEff(:,:,:,k) = hEff(:,:,:,k) / m;
                end
            end

        otherwise
            error('spadEffectivePSFArray:BadNormalizeMode', ...
                'normalizeMode must be ''none'', ''each'', or ''global''.');
    end
end
