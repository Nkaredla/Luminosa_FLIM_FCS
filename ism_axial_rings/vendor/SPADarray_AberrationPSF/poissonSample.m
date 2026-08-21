
function data = poissonSample(mu)
%--------------------------------------------------------------------------
% poissonSample
%
% PURPOSE
%   Generate Poisson-distributed photon count data.
%
% INPUT
%   mu : mean count image / stack
%
% OUTPUT
%   data : noisy sampled image / stack
%
% NOTES
%   If poissrnd is unavailable, a Gaussian approximation is used.
%--------------------------------------------------------------------------

    if exist('poissrnd', 'file') == 2
        data = poissrnd(mu);
    else
        data = max(mu + sqrt(max(mu,1)).*randn(size(mu)), 0);
    end
end