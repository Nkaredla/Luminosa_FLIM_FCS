
function [dataN, beta] = normalizeDataStack(data)
%--------------------------------------------------------------------------
% normalizeDataStack
%
% PURPOSE
%   Estimate a simple background level from the image borders and normalize
%   the data stack after background subtraction.
%
% INPUT
%   data : raw detector stack
%
% OUTPUTS
%   dataN : background-subtracted and normalized stack
%   beta  : estimated background level
%
% DESCRIPTION
%   Background is estimated as the median value over border pixels from all
%   detector channels. After subtraction, the data are clipped at zero and
%   normalized to sum to 1.
%--------------------------------------------------------------------------

    % Collect border pixels from all detector channels
    border = [reshape(data(1,:,:), [], 1); ...
              reshape(data(end,:,:), [], 1); ...
              reshape(data(:,1,:), [], 1); ...
              reshape(data(:,end,:), [], 1)];

    % Median border value as background estimate
    beta = max(median(border), 0);

    % Subtract background and clip negative values
    d = max(data - beta, 0);

    % Normalize to unit total signal
    s = sum(d(:));
    if s <= 0
        dataN = d;
    else
        dataN = d / s;
    end
end