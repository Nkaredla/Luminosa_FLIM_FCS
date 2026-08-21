function pred = predictCollar(model, collar)
%--------------------------------------------------------------------------
% predictCollar
%
% PURPOSE
%   Evaluate the fitted quadratic collar model at arbitrary collar settings.
%
% INPUTS
%   model  : structure from fitQuadraticCollarModel
%   collar : vector of new collar settings
%
% OUTPUT
%   pred : predicted coefficient matrix
%--------------------------------------------------------------------------

    % Shift collar values relative to reference c0
    u = collar(:) - model.c0;

    % Evaluate quadratic model
    pred = model.a0 + u.*model.g1 + (u.^2).*model.g2;
end