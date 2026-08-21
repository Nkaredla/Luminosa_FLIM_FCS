function model = fitQuadraticCollarModel(collar, coeffs, ridge)
%--------------------------------------------------------------------------
% fitQuadraticCollarModel
%
% PURPOSE
%   Fit a quadratic model for Zernike coefficients as a function of collar
%   setting:
%
%       a(c) = a0 + g1*(c-c0) + g2*(c-c0)^2
%
% INPUTS
%   collar : vector of collar settings
%   coeffs : matrix [Nsettings x Nmodes]
%   ridge  : ridge regularization strength
%
% OUTPUT
%   model : structure with fields
%       c0, a0, g1, g2
%--------------------------------------------------------------------------

    if nargin < 3
        ridge = 1e-6;
    end

    % Use median collar as reference c0
    c0 = median(collar(:));

    % Center collar values
    u = collar(:) - c0;

    % Quadratic design matrix
    X = [ones(size(u)) u u.^2];

    % Ridge-regularized least-squares fit
    beta = (X.'*X + ridge*eye(3)) \ (X.'*coeffs);

    % Store fitted parameters
    model.c0 = c0;
    model.a0 = beta(1,:);
    model.g1 = beta(2,:);
    model.g2 = beta(3,:);
end