function [beta,k] = PIRLSnonneg(x,y,max_num_iter)
% Poisson Iterative Reweighted Least Squares
% solves X*beta=Y for Poisson distributed data Y for nonnegative beta
% x             predictors
% y             data
% max_num_iter  maximum number of iterations
% 
% beta          solution vector, beta = X^{-1} Y
% k             actual number of iterations neede
% 
% Sebastian Isbaner, 04/2018

if nargin<3 || isempty(max_num_iter)
    max_num_iter=10;
end
[n,~] = size(y);
TINY = 0.1/n;
w = zeros(n);
beta_last = lsqnonneg(x,y); % initial guess from unweighted regression
for k=1:max_num_iter
    w(1:n+1:end) = 1./max((x*beta_last),TINY); % updated weight matrix (with regularization)
    beta = lsqnonneg(x'*w*x,x'*w*y); % maximum likelihood solution vector
    delta = beta-beta_last;
    %check convergence
    if (sum(delta.^2) < 1e-10), break; end
    beta_last = beta;
end
