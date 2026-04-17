function [beta,k] = PIRLSnonneg(x,y,max_num_iter, TINY)
% Poisson Iterative Reweighted Least Squares
% solves X*beta=Y for Poisson distributed data Y for nonnegative beta
% x             predictors
% y             data
% max_num_iter  maximum number of iterations
% TINY          regularization parameter
% 
% beta          solution vector, beta = X^{-1} Y
% k             actual number of iterations neede
% 
% Sebastian Isbaner, 04/2018

if nargin<3 || isempty(max_num_iter)
    max_num_iter=10;
end
[n,~] = size(y);
if nargin<4 || isempty(TINY)
    TINY = 0.1/n;
end

%% Original
% w = zeros(n);
% beta_last = lsqnonneg(x,y); % initial guess from unweighted regression
% for k=1:max_num_iter
%     w(1:n+1:end) = 1./max((x*beta_last),TINY); % updated weight matrix (with regularization)
%     beta = lsqnonneg(x'*w*x,x'*w*y); % maximum likelihood solution vector
%     delta = beta-beta_last;
%     %check convergence
%     if (sum(delta.^2) < 1e-10), break; end
%     beta_last = beta;
% end

%% Use implicit expansion, inline lsqnonneg and remove option/error handeling
negFlag = false;
[beta_last,negFlag] = lsqnonneg_short(x,y,negFlag); % initial guess from unweighted regression
for k=1:max_num_iter
%     wlin = 1./max((x*beta_last),TINY); % updated weight matrix (with regularization)
%     [beta,negFlag] = lsqnonneg_short((x'.*wlin')*x,(x'.*wlin')*y,negFlag); % maximum likelihood solution vector
    xwlin = (x./max((x*beta_last),TINY))'; % predictors times updated weight matrix (with regularization)
    [beta,negFlag] = lsqnonneg_short(xwlin*x,xwlin*y,negFlag); % maximum likelihood solution vector
    delta = beta-beta_last;
    %check convergence
    if (sum(delta.^2) < 1e-10), break; end
    beta_last = beta;
end
%% Enforce nonneg by resetting to zero. Faster but result deviates by approx. 1% for noisy data.
% beta = x\y; % initial guess from unweighted regression
% % beta(:,1) = lsqnonneg(x,y); % initial guess from unweighted regression
% for k=1:max_num_iter
%     yy = y - x*beta; % residual matrix
%     wlin = 1./max((x*beta),TINY); % updated weight matrix (with regularization)
%     delta = ((x'.*wlin')*x)\((x'.*wlin')*yy); % maximum likelihood solution vector
%     beta = beta + delta; % update betas
%     %check convergence
%     if (sum(delta.^2) < 1e-10) && all(beta>=0), break; end
%     beta(beta<0) = 0;
% end

function [x,negFlag] = lsqnonneg_short(C,d,negFlag)
%LSQNONNEG Linear least squares with nonnegativity constraints.
%   X = LSQNONNEG(C,d) returns the vector X that minimizes NORM(d-C*X)
%   subject to X >= 0. C and d must be real.
%
%   X = LSQNONNEG(C,d,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details. Used
%   options are Display and TolX. (A default tolerance TolX of 
%   10*MAX(SIZE(C))*NORM(C,1)*EPS is used.) 
%   
%   X = LSQNONNEG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the matrix 'C' in PROBLEM.C, the vector 'd' in
%   PROBLEM.d, the options structure in PROBLEM.options, and solver name
%   'lsqnonneg' in PROBLEM.solver. 
%
%   [X,RESNORM] = LSQNONNEG(...) also returns the value of the squared 2-norm of 
%   the residual: norm(d-C*X)^2.
%
%   [X,RESNORM,RESIDUAL] = LSQNONNEG(...) also returns the value of the  
%   residual: d-C*X.
%   
%   [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONNEG(...) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%    1  LSQNONNEG converged with a solution X.
%    0  Iteration count was exceeded. Increasing the tolerance
%       (OPTIONS.TolX) may lead to a solution.
%  
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONNEG(...) returns a structure
%   OUTPUT with the number of steps taken in OUTPUT.iterations, the type of 
%   algorithm used in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONNEG(...) returns 
%   the dual vector LAMBDA  where LAMBDA(i) <= 0 when X(i) is (approximately) 0 
%   and LAMBDA(i) is (approximately) 0 when X(i) > 0.
% 
%   See also LSCOV, SLASH.

%   Copyright 1984-2016 The MathWorks, Inc. 

% Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

% Modified by CT, 2019. Optimised for speed: Removed unnecessary checks and
% options; added precheck for negative components.

% Test whether we really need constrains. Incase we actually do need them 
% this makes the fit ~50% slower.
if nargin>2 && ~negFlag
    x = C\d;
    if all(x>=0)
        return;
    end
end
negFlag = false;

% Default tolerance
tol = 10*eps*norm(C,1)*length(C);

n = size(C,2);
% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
wz = nZeros;

% Initialize set of non-active columns to null
P = false(n,1);
% Initialize set of active columns to all and the initial point to zeros
Z = true(n,1);
x = nZeros;

resid = d - C*x;
w = C'*resid;

% Set up iteration criterion
outeriter = 0;
iter = 0;
itmax = 3*n;

% Outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > tol)
   outeriter = outeriter + 1;
   % Reset intermediate solution z
   z = nZeros; 
   % Create wz, a Lagrange multiplier vector of variables in the zero set.
   % wz must have the same size as w to preserve the correct indices, so
   % set multipliers to -Inf for variables outside of the zero set.
   wz(P) = -Inf;
   wz(Z) = w(Z);
   % Find variable with largest Lagrange multiplier
   [~,t] = max(wz);
   % Move variable t from zero set to positive set
   P(t) = true;
   Z(t) = false;
   % Compute intermediate solution using only variables in positive set
   z(P) = C(:,P)\d;
   % inner loop to remove elements from the positive set which no longer belong
   while any(z(P) <= 0)
       negFlag = true;
       iter = iter + 1;
       if iter > itmax
           msg = getString(message('MATLAB:optimfun:lsqnonneg:IterationCountExceeded'));
           disp(msg)
           x = z;
           return
       end
       % Find indices where intermediate solution z is approximately negative
       Q = (z <= 0) & P;
       % Choose new x subject to keeping new x nonnegative
       alpha = min(x(Q)./(x(Q) - z(Q)));
       x = x + alpha*(z - x);
       % Reset Z and P given intermediate values of x
       Z = ((abs(x) < tol) & P) | Z;
       P = ~Z;
       z = nZeros;           % Reset z
       z(P) = C(:,P)\d;      % Re-solve for z
   end
   x = z;
   resid = d - C*x;
   w = C'*resid;
end
