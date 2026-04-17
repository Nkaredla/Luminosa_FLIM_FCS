function [y] = inverseTransformSampling(p,Nsamples,num)
%y = inverseTransformSampling(p, N)
% inverseTransformSampling performs resampling of a probability distribution 
% by inverse transform sampling.
% The function arguments are:
% p         = 	probability function p(x) to be sampled (will be normalized)
%               p must be a 1D vector
% Nsamples  = 	number of samples in the histogram y (default: 1000)
%
% The return parameters are:
% y         =	resampled histogram
% 
% The program needs no other m-files.
% (c) Sebastian Isbaner and Narain Karedla 2016

if verLessThan('matlab', '7.11')  % for version numbering, see http://en.wikipedia.org/wiki/MATLAB#Release_history
% -- Put code to run under MATLAB 7.10 and earlier here --
else
% -- Put code to run under MATLAB 7.11 and later here --
end 
if nargin < 3
    num = 1;
end
if nargin < 2
    Nsamples=1e3;
end
y = zeros(numel(p),num);
for i = 1:num
    p=p(:)/sum(p);
    p(p==0)=eps; % cdf must be strictly monotonically increasing!
    cdf=cumsum(p);
    u=rand(1,Nsamples);
    x=1:numel(p); x=x(:);
    samples = interp1([0;cdf],[x(1);x],u,'linear');
    
    y(:,i) = histc(samples,x-1);
end

end