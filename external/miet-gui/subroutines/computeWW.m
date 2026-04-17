function [ ww ] = computeWW(mm,apd,epsilon)
% Usage: [ ww ] = computeWW(mm,apd,epsilon) [MEX Function]
%
% Equivalent MATLAB code:
% ww = [repmat(exp(-mm(:,apd)),1,apd-1) + exp(-mm(:,1:apd-1))/(exp(epsilon)-1), exp(-mm(:,apd:end))/(1-exp(-epsilon))];
%
% Author: Simon Christoph Stein
% Date: Jan. 2016
% E-Mail: scstein@phys.uni-goettingen.de
end

