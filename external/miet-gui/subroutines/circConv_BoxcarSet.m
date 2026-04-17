function [ convolutionMatrix ] = circConv_BoxcarSet( signal )
% Usage: [ convolutionMatrix ] = circConv_BoxcarSet( signal )  [MEX Function]
%
% Performs a circular convolution of the signal with boxcar filters 
% (e.g. [1,1,1,1,0,0,0,...]) of length 1 to numel(signal).
%
% -- Equivalent MATLAB code --
% boxcarSet_fft = fft(gallery('triw',numel(signal),1,numel(signal)));
% convolutionMatrix = real(ifft(repmat(fft(signal),[1,numel(signal)]).*boxcarSet_fft));
%
% Input:
%    signal: N-element vector
% Output:
%    convolutionMatrix: NxN matrix where convolutionMatrix(:, boxcarLength)
%                       is the convolution with the boxcar filter made of boxcarLength '1's.
%                       Note: convolutionMatrix(:,1) is the input signal.
%
% Author: Simon Christoph Stein
% Date: Jan. 2016
% E-Mail: scstein@phys.uni-goettingen.de
end

