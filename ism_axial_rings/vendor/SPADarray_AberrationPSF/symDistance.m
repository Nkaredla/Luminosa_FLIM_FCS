function d = symDistance(A, B)
%--------------------------------------------------------------------------
% symDistance
%
% PURPOSE
%   Compute the Euclidean distance between two normalized detector stacks.
%
% INPUTS
%   A, B : image stacks
%
% OUTPUT
%   d : L2 distance between normalized stacks
%
% DESCRIPTION
%   Both stacks are normalized to sum 1 before comparison.
%--------------------------------------------------------------------------

    A = A / sum(A(:));
    B = B / sum(B(:));
    d = norm(A(:) - B(:), 2);
end