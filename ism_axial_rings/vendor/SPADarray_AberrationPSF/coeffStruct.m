function coeffs = coeffStruct(sim, a)
%--------------------------------------------------------------------------
% coeffStruct
%
% PURPOSE
%   Convert either:
%       1) a numeric coefficient vector
%       2) an already-defined coefficient structure
%   into a coefficient structure with named fields.
%
% INPUT
%   sim : simulation structure containing the mode order
%   a   : either a structure or a numeric vector
%
% OUTPUT
%   coeffs : structure, e.g. coeffs.defocus = 0.18
%
% NOTES
%   If a is already a structure, it is returned unchanged.
%   If a is numeric, only nonzero entries are stored as fields.
%--------------------------------------------------------------------------

    coeffs = struct();

    % If already a struct, simply return it
    if isstruct(a)
        coeffs = a;
        return;
    end

    % Convert numeric vector into named fields
    n = min(numel(a), numel(sim.modeOrder));
    for k = 1:n
        if abs(a(k)) > 1e-15
            coeffs.(sim.modeOrder{k}) = a(k);
        end
    end
end
