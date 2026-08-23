function opts = ring_fit_defaults(opts, nComponent)
%RING_FIT_DEFAULTS Shared option defaults for the ring TCSPC fitters.
%
% opts = ring_fit_defaults(opts, nComponent)
%
% Kept in one place so the free-amplitude and height-constrained fitters cannot
% drift apart in their optimiser settings. That matters more than usual here,
% because the whole study compares those two fitters against each other: any
% asymmetry in the optimiser would be indistinguishable from a difference in
% information content.

    defaults = struct( ...
        'restarts', 1, ...
        'seedTauNs', [], ...
        'seedHeightUm', [], ...
        'heightGroups', [], ...
        'evalsPerParam', 400, ...
        'evalBudget', [], ...
        'tolFun', 1e-4, ...
        'tolX', 1e-4, ...
        'tauBoundsNs', [0.05 20]);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end
    if isempty(opts.seedTauNs)
        opts.seedTauNs = logspace(log10(0.3), log10(4), nComponent);
    end
    if numel(opts.seedTauNs) ~= nComponent
        error('ring_fit_defaults:SeedSize', ...
            'seedTauNs must have %d elements.', nComponent);
    end
end
