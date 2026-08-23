function options = ring_search_options(opts, nParam)
%RING_SEARCH_OPTIONS fminsearch options for the ring TCSPC fitters.
%
% options = ring_search_options(opts, nParam)
%
% By default the evaluation budget scales with the number of searched
% parameters (evalsPerParam), which keeps the search effort PER PARAMETER equal
% across models. Set opts.evalBudget to give every model the same ABSOLUTE
% budget instead. The two choices bias the comparison in opposite directions,
% so the study reports which one was used:
%
%   per-parameter  the height-constrained model, having twice as many searched
%                  parameters, gets twice the evaluations - it could win partly
%                  on optimiser effort
%   absolute       the height-constrained model must cover a space of twice the
%                  dimension on the same budget - it is handicapped
%
% If a conclusion survives both, the conclusion is about information and not
% about fminsearch.

    if ~isempty(opts.evalBudget)
        maxEvals = opts.evalBudget;
        maxIter = opts.evalBudget;
    else
        maxEvals = opts.evalsPerParam * nParam;
        maxIter = round(0.5 * opts.evalsPerParam) * nParam;
    end
    options = optimset('Display', 'off', ...
        'MaxFunEvals', maxEvals, 'MaxIter', maxIter, ...
        'TolX', opts.tolX, 'TolFun', opts.tolFun);
end
