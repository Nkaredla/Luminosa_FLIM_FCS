function result = emptyModelOrderResult()
    import membrane_tracking.focused_ism.*

    result = struct('nEmitters', 0, 'posteriorProbability', NaN, ...
        'entropy', NaN, 'candidateEmitterCounts', zeros(0,1), ...
        'posteriorWeights', zeros(0,1));
end
