function spacing = estimateStepObservationSpacing(stepTable)
%ESTIMATESTEPOBSERVATIONSPACING Estimate time between independent HMM steps.

    values = zeros(0, 1);
    sequences = unique(stepTable.sequenceId, 'stable').';
    for sequence = sequences
        time = stepTable.timeCenterS(stepTable.sequenceId == sequence);
        if numel(time) > 1
            values = [values; diff(time)]; %#ok<AGROW>
        end
    end
    if isempty(values)
        spacing = 2 * median(stepTable.stepDtS);
    else
        spacing = median(values);
    end
end
