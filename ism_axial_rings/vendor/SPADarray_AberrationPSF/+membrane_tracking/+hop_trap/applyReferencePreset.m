function custom = applyReferencePreset(custom, usingDefaults)
%APPLYREFERENCEPRESET Fill model defaults from a named literature preset.

    mode = lower(char(string(custom.diffusionMode)));
    preset = lower(char(string(custom.referencePreset)));
    if strcmp(preset, 'auto')
        switch mode
            case 'hop'
                preset = 'acsnanohop';
            case 'trap'
                preset = 'stedflcstrap';
            case 'mixed'
                preset = 'acsnanohop';
            otherwise
                preset = 'custom';
        end
    end

    switch preset
        case 'acsnanohop'
            values = struct('freeDiffusionUm2PerS', 0.20, ...
                'hopProbability', 0.01, 'compartmentSizeUm', 0.10);
            resolved = 'acsNanoHop';
        case 'stedfcsnrkhop'
            values = struct('freeDiffusionUm2PerS', 0.80, ...
                'hopProbability', 0.10, 'compartmentSizeUm', 0.08);
            resolved = 'stedFcsNrkHop';
        case 'stedfcsia32hop'
            values = struct('freeDiffusionUm2PerS', 0.80, ...
                'hopProbability', 0.10, 'compartmentSizeUm', 0.15);
            resolved = 'stedFcsIa32Hop';
        case 'stedflcstrap'
            values = struct('freeDiffusionUm2PerS', 0.80, ...
                'trappedDiffusionUm2PerS', 1e-10, ...
                'meanFreeTimeS', 0.030, 'meanTrappedTimeS', 0.010);
            resolved = 'stedFlcsTrap';
        case 'custom'
            values = struct();
            resolved = 'custom';
        otherwise
            error('simulateHopTrapMembraneTracking:BadPreset', ...
                ['referencePreset must be auto, acsNanoHop, ' ...
                'stedFcsNrkHop, stedFcsIa32Hop, stedFlcsTrap, or custom.']);
    end

    names = fieldnames(values);
    for k = 1:numel(names)
        name = names{k};
        if any(strcmp(usingDefaults, name))
            custom.(name) = values.(name);
        end
    end
    custom.referencePresetResolved = resolved;
end
