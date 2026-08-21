function writeOutputs(outputDir, out)
%WRITEOUTPUTS Save reproducible MAT and tabular hop/trap results.

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    % Graphics handles are session-specific and make MAT files needlessly
    % large; the exported PNG is the persistent figure artifact.
    out.figureSummary = [];
    save(fullfile(outputDir, 'hop_trap_tracking_results.mat'), ...
        'out', '-v7.3');
    if ~isempty(out.localizationTable)
        writetable(out.localizationTable, fullfile(outputDir, ...
            'hop_trap_localizations.csv'));
    end
    if ~isempty(out.trackTable)
        writetable(out.trackTable, fullfile(outputDir, ...
            'hop_trap_tracks.csv'));
    end
    if ~isempty(out.stepTable)
        writetable(out.stepTable, fullfile(outputDir, ...
            'hop_trap_independent_steps.csv'));
    end
    if ~isempty(out.msdTable)
        writetable(out.msdTable, fullfile(outputDir, ...
            'hop_trap_msd.csv'));
    end
end
