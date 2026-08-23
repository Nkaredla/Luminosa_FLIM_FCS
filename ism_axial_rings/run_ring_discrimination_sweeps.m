function out = run_ring_discrimination_sweeps(opts)
%RUN_RING_DISCRIMINATION_SWEEPS Photon/abundance and displacement sweeps.
%
% out = run_ring_discrimination_sweeps()
%
% Two sweeps of study_ring_height_discrimination, sharing one forward-model
% evaluation so the ring weight table is computed once:
%
%   A  how few photons, and how small a third-pool abundance, still support a
%      displaced-versus-coplanar decision at a fixed displacement
%   B  the minimum detectable displacement at a fixed photon count and
%      abundance - the number that decides whether this is useful on the real
%      data, given that the ring pattern saturates above roughly 0.3 um
%
% Sweep B is the one that matters for interpretation. MIET lifetime alone goes
% flat beyond about 0.2 um, so the interesting question is whether the rings
% cover the gap between there and the point where their own pattern saturates.

    thisDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(opts); opts = struct(); end
    defaults = struct( ...
        'repeats', 60, 'restarts', 3, ...
        'fractions', [0.02 0.05 0.10], ...
        'photonTotals', [1e3 3e3 1e4 3e4], ...
        'displacements', [0.08 0.15 0.25 0.40 0.60], ...
        'displacementPhotons', 1e4, 'displacementFraction', 0.10, ...
        'outputDir', thisDir);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
            opts.(names{k}) = defaults.(names{k});
        end
    end

    fprintf('computing the ring weight table once for both sweeps...\n');
    ring = simulate_ring_weight_vs_height(struct('makeFigure', false, ...
        'heightsUm', 0:0.02:1.0, 'nx', 41, 'fovUm', 2.6));
    shared = struct('ringWeight', ring.ringWeight, 'heightsUm', ring.heightsUm);

    fprintf('\n================ SWEEP A: photons and abundance ================\n');
    sweepA = study_ring_height_discrimination(struct( ...
        'ringWeight', shared, 'repeats', opts.repeats, ...
        'restarts', opts.restarts, ...
        'thirdPoolFractions', opts.fractions, ...
        'photonTotals', opts.photonTotals, ...
        'outputDir', opts.outputDir, 'makeFigure', true));

    fprintf('\n============= SWEEP B: minimum detectable displacement =============\n');
    rows = struct([]);
    for k = 1:numel(opts.displacements)
        displacement = opts.displacements(k);
        fprintf('\n  --- third pool at z = %.3f um ---\n', displacement);
        result = study_ring_height_discrimination(struct( ...
            'ringWeight', shared, 'repeats', opts.repeats, ...
            'restarts', opts.restarts, ...
            'displacedHeightUm', displacement, ...
            'thirdPoolFractions', opts.displacementFraction, ...
            'photonTotals', opts.displacementPhotons, ...
            'makeFigure', false, 'outputDir', opts.outputDir));
        entry = table2struct(result.summary(1, :));
        entry.displacedHeightUm = displacement;
        if isempty(rows); rows = entry; else; rows(end+1) = entry; end %#ok<AGROW>
    end
    sweepB = struct2table(rows);

    fprintf('\n  minimum detectable displacement, N=%d, third pool %.0f%% of photons\n', ...
        opts.displacementPhotons, 100 * opts.displacementFraction);
    fprintf('      z(um)   ring power   summed power   ring pattern gap\n');
    for k = 1:height(sweepB)
        fprintf('      %5.3f   %10.2f   %12.2f   %16.4f\n', ...
            sweepB.displacedHeightUm(k), sweepB.ringPower(k), ...
            sweepB.summedPower(k), sweepB.ringPatternGap(k));
    end

    csvFile = fullfile(opts.outputDir, 'ring_displacement_sweep.csv');
    writetable(sweepB, csvFile);
    fprintf('\n  wrote %s\n', csvFile);

    out = struct('sweepA', sweepA, 'sweepB', sweepB, 'opts', opts, ...
        'ringWeight', shared);
end
