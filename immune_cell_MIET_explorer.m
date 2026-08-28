function fig = immune_cell_MIET_explorer(analysisMat)
%IMMUNE_CELL_MIET_EXPLORER Click a pixel or region; see its TCSPC and fits.
%
% fig = immune_cell_MIET_explorer()                 % pick a MAT interactively
% fig = immune_cell_MIET_explorer(analysisMat)      % path to an analysis MAT
% fig = immune_cell_MIET_explorer(acquisitionFolder)
%
% WHAT IT IS FOR
%
% The pipeline decides per pixel between a fixed-SLB, a biexponential and a
% triexponential model, and earlier work in this project established that the
% likelihood alone is a poor guide: the false three-component rate RISES with
% photon count, from 9% to 54%. The maps therefore cannot be taken at face
% value, and the only way to judge a pixel is to look at its decay and at the
% residuals of both models. This window makes that a click.
%
% HOW TO USE IT
%
%   - choose any map from the dropdown: red or blue mean FLIM, the ISM
%     intensities, the fitted lifetime and probability maps, or the 2x2 and 4x4
%     binned versions of those
%   - click a pixel, or type row/column, and set the region size to average a
%     small patch (1x1 up to 15x15)
%   - the lower panels show the pooled TCSPC with BOTH the biexponential and
%     triexponential fits, and their residuals
%
% Changing the map does not refit anything - the decay panel keeps showing the
% region you selected while you flip between maps looking at it.
%
% THE RESIDUALS ARE THE POINT, NOT THE CHI-SQUARE
%
% Three exponentials can only reduce the residual, so a lower chi-square is not
% evidence for a third component. What distinguishes the models is STRUCTURE in
% the biexponential residual. If that residual is flat within noise, the third
% component is not needed for that pixel whatever the maps say. The window
% reports both chi-squares without interpreting them and draws both residuals on
% the same axes so they can be compared directly.
%
% WHERE THE PHOTONS COME FROM
%
% Always the NATIVE per-pixel cube (tcspc_pix), whichever map is displayed. The
% binned maps are anchored sliding sums of that same cube, so pooling the native
% cube over the selected region reproduces them exactly and there is only one
% code path that can be wrong. The cube is read one region at a time with
% matfile rather than loaded whole: it is about 113 MB, and this machine runs
% with C: chronically full, which surfaces as MATLAB "Out of memory" rather than
% as a disk error.
%
% Requires the analysis to have been saved with saveTcspcPix = true.

    if nargin < 1 || isempty(analysisMat)
        [name, folder] = uigetfile({'*.mat', 'Analysis MAT'}, ...
            'Select an immune_cell_MIET analysis MAT');
        if isequal(name, 0); fig = []; return; end
        analysisMat = fullfile(folder, name);
    end
    analysisMat = char(analysisMat);
    if isfolder(analysisMat)
        found = dir(fullfile(analysisMat, '**', ...
            'immune_cell_MIET_640nm_red_analysis.mat'));
        if isempty(found)
            error('immune_cell_MIET_explorer:NoMat', ...
                'No analysis MAT found below %s', analysisMat);
        end
        analysisMat = fullfile(found(1).folder, found(1).name);
    end
    if ~isfile(analysisMat)
        error('immune_cell_MIET_explorer:NoFile', ...
            'Not a file: %s', analysisMat);
    end

    fprintf('immune_cell_MIET_explorer: loading %s\n', analysisMat);
    stored = whos('-file', analysisMat);
    names = {stored.name};

    % A biexponential result file holds "out", not "result", and cannot be
    % browsed directly. Rather than error at the user for opening the
    % obvious-looking file, accept it: use the companion written beside it, or
    % create one. Two files sitting in a folder where only one of them opens is
    % a design problem, not a user mistake.
    if ~ismember('result', names) && ismember('out', names)
        companion = strrep(analysisMat, '.mat', '_explorer.mat');
        if isfile(companion)
            fprintf(['  this is a biexp result; opening its explorer ' ...
                'companion instead:\n  %s\n'], companion);
        else
            fprintf(['  this is a biexp result with no companion; building ' ...
                'one now\n']);
            if exist('export_biexp_for_explorer', 'file') ~= 2
                error('immune_cell_MIET_explorer:NoExporter', ...
                    ['%s holds a biexp result and needs ' ...
                     'export_biexp_for_explorer on the path to be viewed.'], ...
                    analysisMat);
            end
            companion = export_biexp_for_explorer(analysisMat);
        end
        analysisMat = companion;
        stored = whos('-file', analysisMat);
        names = {stored.name};
    end

    if ~ismember('result', names)
        error('immune_cell_MIET_explorer:NoResult', ...
            ['This MAT has no "result" variable: %s\n' ...
             'Expected an analysis MAT from immune_cell_MIET, or a biexp ' ...
             'result (which holds "out").'], analysisMat);
    end
    loaded = load(analysisMat, 'result');
    result = loaded.result;

    cubeVariable = '';
    for candidate = {'tcspc_pix', 'tcspcPix'}
        if ismember(candidate{1}, names); cubeVariable = candidate{1}; break; end
    end
    % A biexp explorer companion deliberately does NOT carry the cube - it is
    % 113 MB per acquisition and this machine runs with C: chronically full, so
    % copying it per binning would cost gigabytes. It records where the cube
    % lives instead, and the decay panel reads it from there through matfile
    % exactly as it would from an analysis MAT.
    cubeFile = analysisMat;
    if isempty(cubeVariable)
        source = immune_cell_MIET_explorer_field(result, ...
            'explorerSource.analysisMat');
        if ~isempty(source)
            source = char(source);
            if isfile(source)
                sourceVars = {whos('-file', source).name};
                for candidate = {'tcspc_pix', 'tcspcPix'}
                    if ismember(candidate{1}, sourceVars)
                        cubeVariable = candidate{1};
                        cubeFile = source;
                        break;
                    end
                end
                if ~isempty(cubeVariable)
                    fprintf('  TCSPC cube read from the source analysis:\n    %s\n', ...
                        cubeFile);
                end
            else
                fprintf(['  NOTE: the recorded source analysis is missing, so ' ...
                    'no decay can be shown:\n    %s\n'], source);
            end
        end
    end
    if isempty(cubeVariable)
        error('immune_cell_MIET_explorer:NoCube', ...
            ['This MAT has no per-pixel TCSPC cube, so no decay can be ' ...
             'shown. Re-run the analysis with cfg.saveTcspcPix = true.']);
    end

    dtNs = immune_cell_MIET_explorer_field(result, 'bayesian.compact.dtNs');
    if isempty(dtNs)
        dtNs = immune_cell_MIET_explorer_field(result, 'channel.dtNs');
    end
    periodNs = immune_cell_MIET_explorer_field(result, ...
        'bayesian.compact.pulsePeriodNs');
    irf = immune_cell_MIET_explorer_field(result, 'irf.curve');
    if isempty(dtNs) || isempty(periodNs) || isempty(irf)
        error('immune_cell_MIET_explorer:NoTiming', ...
            ['The result is missing dtNs, pulsePeriodNs or the IRF, so the ' ...
             'decay cannot be fitted.']);
    end

    catalogue = immune_cell_MIET_explorer_maps(result);
    imageSize = size(catalogue(1).data);
    fprintf('  %d map(s) available, image %d x %d, dt %.4f ns, period %.2f ns\n', ...
        numel(catalogue), imageSize(1), imageSize(2), dtNs, periodNs);

    % ---- layout ---------------------------------------------------------
    fig = uifigure('Name', sprintf('immune_cell_MIET explorer - %s', ...
        analysisMat), 'Position', [80 80 1500 860]);
    outer = uigridlayout(fig, [2 2]);
    outer.RowHeight = {'1x', '1x'};
    outer.ColumnWidth = {'1.1x', '1x'};

    mapPanel = uipanel(outer, 'Title', 'map');
    mapPanel.Layout.Row = [1 2]; mapPanel.Layout.Column = 1;
    mapGrid = uigridlayout(mapPanel, [2 1]);
    mapGrid.RowHeight = {34, '1x'};

    bar = uigridlayout(mapGrid, [1 8]);
    bar.ColumnWidth = {60, '1x', 45, 60, 45, 60, 70, 70};
    bar.Padding = [2 2 2 2];

    uilabel(bar, 'Text', 'map');
    mapMenu = uidropdown(bar, 'Items', {catalogue.label}, ...
        'Value', catalogue(1).label);
    uilabel(bar, 'Text', 'row');
    rowBox = uispinner(bar, 'Limits', [1 imageSize(1)], ...
        'Value', round(imageSize(1) / 2), 'RoundFractionalValues', 'on');
    uilabel(bar, 'Text', 'col');
    colBox = uispinner(bar, 'Limits', [1 imageSize(2)], ...
        'Value', round(imageSize(2) / 2), 'RoundFractionalValues', 'on');
    uilabel(bar, 'Text', 'region');
    regionSize = uispinner(bar, 'Limits', [1 15], 'Value', 1, ...
        'Step', 2, 'RoundFractionalValues', 'on');

    mapAxes = uiaxes(mapGrid);

    decayPanel = uipanel(outer, 'Title', 'pooled TCSPC with both fits');
    decayPanel.Layout.Row = 1; decayPanel.Layout.Column = 2;
    decayAxes = uiaxes(decayPanel);
    decayAxes.Units = 'normalized';
    decayAxes.Position = [0.10 0.14 0.86 0.80];

    lower = uigridlayout(outer, [2 1]);
    lower.Layout.Row = 2; lower.Layout.Column = 2;
    lower.RowHeight = {'1.4x', '1x'};

    residualPanel = uipanel(lower, 'Title', ...
        'residuals - look for STRUCTURE in the 2-exp trace');
    residualAxes = uiaxes(residualPanel);
    residualAxes.Units = 'normalized';
    residualAxes.Position = [0.10 0.18 0.86 0.74];

    readoutPanel = uipanel(lower, 'Title', 'fits');
    readoutGrid = uigridlayout(readoutPanel, [2 1]);
    readoutGrid.RowHeight = {'1x', 22};
    readout = uitextarea(readoutGrid, 'Editable', 'off', 'Value', ...
        {'Click the map, or set row/col, to pool a region and fit it.'});
    status = uilabel(readoutGrid, 'Text', 'ready');

    % ---- state ----------------------------------------------------------
    state = struct( ...
        'analysisMat', analysisMat, 'cubeVariable', cubeVariable, ...
        'cube', matfile(cubeFile), 'catalogue', catalogue, ...
        'imageSize', imageSize, 'dtNs', double(dtNs), ...
        'periodNs', double(periodNs), 'irf', double(irf), ...
        'mapMenu', mapMenu, 'mapAxes', mapAxes, 'decayAxes', decayAxes, ...
        'residualAxes', residualAxes, 'readout', readout, ...
        'status', status, 'rowBox', rowBox, 'colBox', colBox, ...
        'regionSize', regionSize, 'marker', gobjects(1));
    fig.UserData = state;

    mapMenu.ValueChangedFcn = @(src, evt) immune_cell_MIET_explorer_draw(fig);
    rowBox.ValueChangedFcn = @(src, evt) immune_cell_MIET_explorer_select( ...
        fig, rowBox.Value, colBox.Value);
    colBox.ValueChangedFcn = @(src, evt) immune_cell_MIET_explorer_select( ...
        fig, rowBox.Value, colBox.Value);
    regionSize.ValueChangedFcn = @(src, evt) immune_cell_MIET_explorer_select( ...
        fig, rowBox.Value, colBox.Value);

    immune_cell_MIET_explorer_draw(fig);
    immune_cell_MIET_explorer_select(fig, rowBox.Value, colBox.Value);
    fprintf('  ready - click the map to select a region\n');
end
