function makeSensitivityFigure(sim, outFile)
%--------------------------------------------------------------------------
% makeSensitivityFigure
%
% PURPOSE
%   Quantify how distinguishable different aberration-induced detector stacks
%   are from one another.
%
% FIGURE CONTENT
%   Left panel:
%       distance from the unaberrated stack as mode amplitude increases
%
%   Right panel:
%       pairwise distance matrix between representative mode stacks at
%       0.18 waves RMS
%--------------------------------------------------------------------------

    % Aberration modes to compare
    modes = {'defocus','astig_x','coma_x','spherical'};
    labels = {'none','defocus','astig. x','coma x','spherical'};

    % Amplitudes used in the sensitivity sweep
    amps = linspace(0, 0.18, 7);

    % Baseline stack with no aberration
    base = normalizedStack(sim, struct());

    % Distance curves versus amplitude
    curves = zeros(numel(amps), numel(modes));
    for j = 1:numel(modes)
        for i = 1:numel(amps)
            coeffs = struct(modes{j}, amps(i));
            st = normalizedStack(sim, coeffs);
            curves(i,j) = symDistance(base, st);
        end
    end

    % Representative stacks at 0.18 waves RMS
    stacks = cell(1, numel(labels));
    stacks{1} = base;
    for j = 1:numel(modes)
        stacks{j+1} = normalizedStack(sim, struct(modes{j}, 0.18));
    end

    % Pairwise distance matrix
    D = zeros(numel(labels));
    for i = 1:numel(labels)
        for j = 1:numel(labels)
            D(i,j) = symDistance(stacks{i}, stacks{j});
        end
    end

    % Plot figure
    figure('Color','w','Position',[100 100 900 360]);

    subplot(1,2,1);
    plot(amps, curves, 'LineWidth', 1.7);
    grid on;
    xlabel('mode amplitude (waves RMS)');
    ylabel('symmetric normalized distance');
    title('Sensitivity to aberration amplitude');
    legend({'defocus','astig. x','coma x','spherical'}, 'Location', 'northwest', 'Box', 'off');

    subplot(1,2,2);
    imagesc(D);
    axis image;
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
             'YTick', 1:numel(labels), 'YTickLabel', labels);
    xtickangle(35);
    colorbar;
    title('Pairwise distances at 0.18 waves RMS');

    exportgraphics(gcf, outFile, 'Resolution', 180);
    close(gcf);
end