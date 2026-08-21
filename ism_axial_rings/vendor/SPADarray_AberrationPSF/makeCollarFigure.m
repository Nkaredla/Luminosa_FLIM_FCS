function makeCollarFigure(outFile)
%--------------------------------------------------------------------------
% makeCollarFigure
%
% PURPOSE
%   Generate an illustrative figure for correction-collar calibration.
%
% FIGURE CONTENT
%   Left panel:
%       synthetic Zernike coefficients versus collar setting, along with a
%       quadratic fit
%
%   Right panel:
%       schematic workflow for experimental collar calibration
%
% IMPORTANT
%   This figure is illustrative only. It does not use measured data.
%--------------------------------------------------------------------------

    % Example collar settings [mm]
    collar = [0.13; 0.15; 0.17; 0.19; 0.21];

    % Center collar values about nominal reference
    u = collar - 0.17;

    % Fix RNG for reproducibility
    rng(4);

    % Synthetic coefficient trends versus collar
    coeffs = [...
        -0.12*u + 0.50*u.^2, ...
        0.00*u + 0.06*u.^2, ...
        0.03*u - 0.03*u.^2, ...
        0.02*u, ...
        -0.01*u, ...
        0.35*u + 1.10*u.^2];

    % Add small random scatter to mimic fit noise
    coeffs = coeffs + 0.01*randn(size(coeffs));

    % Fit quadratic collar model
    model = fitQuadraticCollarModel(collar, coeffs, 1e-6);

    % Dense collar values for smooth curves
    dense = linspace(min(collar), max(collar), 200).';

    % Predict fitted trajectories
    pred = predictCollar(model, dense);

    % Create figure
    figure('Color','w','Position',[100 100 980 360]);

    %---------------- Left panel: coefficient trends ----------------%
    subplot(1,2,1);
    hold on;

    cols = lines(6);
    names = {'defocus','astig_x','astig_y','coma_x','coma_y','spherical'};

    for j = 1:6
        plot(collar, coeffs(:,j), 'o', 'Color', cols(j,:), 'LineWidth', 1.2, 'MarkerFaceColor', cols(j,:));
        plot(dense, pred(:,j), '-', 'Color', cols(j,:), 'LineWidth', 1.5);
    end

    grid on;
    xlabel('collar setting (mm)');
    ylabel('coefficient (waves RMS)');
    title('Illustrative fitted Zernike trajectory');
    legend(names, 'Location', 'eastoutside', 'Box', 'off');

    %---------------- Right panel: workflow cartoon ----------------%
    subplot(1,2,2);
    axis off;
    xlim([0 1]);
    ylim([0 1]);

    % Boxes and labels
    text(0.08, 0.82, 'Bead scan', 'FontSize', 12, 'HorizontalAlignment', 'center');
    rectangle('Position', [0.03 0.73 0.10 0.08], 'Curvature', 0.05);

    text(0.30, 0.82, 'Fit raw detector-channel stack', 'FontSize', 12, 'HorizontalAlignment', 'center');
    rectangle('Position', [0.22 0.73 0.16 0.08], 'Curvature', 0.05);

    text(0.58, 0.82, 'Estimate Zernike vector', 'FontSize', 12, 'HorizontalAlignment', 'center');
    rectangle('Position', [0.49 0.73 0.18 0.08], 'Curvature', 0.05);

    text(0.84, 0.82, 'Fit a(c)', 'FontSize', 12, 'HorizontalAlignment', 'center');
    rectangle('Position', [0.78 0.73 0.10 0.08], 'Curvature', 0.05);

    % Arrows between stages
    annotation('arrow', [0.24 0.36], [0.74 0.74]);
    annotation('arrow', [0.48 0.58], [0.74 0.74]);
    annotation('arrow', [0.72 0.80], [0.74 0.74]);

    % Additional descriptive text
    text(0.17, 0.58, 'repeat for multiple', 'FontSize', 11);
    text(0.17, 0.51, 'collar settings c_i', 'FontSize', 11);

    text(0.63, 0.58, 'quadratic model:', 'FontSize', 11);
    text(0.57, 0.49, 'a(c) = a_0 + g_1(c-c_0) + g_2(c-c_0)^2', 'FontSize', 11, 'Interpreter', 'tex');

    text(0.53, 0.30, 'Insert experimental fitted points here', 'FontSize', 12, 'FontAngle', 'italic');
    title('Experimental workflow');

    exportgraphics(gcf, outFile, 'Resolution', 180);
    close(gcf);
end
