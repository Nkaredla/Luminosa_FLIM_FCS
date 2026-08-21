function makeRetrievalFigure(sim, outFile)
%--------------------------------------------------------------------------
% makeRetrievalFigure
%
% PURPOSE
%   Validate single-mode recovery from synthetic noisy detector-channel data.
%
% FIGURE CONTENT
%   Bar chart comparing:
%       - true aberration amplitude
%       - retrieved aberration amplitude
%
% NOTES
%   Even-symmetry modes are plotted in magnitude because sign ambiguity can
%   occur in intensity-only retrieval.
%--------------------------------------------------------------------------

    % Modes tested
    modes = {'defocus','astig_x','coma_x','spherical'};

    % True amplitude used in each test
    amp = 0.03;

    % Allocate retrieved amplitudes
    ests = zeros(size(modes));

    % Run single-mode retrieval for each mode
    for k = 1:numel(modes)
        res = singleModeFit(sim, modes{k}, amp);

        % Take magnitude for modes that can be sign-ambiguous
        if ismember(modes{k}, {'defocus','astig_x','spherical'})
            ests(k) = abs(res.estAmp);
        else
            ests(k) = res.estAmp;
        end
    end

    % Make bar plot
    figure('Color','w','Position',[100 100 650 380]);
    xs = 1:numel(modes);

    bar(xs-0.17, repmat(amp, size(xs)), 0.34, 'FaceColor', [0.85 0.2 0.2]);
    hold on;
    bar(xs+0.17, ests, 0.34, 'FaceColor', [0.2 0.3 0.85]);

    set(gca, 'XTick', xs, 'XTickLabel', {'defocus','astig. x','coma x','spherical'});
    ylabel('aberration amplitude (waves RMS)');
    title(sprintf('Single-mode recovery from the raw %d-channel stack', size(sim.detXY,1)));
    legend({'true','retrieved'}, 'Box', 'off', 'Location', 'northwest');
    grid on;

    exportgraphics(gcf, outFile, 'Resolution', 180);
    close(gcf);
end
