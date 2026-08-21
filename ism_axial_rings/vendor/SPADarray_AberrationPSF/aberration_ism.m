function aberration_ism(outDir)
%--------------------------------------------------------------------------
% generate_all_latex_figures_matlab
%
% PURPOSE
%   Generate all figure assets used by the manuscript LaTeX file.
%
% OUTPUT FILES
%   figure1_combined.png
%   retrieval_validation.png
%   sensitivity_and_distance_symmetric.png
%   collar_calibration.png
%
% DESCRIPTION
%   This script is self-contained. It builds a simplified detector-resolved
%   ISM forward model using:
%       1) a 3D fluorescent bead object,
%       2) aberrated excitation and detection PSFs,
%       3) detector offsets for a 5x5 array,
%       4) Bessel / azimuthal harmonic evaluation of the PSF.
%
%   It then generates:
%       - representative mode signatures,
%       - a single-mode inverse-retrieval validation figure,
%       - a sensitivity / separability figure,
%       - an illustrative correction-collar calibration figure.
%--------------------------------------------------------------------------

    % If the user did not provide an output directory, use 'assets'
    if nargin < 1 || isempty(outDir)
        outDir = 'assets';
    end

    % Create the output folder if it does not already exist
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Build the main simulation parameter structure
    sim = defaultParams();

    % Generate each manuscript figure
%     makeModeFigure(sim, fullfile(outDir, 'figure1_combined.png'));
    makeModeFigureExplicitDetector(sim, fullfile(outDir,'figure1_explicit.png'));
    makeRetrievalFigure(sim, fullfile(outDir, 'retrieval_validation.png'));
    makeSensitivityFigure(sim, fullfile(outDir, 'sensitivity_and_distance_symmetric.png'));
    makeCollarFigure(fullfile(outDir, 'collar_calibration.png'));

    % Report completion to the command window
    fprintf('Wrote assets to %s\n', outDir);
end
