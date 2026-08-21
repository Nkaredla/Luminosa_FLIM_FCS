function result = fitCenterCoordinateISMWaterAirInterfaceBlurred(stackInput, varargin)
%FITCENTERCOORDINATEISMWATERAIRINTERFACEBLURRED Compare with model-only blur.
%
%   result = fitCenterCoordinateISMWaterAirInterfaceBlurred()
%   result = fitCenterCoordinateISMWaterAirInterfaceBlurred(stackInput)
%
%   This keeps the physical 60x 1.2 NA water/glass/air model and the fixed
%   regular PDA23 detector honeycomb, but broadens only the recovered 3D
%   model used for measured-vs-fitted projection figures. The fitted
%   aberration coefficients and detector-center fit are not blurred.
%
%   Default empirical broadening:
%       sigmaXY = 0.08 um
%       sigmaZ  = 0.15 um

    if nargin < 1
        stackInput = [];
    end

    root = fileparts(fileparts(mfilename('fullpath')));
    outDir = fullfile(root, 'output_matlab', ...
        'water_air_center_coordinate_ism', ...
        'test_20260515_144001_model_blur_xy80nm_z150nm');

    defaultBlurArgs = { ...
        'modelBlurSigmaXYUm', 0.08, ...
        'modelBlurSigmaZUm', 0.15, ...
        'outputPrefix', 'water_air_center_coordinate_model_blur_xy80nm_z150nm', ...
        'reportTitle', ...
            'Water/glass/air center-coordinate ISM fit, model blur XY 80 nm / Z 150 nm', ...
        'outputDir', outDir};

    result = fitCenterCoordinateISMWaterAirInterface( ...
        stackInput, defaultBlurArgs{:}, varargin{:});
end
