function out = showDetectorIntensityVsAberration93(varargin)
%SHOWDETECTORINTENSITYVSABERRATION93 Aberration-mode intensity sweep, 93-pixel array.
%
%   out = showDetectorIntensityVsAberration93();
%   out = showDetectorIntensityVsAberration93('amplitudeWaves', 0.2);
%   out = showDetectorIntensityVsAberration93('modes', {'defocus','spherical'});
%
%   Identical physics and options to showDetectorIntensityVsAberration.m
%   (single on-axis point emitter, vectorial Richards-Wolf detection PSF
%   via psfBessel, explicit finite-pixel integration via
%   detectorCollectionEfficiencyExplicit), but run on the 93-channel
%   honeycomb array (detectorLayout('honeycomb93', ...)) instead of the
%   23-channel array. 'defocus' is included in the default mode list, so
%   this also covers a defocus sweep across amplitudeWaves/zDefocusUm.
%
%   All name-value options are passed through to
%   showDetectorIntensityVsAberration; see that function's help for the
%   full list (NA, lamEm, nMedium, amplitudeWaves, zDefocusUm, modes, ...).
%
%   By default the 93-array pixel pitch is set (via the magnification) so
%   the array footprint spans 'arrayDiameterAU' Airy units (default 1.7),
%   using the same NA/lamEm the sweep runs at (defaults 1.2 / 0.520 um,
%   matching showDetectorIntensityVsAberration). Pass an explicit
%   'detPitchUm' to pin the pitch instead. 'detXY' is not accepted here
%   since the layout is fixed to the 93-array.

    addpath(fileparts(mfilename('fullpath')));

    args = varargin;
    detPitchUm       = getNameValue(args, 'detPitchUm', []);
    arrayDiameterAU  = getNameValue(args, 'arrayDiameterAU', 1.7);
    % NA / lamEm / detFillRatio default to the same values
    % showDetectorIntensityVsAberration uses, so the AU footprint is
    % computed against the optics the sweep actually runs at.
    NA           = getNameValue(args, 'NA', 1.2);
    lamEm        = getNameValue(args, 'lamEm', 0.520);
    detFillRatio = getNameValue(args, 'detFillRatio', 1.0);
    if any(strcmpi(args(1:2:max(numel(args)-1,0)), 'detXY'))
        error('showDetectorIntensityVsAberration93:FixedLayout', ...
            '''detXY'' is fixed to the 93-pixel honeycomb layout in this wrapper; call showDetectorIntensityVsAberration directly to supply custom detector centers.');
    end

    if isempty(detPitchUm)
        detPitchUm = detectorPitchForDiameterAU('honeycomb93', ...
            arrayDiameterAU, NA, lamEm, detFillRatio);
    end
    args = stripNameValue(args, {'detPitchUm', 'arrayDiameterAU'});

    detXY93 = detectorLayout('honeycomb93', detPitchUm);
    out = showDetectorIntensityVsAberration('detXY', detXY93, args{:});
end

% ----------------------------------------------------------------------------
function val = getNameValue(args, name, default)
    val = default;
    for k = 1:2:numel(args)-1
        if strcmpi(args{k}, name)
            val = args{k+1};
        end
    end
end

% ----------------------------------------------------------------------------
function args = stripNameValue(args, names)
    keep = true(1, numel(args));
    for k = 1:2:numel(args)-1
        if any(strcmpi(args{k}, names))
            keep(k:k+1) = false;
        end
    end
    args = args(keep);
end
