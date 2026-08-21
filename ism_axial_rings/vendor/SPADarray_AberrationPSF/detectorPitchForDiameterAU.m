function [detPitchUm, info] = detectorPitchForDiameterAU(layoutName, diameterAU, NA, lamEmUm, detFillRatio)
%DETECTORPITCHFORDIAMETERAU Detector pitch giving a target array diameter in Airy units.
%
%   [detPitchUm, info] = detectorPitchForDiameterAU(layoutName, diameterAU, NA, lamEmUm, detFillRatio)
%
%   Returns the sample-equivalent detector pitch [um] such that the full
%   active-area extent (edge-to-edge, along the widest direction) of the
%   named honeycomb layout equals diameterAU Airy units, using the pinhole
%   convention 1 AU = 1.22*lamEm/NA. This is the knob for objective/relay
%   magnification: the higher the microscope magnification, the smaller the
%   projected (sample-equivalent) pitch and the smaller the array footprint
%   in AU. detFillRatio is detSize/detPitch (default 1, touching hex cells).
%
%   The extent is measured edge-to-edge, i.e. the center-to-center span of
%   the widest row plus one pixel width (half a pixel beyond each outer
%   center). For point-up hex pixels the horizontal flat-to-flat width is
%   detSize and the vertical vertex-to-vertex height is 2*detSize/sqrt(3).
%
%   info fields:
%     airyUnitUm       1 AU in sample space [um]
%     targetDiameterUm target edge-to-edge extent [um]
%     extentInPitches  edge-to-edge extent at unit pitch (pitch multiplier)
%     diameterAU       requested diameter [AU]
%     layoutName       layout used

    if nargin < 5 || isempty(detFillRatio)
        detFillRatio = 1.0;
    end

    airyUnitUm = 1.22 * lamEmUm / NA;
    targetDiameterUm = diameterAU * airyUnitUm;

    % Unit-pitch layout; measure the edge-to-edge active-area extent.
    detXY1 = detectorLayout(layoutName, 1);
    xCenterSpan = max(detXY1(:,1)) - min(detXY1(:,1));
    yCenterSpan = max(detXY1(:,2)) - min(detXY1(:,2));

    detSize1   = detFillRatio;                 % hex flat-to-flat width at unit pitch
    hexHeight1 = 2 * detSize1 / sqrt(3);       % point-up hex vertex-to-vertex height
    spanX = xCenterSpan + detSize1;            % add half a pixel width on each side
    spanY = yCenterSpan + hexHeight1;
    extentInPitches = max(spanX, spanY);

    detPitchUm = targetDiameterUm / extentInPitches;

    info.airyUnitUm       = airyUnitUm;
    info.targetDiameterUm = targetDiameterUm;
    info.extentInPitches  = extentInPitches;
    info.diameterAU       = diameterAU;
    info.layoutName       = char(layoutName);
end
