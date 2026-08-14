function [floatMaps, modelOrderMap, posteriorMap] = insertPatch( ...
        floatMaps, modelOrderMap, posteriorMap, floatPatch, modelPatch, ...
        posteriorPatch, xRange, yRange, regionMask)
%INSERTPATCH Insert one cropped FLIM result without erasing overlapping boxes.

    if isempty(xRange)
        return;
    end
    croppedMask = regionMask(xRange, yRange);
    for layer = 1:size(floatMaps, 3)
        target = floatMaps(:, :, layer);
        targetCrop = target(xRange, yRange);
        source = floatPatch(:, :, layer);
        targetCrop(croppedMask) = source(croppedMask);
        target(xRange, yRange) = targetCrop;
        floatMaps(:, :, layer) = target;
    end
    modelCrop = modelOrderMap(xRange, yRange);
    modelCrop(croppedMask) = modelPatch(croppedMask);
    modelOrderMap(xRange, yRange) = modelCrop;
    posteriorCrop = posteriorMap(xRange, yRange);
    posteriorCrop(croppedMask) = posteriorPatch(croppedMask);
    posteriorMap(xRange, yRange) = posteriorCrop;
end
