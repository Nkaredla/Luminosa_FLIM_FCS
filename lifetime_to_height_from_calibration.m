function heightNm = lifetime_to_height_from_calibration(lifetimeNs, calibOrFile)
% LIFETIME_TO_HEIGHT_FROM_CALIBRATION
% Convert lifetime [ns] to height [nm] using a calibration struct or MAT file.

    if nargin < 2 || isempty(calibOrFile)
        calib = load_lifetime_height_calibration('calibrationCurve.mat');
    elseif ischar(calibOrFile) || isstring(calibOrFile)
        calib = load_lifetime_height_calibration(char(calibOrFile));
    else
        calib = calibOrFile;
    end

    heightNm = nan(size(lifetimeNs), 'double');
    valid = isfinite(lifetimeNs);
    if ~any(valid(:))
        return;
    end

    heightNm(valid) = interp1( ...
        calib.lifetimeNs, calib.heightNm, double(lifetimeNs(valid)), ...
        calib.method, NaN);
end
