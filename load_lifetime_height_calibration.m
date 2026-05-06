function calib = load_lifetime_height_calibration(matFile)
% LOAD_LIFETIME_HEIGHT_CALIBRATION
% Load a MIET-style calibration MAT file with columns:
%   column 1 = height [nm]
%   column 2 = lifetime [ns]

    if nargin < 1 || isempty(matFile)
        matFile = 'calibrationCurve.mat';
    end

    s = load(matFile);
    curve = [];

    if isfield(s, 'calibrationCurve')
        curve = s.calibrationCurve;
    else
        f = fieldnames(s);
        for k = 1:numel(f)
            v = s.(f{k});
            if isnumeric(v) && ismatrix(v) && all(size(v) >= 2)
                curve = v;
                break;
            end
        end
    end

    if isempty(curve)
        error('No numeric calibration curve found in %s.', matFile);
    end

    curve = double(curve);
    if size(curve, 2) < 2 && size(curve, 1) >= 2
        curve = curve.';
    end
    if size(curve, 2) < 2
        error('Calibration curve in %s must have at least two columns.', matFile);
    end

    heightNm = curve(:, 1);
    lifetimeNs = curve(:, 2);
    valid = isfinite(heightNm) & isfinite(lifetimeNs);
    heightNm = heightNm(valid);
    lifetimeNs = lifetimeNs(valid);

    [lifetimeNs, ord] = sort(lifetimeNs(:), 'ascend');
    heightNm = heightNm(ord);
    [lifetimeNs, uniqIdx] = unique(lifetimeNs, 'stable');
    heightNm = heightNm(uniqIdx);

    if numel(lifetimeNs) < 2
        error('Calibration curve in %s must contain at least two unique lifetime values.', matFile);
    end

    if numel(lifetimeNs) >= 3
        method = 'pchip';
    else
        method = 'linear';
    end

    calib = struct();
    calib.file = matFile;
    calib.heightNm = heightNm(:);
    calib.lifetimeNs = lifetimeNs(:);
    calib.method = method;
    calib.curve = [heightNm(:) lifetimeNs(:)];
end
