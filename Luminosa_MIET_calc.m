function [z, life] = Luminosa_MIET_calc(al_res, lamem, n0, n, n1, d0, d, d1, qyield, tau_free, fig, curveType, orientationMode)
% Calculate local MIET calibration curves without modifying the vendored code.
% All length inputs are in nm.

    if nargin < 11 || isempty(fig)
        fig = 0;
    else
        fig = 1;
    end
    if nargin < 12 || isempty(curveType)
        curveType = 'minimum';
    end
    if nargin < 13
        orientationMode = '';
    end

    orientationMode = luminosa_miet_normalize_orientation(orientationMode);
    useOrientationMode = ~isempty(orientationMode) && isOrientationPlaceholder(al_res);

    if qyield > 1
        qyield = qyield / 100;
    end

    maxCurve = strcmpi(curveType, 'maximum');
    manualCurve = strcmpi(curveType, 'manual');

    if isstruct(lamem)
        polychrome = true;
        metalsFile = luminosa_miet_vendor_file('metals.mat');
        if isempty(metalsFile)
            error('MIET-GUI support file metals.mat was not found in the vendored MIET tree.');
        end
        metals = load(metalsFile);
        spectrum = load(lamem.SpectrumFile);
        if spectrum(1,1) < 2
            spectrum(:,1) = spectrum(:,1) * 1000;
        end
        wavelSmall = lamem.Wavel_Small;
        wavelLarge = lamem.Wavel_Large;
        spectrum = spectrum(spectrum(:,1) >= wavelSmall,:);
        spectrum = spectrum(spectrum(:,1) <= wavelLarge,:);
        spectrum(:,1) = round(spectrum(:,1) / 5);
        spectrIntensity = accumarray(spectrum(:,1) - spectrum(1,1) + 1, spectrum(:,2));
        spectrum = [unique(spectrum(:,1) * 5) spectrIntensity / sum(spectrIntensity)];
    else
        polychrome = false;
        spectrum = [lamem 1];
        metals = struct([]);
    end

    dVar = 100;
    if ~ischar(d)
        if numel(n1) > 1 || n1 ~= n
            if d < 100
                z = 1:(d - 1) / 100:d - 1;
            else
                z = 1:1:min(d, dVar) - 1;
            end
        else
            z = 1:1:dVar - 1;
            d = dVar;
        end
    else
        z = 0.5;
    end

    n0Backup = n0;
    n1Backup = n1;
    if useOrientationMode || isOrientationPlaceholder(al_res)
        life = zeros(numel(z), 1);
    elseif numel(al_res) == 1
        life = zeros(numel(z), numel(90:-al_res:0));
    else
        life = zeros(numel(z), numel(al_res));
    end

    for wavelCount = 1:size(spectrum, 1)
        lamemNow = spectrum(wavelCount, 1);
        fac = 2 * pi / lamemNow;
        if polychrome
            for nIndex = 1:numel(n0)
                if n0Backup(nIndex) >= 10
                    n0(nIndex) = refrIndex(n0Backup(nIndex), lamemNow, metals);
                end
            end
            for nIndex = 1:numel(n1)
                if n1Backup(nIndex) >= 10
                    n1(nIndex) = refrIndex(n1Backup(nIndex), lamemNow, metals);
                end
            end
        end

        if ~ischar(d)
            [~,~,~,~,qvd,qvu,qpd,qpu] = LifetimeL(fac * z, n0, n, n1, fac * d0, fac * d, fac * d1);
        else
            tmpd = 0.5;
            z = 1:dVar - 1;
            qvd = zeros(size(z));
            qvu = zeros(size(z));
            qpd = zeros(size(z));
            qpu = zeros(size(z));
            for idx = 1:numel(z)
                [~,~,~,~,qvd(idx),qvu(idx),qpd(idx),qpu(idx)] = LifetimeL(fac * tmpd, n0, n, n1, fac * [d0 z(idx)], fac * 0, fac * d1);
            end
        end

        sv = qvu + qvd;
        sp = qpu + qpd;

        if useOrientationMode
            wavelengthLife = orientLifetimeCurve(sv, sp, n, qyield, tau_free, orientationMode);
            wavelengthLife = limitCurve(wavelengthLife, maxCurve, manualCurve);
        elseif isOrientationPlaceholder(al_res)
            wavelengthLife = tau_free ./ (1 - qyield + qyield * (sv.' + 2 .* sp.') / (4 * n));
            wavelengthLife = limitCurve(wavelengthLife, maxCurve, manualCurve);
        else
            if numel(al_res) == 1
                theta = (90:-al_res:0) .* pi ./ 180;
            else
                theta = al_res;
            end
            wavelengthLife = zeros(numel(z), numel(theta));
            for idx = 1:numel(theta)
                sr = sv .* (cos(theta(idx)).^2) + sp .* (sin(theta(idx)).^2);
                wavelengthLife(:,idx) = tau_free ./ ((1 - qyield) + sr / (4/3 * n) * qyield);
                wavelengthLife(:,idx) = limitCurve(wavelengthLife(:,idx), maxCurve, manualCurve);
            end
        end

        life = life + spectrum(wavelCount, 2) * wavelengthLife;
    end

    if fig
        figure;
        set(gcf, 'NumberTitle', 'off');
        if useOrientationMode
            set(gcf, 'Name', sprintf('MIET Calibration Curve (%s)', strrep(orientationMode, '_', ' ')));
            plot(z, life, 'LineWidth', 2);
        elseif isOrientationPlaceholder(al_res)
            set(gcf, 'Name', 'MIET Calibration Curve for Random Orientation');
            plot(z, life, 'LineWidth', 2);
        else
            set(gcf, 'Name', 'MIET Calibration Curves for Various Orientations');
            c = hsv(size(life, 2));
            handles = gobjects(1, size(life, 2));
            labels = cell(1, size(life, 2));
            for idx = 1:size(life, 2)
                handles(idx) = plot(z, life(:,idx), 'Color', c(idx,:), 'LineWidth', 2);
                labels{idx} = sprintf('polar angle = %g^{\\circ}', theta(idx) * 180 / pi);
                hold on;
            end
            legend(handles, labels, 'Location', 'SouthEast');
            hold off;
        end
        xlabel('distance from surface (nm)');
        ylabel('lifetime (ns)');
        axis tight;
    end

    function tf = isOrientationPlaceholder(value)
        tf = isempty(value) || (isscalar(value) && isnumeric(value) && isnan(value));
    end

    function lifeVec = orientLifetimeCurve(svIn, spIn, refrIndex, quantumYield, tauFree, mode)
        switch mode
            case 'vertical'
                sr = svIn(:);
                lifeVec = tauFree ./ ((1 - quantumYield) + quantumYield * sr / (4/3 * refrIndex));
            case 'parallel'
                sr = spIn(:);
                lifeVec = tauFree ./ ((1 - quantumYield) + quantumYield * sr / (4/3 * refrIndex));
            case 'random_fixed'
                thetaGrid = linspace(0, pi/2, 181);
                weights = sin(thetaGrid);
                sr = svIn(:) * (cos(thetaGrid).^2) + spIn(:) * (sin(thetaGrid).^2);
                lifeByTheta = tauFree ./ ((1 - quantumYield) + quantumYield * sr / (4/3 * refrIndex));
                lifeVec = trapz(thetaGrid, lifeByTheta .* weights, 2) ./ trapz(thetaGrid, weights);
            otherwise
                lifeVec = tauFree ./ (1 - quantumYield + quantumYield * (svIn(:) + 2 .* spIn(:)) / (4 * refrIndex));
        end
    end

    function lifeVec = limitCurve(lifeVec, keepMaximum, keepManual)
        if keepManual
            return;
        end
        peak = find(diff(lifeVec) < 0, 1);
        if isempty(peak)
            return;
        end
        if keepMaximum
            lifeVec(peak:end) = NaN;
            return;
        end
        limitLT = min(lifeVec(peak:end));
        limitHeight = find(lifeVec > limitLT, 1);
        if ~isempty(limitHeight)
            lifeVec(limitHeight:end) = NaN;
        end
    end

    function value = refrIndex(material, wavelengthNm, metalsData)
        switch material
            case 10
                value = metalsData.silver(metalsData.wavelength == wavelengthNm);
            case 20
                value = metalsData.gold(metalsData.wavelength == wavelengthNm);
            case 30
                value = metalsData.platinum(metalsData.wavelength == wavelengthNm);
            case 40
                value = metalsData.palladium(metalsData.wavelength == wavelengthNm);
            case 50
                value = metalsData.copper(metalsData.wavelength == wavelengthNm);
            case 60
                value = metalsData.aluminum(metalsData.wavelength == wavelengthNm);
            case 70
                value = metalsData.chromium(metalsData.wavelength == wavelengthNm);
            case 80
                value = metalsData.titan(metalsData.wavelength == wavelengthNm);
            case 90
                value = metalsData.tungsten(metalsData.wavelength == wavelengthNm);
            case 100
                value = metalsData.nickel(metalsData.wavelength == wavelengthNm);
            case 110
                value = metalsData.beryllium(metalsData.wavelength == wavelengthNm);
            case 120
                value = metalsData.ito(metalsData.wavelength == wavelengthNm);
            otherwise
                error('The chosen material does not exist.');
        end
    end
end
