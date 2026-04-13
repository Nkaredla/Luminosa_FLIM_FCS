function Luminosa_MIET_Analysis(name, IRF, dye, layers, mic, flag, flagload, pic)
% Local MIET analysis wrapper that keeps upstream vendored code untouched.

    if nargin < 8
        pic = [];
    end
    if nargin < 7
        flagload = [];
    end

    dipoleOrientation = 'fast_rotating';
    if nargin >= 3 && isstruct(dye) && isfield(dye, 'DipoleOrientation')
        dipoleOrientation = luminosa_miet_normalize_orientation(dye.DipoleOrientation);
    end

    if strcmpi(flag, '3D_MIET')
        lamem = dye.lamem;
        lamex = dye.lamex;
        qy = dye.qy;
        tauFree = dye.tau_free;
        curveType = dye.CurveType;

        n0 = layers.n0;
        n1 = layers.n1;
        n = layers.n;
        d0 = layers.d0;
        d1 = layers.d1;
        d = layers.d;

        NA = mic.NA;
        focpos = mic.focpos;
        if isfield(mic, 'pattern')
            pattern = mic.pattern;
        else
            pattern = [];
        end
        dirname = 'W:\Arindam\20200721_Nanodisc_TestMeasurement.sptw\MIET_analysis\';
        fnames = dir([dirname '*.mat']);

        for idx = 1:numel(fnames)
            fileName = fnames(idx).name;
            Luminosa_MIET_Analysis([dirname fileName], [], dye, layers, mic, flag, 1, 1);
        end

        [head, FLIM, INT, SM] = fastSM_FLIM_MLE(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic); %#ok<NASGU,ASGLU>

    elseif strcmpi(flag, 'MIET') || strcmpi(flag, 'ROI MIET')
        lamem = dye.lamem;
        lamex = dye.lamex;
        qy = dye.qy;
        tauFree = dye.tau_free;
        curveType = dye.CurveType;

        n0 = layers.n0;
        n1 = layers.n1;
        n = layers.n;
        d0 = layers.d0;
        d1 = layers.d1;
        d = layers.d;
        [lamemCalc, d0Calc, dCalc, d1Calc] = convertCalcInputs(lamem, d0, d, d1);

        NA = mic.NA;
        focpos = mic.focpos;
        if isfield(mic, 'pattern')
            pattern = mic.pattern;
        else
            pattern = [];
        end
        if strcmpi(flag, 'MIET')
            if flagload
                load([name(1:end-4), 'GOLD']);
            else
                [head, FLIM, INT, SM] = fastSM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic); %#ok<ASGLU>
            end
        else
            [head, FLIM, INT] = ROI_FLIM(name, IRF, pic); %#ok<ASGLU>
        end

    elseif strcmpi(flag, 'pattern match') || strcmpi(flag, 'lifetime') || strcmpi(flag, 'random') || strcmpi(flag, 'combined')
        lamex = dye.lamex;
        n0 = layers.n0;
        n1 = layers.n1;
        n = layers.n;
        d0 = layers.d0;
        d1 = layers.d1;
        d = layers.d;

        NA = mic.NA;
        focpos = mic.focpos;
        if isfield(mic, 'pattern')
            pattern = mic.pattern;
        else
            pattern = [];
        end
        [head, FLIM, INT, SM] = fastSM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic); %#ok<ASGLU>

    elseif strcmpi(flag, 'intensity')
        [head, FLIM, INT, SM] = SM_FLIM(name, IRF, [], [], [], [], [], [], [], [], [], flag, 0, pic); %#ok<ASGLU>

    elseif strcmpi(flag, 'show MIET curves')
        lamem = dye.lamem;
        qy = dye.qy;
        tauFree = dye.tau_free;
        curveType = dye.CurveType;

        n0 = layers.n0;
        n1 = layers.n1;
        n = layers.n;
        d0 = layers.d0;
        d1 = layers.d1;
        d = layers.d;
        [lamemCalc, d0Calc, dCalc, d1Calc] = convertCalcInputs(lamem, d0, d, d1);

        [z, lifecurve] = Luminosa_MIET_calc(NaN, lamemCalc, n0, n, n1, d0Calc, dCalc, d1Calc, qy, tauFree, [], curveType, dipoleOrientation);
        valid = ~isnan(lifecurve);
        figure;
        set(gcf, 'Name', sprintf('MIET Calibration Curve (%s)', strrep(dipoleOrientation, '_', ' ')), 'NumberTitle', 'off');
        plot(z(valid), lifecurve(valid), 'LineWidth', 2);
        xlabel('distance from surface (nm)');
        ylabel('lifetime (ns)');
        grid on;
        return

    elseif strcmpi(flag, 'ROI FLIM')
        [head, FLIM, INT] = ROI_FLIM(name, IRF, pic); %#ok<ASGLU>
    else
        error('Unsupported MIET analysis flag: %s', flag);
    end

    if strcmpi(flag, 'MIET')
        display('Determining axial distances...')
        maxch = size(INT.tag,3);
        al_res = unique(SM.al_res);
        [z, lifecurve] = Luminosa_MIET_calc(al_res, lamemCalc, n0, n, n1, d0Calc, dCalc, d1Calc, qy, tauFree, [], curveType);
        for ch = 1:maxch
            if iscell(FLIM.life_matav)
                lifetimes = FLIM.life_matav{ch};
            else
                if size(FLIM.life_matav,1) < 4
                    lifetimes = FLIM.life_matav(ch,:);
                else
                    lifetimes = FLIM.life_matav(:);
                end
            end
            if iscell(SM.theta)
                theta = SM.theta{ch};
            else
                if size(SM.theta,1) < 4
                    theta = SM.theta(:);
                else
                    theta = SM.theta(:,ch);
                end
            end
            if all(isnan(theta(:)))
                theta = NaN * lifetimes;
            end
            ax_pos = Axial_pos(theta, lifetimes, z, lifecurve);
            if iscell(INT.field)
                field = INT.field{ch};
            else
                field = INT.field;
            end
            posimm = zeros(size(field));
            for idx = 1:size(field,3)
                posimm(:,:,idx,ch) = field(:,:,idx) .* ax_pos(idx);
            end
            posimm_c = zeros(size(posimm,1), size(posimm,2), maxch);
            for x = 1:size(posimm,1)
                for y = 1:size(posimm,2)
                    posimm_c(x,y,ch) = pmean(posimm(x,y,:,ch));
                end
            end
        end
        if isequal(pic, 1)
            figure
            set(gcf, 'name', 'Height Image', 'NumberTitle', 'off')
            for ch = 1:maxch
                x_cord = head.ImgHdr.X0;
                y_cord = head.ImgHdr.Y0;
                subplot(ch, 2, maxch)
                imagesc(x_cord, y_cord, posimm_c(:,:,ch).')
                colormap('copper');
                colorbar
            end
        end
        save([name(1:end-4) '_miet_analysis.mat'], 'posimm', 'posimm_c', 'dye', 'mic', 'layers', 'z', 'lifecurve');
    end

    if strcmpi(flag, 'ROI MIET')
        display('Determining axial distances...')
        maxch = size(INT.tag,3);
        [lamemCalc, d0Calc, dCalc, d1Calc] = convertCalcInputs(lamem, d0, d, d1);
        [z, lifecurve] = Luminosa_MIET_calc(NaN, lamemCalc, n0, n, n1, d0Calc, dCalc, d1Calc, qy, tauFree, [], curveType, dipoleOrientation);
        limitLTInd = find(isnan(lifecurve), 1);
        if isempty(limitLTInd)
            limitLT = lifecurve(end);
        else
            limitLT = lifecurve(limitLTInd - 1);
        end
        for ch = 1:maxch
            lifetimes = FLIM.life_matav(:,ch);
            coord = INT.coord(:,ch);
            ax_pos = zeros(numel(lifetimes),1);
            heightImage = zeros(size(INT.tag));
            for idx = 1:numel(lifetimes)
                if lifetimes(idx) > limitLT || lifetimes(idx) < lifecurve(1)
                    ax_pos(idx) = NaN;
                else
                    for y = 2:length(lifecurve)
                        if lifecurve(y) >= lifetimes(idx)
                            ax_pos(idx) = z(y-1) + (lifetimes(idx) - lifecurve(y-1)) * ...
                                (z(y) - z(y-1)) / (lifecurve(y) - lifecurve(y-1));
                            break;
                        end
                    end
                end
                for coordIdx = 1:size(coord{idx},1)
                    heightImage(coord{idx}(coordIdx,2), coord{idx}(coordIdx,1)) = ax_pos(idx);
                end
            end
        end
        if isequal(pic, 1)
            figure
            set(gcf, 'name', 'Height Image', 'NumberTitle', 'off')
            for ch = 1:maxch
                x_cord = head.ImgHdr.Xcord;
                y_cord = head.ImgHdr.Ycord;
                subplot(ch, 1, maxch)
                imagesc(x_cord, y_cord, heightImage)
                colormap('hot');
                colorbar
            end
        end
        save([name(1:end-4) '_ROImiet_analysis.mat'], 'heightImage', 'dye', 'mic', 'layers', 'z', 'lifecurve');
    end

    function [lamemOut, d0Out, dOut, d1Out] = convertCalcInputs(lamemIn, d0In, dIn, d1In)
        lamemOut = lamemIn;
        if ~isstruct(lamemOut)
            lamemOut = lamemOut * 1e3;
        end
        d0Out = d0In * 1e3;
        d1Out = d1In * 1e3;
        if ischar(dIn)
            dOut = dIn;
        else
            dOut = dIn * 1e3;
        end
    end
end
