function MIET_Analysis(name, IRF, dye, layers, mic, flag, flagload, pic)
%  MIET_Analysis(name, IRF, dye, layers, mic, flag, pic)
%  This function initiates the analysis of the data file (.ht3 format)
%
%  Input Parameters:
%  name              : Name and complete address of the data .ht3 file
%  IRF               : Name and complete address of the IRF .ht3 file
%  dye               : Structure that contains dye properties
%     dye.lamem      : The emission maximum of the the dye (in nm)
%     dye.lamex      : The excitation wavelength of the laser (in ?m)
%     dye.qy         : The quantum yield of the dye in free solution or
%                      medium
%     dye.tau_free   : The free solution lifetime of the dye
%     dye.CurveType  : Calibration curve calculated up to 1st maximum
%                      ('maximum') or to last unambiguous value ('minimum')
%  layers            : Structure that contains the layer properties of the
%                      sample
%     layers.n0      : vector of refractive indices of layers below the
%                      sample from
%                      bottom to top
%     layers.n1      : vector of refractive indices of layers above the
%                      surface from
%                      bottom to top
%     layers.n       : the refractive index of the medium containing the
%                      sample molecules
%     layers.d0      : vector of the thicknesses of the layers below the
%                      sample layer in ?m.
%                      length (d0)= length(n0)-1
%     layers.d1      : d1 is the vector of the thicknesses of the layers
%                      above the sample layer in ?m.
%                      length (d1)= length(n1)-1
%     layers.d       : thickness of the layer containing the sample
%                      molecules in ?m.
%  mic               : Structure that contains the microscope properties
%     mic.NA         : Numerical Aperture of the objective used
%     mic.focpos     : The position of the focus with respect to the
%                      coverslip.
%     mic.pattern    : Mentions the polarization of the laser used for the
%                      experiment.
% flag               : Is the task flag for the user to mention what to
%                      calculate
% flag = 'intensity'
%                      Tells the program to compute only the intensity
%                      image
% flag = 'pattern match'
%                      Tells the program to compute the intensity image and
%                      then match with computed patterns based on
%                      mic.pattern
% flag = 'lifetime'
%                      Tells the program to identify the patterns and fit
%                      the lifetime of the photons for the identified
%                      molecules pixels.
% flag = 'MIET'
%                      Tells the program to perform complete MIET analysis
% flag = 'show MIET curves'
%                      Commands the program to present only the calibration
%                      lifetime vs axial position curves for various polar
%                      angles.
% flag = 'ROI FLIM'
%                      Tells the program to compute the intensity image,
%                      and then acquire user input to obtain the region in
%                      the image where the lifetime analysis must be
%                      performed.
% flag = 'ROI MIET'
%                      Tells the program to do MIET analysis on top of the
%                      ROI FLIM process.
% flagauto             Tells the program to select from the fitted patterns
%                      itself so that the program is fully automatic.

% (c) Narain Karedla & Daja Ruhlandt, 2014
if  strcmpi(flag,'3D_MIET')
    lamem    = dye.lamem;   % nm
    lamex    = dye.lamex;   % ?m
    qy       = dye.qy;
    tau_free = dye.tau_free;
    curveType= dye.CurveType;
    
    n0       = layers.n0;
    n1       = layers.n1;
    n        = layers.n;
    d0       = layers.d0;   % ?m
    d1       = layers.d1;   % ?m
    d        = layers.d;    % ?m
    
    
    NA       = mic.NA;
    focpos   = mic.focpos;
    if isfield(mic,'pattern')
        pattern = mic.pattern;
    else
        pattern=[];
    end
    dirname = 'W:\Arindam\20200721_Nanodisc_TestMeasurement.sptw\MIET_analysis\';
    fnames = dir([dirname '*.mat']);
    
    for j=1:numel(fnames)
        name = fnames(j).name
        %     IRF='U:\Narain\140305\Point_011.ht3';
        IRF=[];
        MIET_Analysis([dirname name], IRF, dye, layers, mic, flag,1,1);
        %MIET_Analysis(fnames, IRF, dye, layers, mic, flag,1,1);
    end
    
    [head, FLIM, INT, SM]= fastSM_FLIM_MLE(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic);
    
elseif  strcmpi(flag,'MIET') || strcmpi(flag,'ROI MIET')
    lamem    = dye.lamem;   % nm
    lamex    = dye.lamex;   % ?m
    qy       = dye.qy;
    tau_free = dye.tau_free;
    curveType= dye.CurveType;
    
    n0       = layers.n0;
    n1       = layers.n1;
    n        = layers.n;
    d0       = layers.d0;   % ?m
    d1       = layers.d1;   % ?m
    d        = layers.d;    % ?m
    
    NA       = mic.NA;
    focpos   = mic.focpos;
    
    if isfield(mic,'pattern')
        pattern = mic.pattern;
    else
        pattern=[];
    end
    if strcmpi(flag,'MIET')
        if flagload
            load([name(1:end-4),'GOLD']);
        else
            [head, FLIM, INT, SM]= fastSM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic);
            
            %             [head, FLIM, INT, SM] = SM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, flagauto, pic);
        end
    else
        [head, FLIM, INT] = ROI_FLIM(name, IRF, pic);
    end
elseif strcmpi(flag,'pattern match')|| strcmpi(flag,'lifetime')|| strcmpi(flag,'random')|| strcmpi(flag,'combined')
    lamem       = [];
    lamex       = dye.lamex;    % ?m
    qy          = [];
    tau_free    = [];
    
    n0       = layers.n0;
    n1       = layers.n1;
    n        = layers.n;
    d0       = layers.d0;       % ?m
    d1       = layers.d1;       % ?m
    d        = layers.d;        % ?m
    
    NA       = mic.NA;
    focpos   = mic.focpos;
    if isfield(mic,'pattern')
        pattern = mic.pattern;
    else
        pattern = [];
    end
    %     [head, FLIM, INT, SM] = SM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, flagauto, pic);
    [head, FLIM, INT, SM]= fastSM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, 1, pic);
elseif strcmpi(flag,'intensity')
    lamex=[];
    n0=[];
    n=[];
    n1=[];
    d0=[];
    d=[];
    d1=[];
    NA=[];
    focpos=[];
    pattern=[];
    [head, FLIM, INT, SM] = SM_FLIM(name, IRF, pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, 0, pic);
elseif strcmpi(flag,'show MIET curves')
    lamem    = dye.lamem;
    lamex    = dye.lamex;
    qy       = dye.qy;
    tau_free = dye.tau_free;
    curveType= dye.CurveType;MI
    
    n0       = layers.n0;
    n1       = layers.n1;
    n        = layers.n;
    d0       = layers.d0;
    d1       = layers.d1;
    d        = layers.d;
    al_res   = 10;
    if ~ischar(d)
        [z, lifecurve] = MIET_calc(al_res, lamem*1e3, n0, n, n1, d0*1e3, d*1e3, d1*1e3, qy, tau_free, [], curveType); % all lengths in nm
    else
        [z, lifecurve] = MIET_calc(al_res, lamem*1e3, n0, n, n1, d0*1e3, d, d1*1e3, qy, tau_free, [], curveType); % all lengths in nm
    end
elseif strcmpi(flag,'ROI FLIM')
    [head, FLIM, INT] = ROI_FLIM(name, IRF, pic);
end


if strcmpi(flag,'MIET')
    display('Determining axial distances...')
    maxch   = size(INT.tag,3);
    al_res  = unique(SM.al_res);
    if ~ischar(d)
        [z, lifecurve] = MIET_calc(al_res, lamem*1e3, n0, n, n1, d0*1e3, d*1e3, d1*1e3, qy, tau_free, [], curveType); % all lengths in nm
    else
        [z, lifecurve] = MIET_calc(al_res, lamem*1e3, n0, n, n1, d0*1e3, d, d1*1e3, qy, tau_free, [], curveType); % all lengths in nm
    end
    for ch=1:maxch
       
        if iscell(FLIM.life_matav)
            lifetimes = FLIM.life_matav{ch};
        else
            if size(FLIM.life_matav,1)<4
                lifetimes = FLIM.life_matav(ch,:);
            else
                lifetimes = FLIM.life_matav(:);
            end
            %             for num = 1:numel(FLIM.life_mat)
            %
            %                 lifetimes(num) = FLIM.life_mat{num}(end);
            %             end
        end
         if iscell(SM.theta)
            theta = SM.theta{ch};
         else
            if size( SM.theta,1)<4
                theta     = SM.theta(:);
            else
                theta = SM.theta(:,ch);
            end
        end
        if isnan(theta)
            theta = NaN*lifetimes;
        end
        [ax_pos]  = Axial_pos(theta, lifetimes,z, lifecurve);
        if iscell(INT.field)
            field     = INT.field{ch};
        else
            field = INT.field;
        end
        posimm=zeros(size(field));
        for i=1:size(field,3)
            posimm(:,:,i,ch) = field(:,:,i).*ax_pos(i);
        end
        posimm_c=zeros(size(posimm,1),size(posimm,2),maxch);
        for x=1:size(posimm,1)
            for y=1:size(posimm,2)
                posimm_c(x,y,ch)=pmean(posimm(x,y,:,ch));
            end
        end
        
    end
    if pic==1
        figure
        set(gcf,'name','Height Image','NumberTitle','off')
        for ch=1:maxch
            x_cord=head.ImgHdr.X0;
            y_cord=head.ImgHdr.Y0;
            subplot(ch,2,maxch)
            imagesc(x_cord, y_cord, posimm_c(:,:,ch).')
            colormap('copper');
            colorbar
        end
    end
    save([name(1:end-4) '_miet_analysis.mat'],'posimm','posimm_c','dye','mic','layers','z','lifecurve');
end

if strcmpi(flag,'ROI MIET')
    display('Determining axial distances...')
    maxch   = size(INT.tag,3);
    [z, lifecurve] = MIET_calc(NaN, lamem, n0, n, n1, d0, d, d1, qy, tau_free, [], curveType);
    limit_LT_ind = find(isnan(lifecurve),1); % index of 1st invalid height-lifetime-pair
    if isempty(limit_LT_ind) % if we didn't simulate up to oscillations:
        limit_LT = lifecurve(end);             % use the whole curve
    else % if we know the first invalid height-lifetime-pair:
        limit_LT = lifecurve(limit_LT_ind-1);  % use curve up to last valid height-lifetime-pair
    end
    for ch=1:maxch
        lifetimes = FLIM.life_matav(:,ch);  % lifetimes of different ROIs
        coord     = INT.coord(:,ch);        % (x,y) coordinates of ROIs
        ax_pos=zeros(numel(lifetimes),1);   % will store heights of ROIs
        heightImage = zeros(size(INT.tag)); % will store final image of whole field of view
        for i=1:numel(lifetimes)
            % Test if lifetime is above unique range or smaller
            % than smallest value of calibration curve
            if lifetimes(i)>limit_LT || lifetimes(i)<lifecurve(1)
                ax_pos(i)=NaN;
            else
                % Find 1st element of curve that is larger than current matrix
                % element, interpolate between it and the previous curve element.
                for y=2:length(lifecurve)
                    if lifecurve(y) >= lifetimes(i)
                        ax_pos(i)=z(y-1)+(lifetimes(i)-lifecurve(y-1))*...
                            (z(y)-z(y-1))/(lifecurve(y)-lifecurve(y-1));
                        break;
                    end
                end
            end
            for j=1:size(coord{i},1)
                heightImage(coord{i}(j,2),coord{i}(j,1)) = ax_pos(i);
            end
        end
    end
    if pic==1
        figure
        set(gcf,'name','Height Image','NumberTitle','off')
        for ch=1:maxch
            x_cord=head.ImgHdr.Xcord;
            y_cord=head.ImgHdr.Ycord;
            subplot(ch,1,maxch)
            imagesc(x_cord, y_cord, heightImage)
            colormap('hot');
            colorbar
        end
    end
    save([name(1:end-4) '_ROImiet_analysis.mat'],'heightImage','dye','mic','layers','z','lifecurve');
end

end
