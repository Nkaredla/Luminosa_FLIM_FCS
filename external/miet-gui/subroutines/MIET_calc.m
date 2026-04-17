function [z, life] = MIET_calc(al_res, lamem, n0, n, n1, d0, d, d1, qyield, tau_free, fig, curveType)
% MIET_Calc calculates the calibration curves for different angles given
% the single molecule's quantum yield, free lifetime and the wavelength.
% The program needs the layer information as inputs. All lengths in nm.

% (c) Narain Karedla, 2018

if nargin<10 || isempty(fig)
    fig = 0;
else
    fig=1;
end

if qyield>1
    qyield = qyield/100;
end

if ~isempty(curveType) && strcmp(curveType,'maximum')
    maxCurve = true;    % calculate MIET calibration curve to 1st maximum
else
    maxCurve = false;   % calculate MIET calibration curve to highest unambiguous value (= 1st minimum)
end

if isstruct(lamem) % use emission spectrum of dye?
    polychrome = true;
    metals=load('metals.mat');        % refr. indices depend on wavelength
    spectrum = load(lamem.SpectrumFile);
    if spectrum(1,1)<2  % wavelength in micrometers?
        spectrum(:,1)=spectrum(:,1)*1000; % convert to nm
    end
    wavel_small = lamem.Wavel_Small;        % lower boundary of bandpass filter [nm]
    wavel_large = lamem.Wavel_Large;        % upper boundary of bandpass filter [nm}
    spectrum=spectrum(spectrum(:,1)>=wavel_small,:);    % crop to wavelengths...
    spectrum=spectrum(spectrum(:,1)<=wavel_large,:);    % ... fitting the filters
    spectrum(:,1)=round(spectrum(:,1)/5);   % if you multiply this by 5, you get wavelengths grouped in 5nm steps
    spectrIntensity=accumarray(spectrum(:,1)-spectrum(1,1)+1, spectrum(:,2)); % sum intensities of grouped wavelengths
    spectrum=[unique(spectrum(:,1)*5) spectrIntensity/sum(spectrIntensity)];  % normalise spectrum to get probabilities
else
    polychrome = false;
    spectrum=[lamem 1]; % only use one wavelength, given as lamem [nm]
end

d_var=100;              % it doesn't make sense to calculate MIET curve higher than 500nm
if ~ischar(d) % d is numeric
    if numel(n1)>1 || n1~=n     % is there a stack of layers above the molecule's layer?
        if d<100
            z = 1:(d-1)/100:d-1;
        else
            z = 1:1:min(d,d_var)-1; % then calculate MIET calibration curve only inside molecule's layer
        end
        
    else
        z = 1:1:d_var-1;    % else, calculate MIET calibration curve up to a height of 500nm
        d = d_var;          % if molecule's layer and layer above it are identical, user may input arbitrary values for d
        % -> if d<d_var, LifetimeLSimpsExp will produce NaN at z=d, so set d to d=d_var
    end
else % d is a character, indication for variable silica spacer thickness, molecule spin coated directly on top of the spacer
    z = 0.5;
end


n0_backup=n0; n1_backup=n1;  % refractive indices are overwritten for different wavelengths in polychromatic mode
if isnan(al_res)
    life_avg=zeros(numel(z),1);
elseif numel(al_res)==1
    life_avg=zeros(numel(z),numel((90:-al_res:0)));
else
    life_avg = zeros(numel(z),numel(al_res)); % al_res contains now a vector of the theta angles
end

for wavelCount=1:size(spectrum,1) % loop over different wavelengths
    lamem = spectrum(wavelCount,1);
    fac = 2*pi/lamem;
    if polychrome
        for n_index=1:numel(n0) % replace placeholders by correct refr. indices
            if(n0_backup(n_index)>=10)
                n0(n_index)=refrIndex(n0_backup(n_index));
            end
        end
        for n_index=1:numel(n1) % replace placeholders by correct refr. indices
            if(n1_backup(n_index)>=10)
                n1(n_index)=refrIndex(n1_backup(n_index));
            end
        end
    end
    %         [~,~,~,~,qvd,qvu,qpd,qpu,~,~] = LifetimeLSimpsExp(fac*z,n0,n,n1,fac*d0,fac*d,fac*d1);
    if ~ischar(d)
        [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(fac*z,n0,n,n1,fac*d0,fac*d,fac*d1);
        %     [~,~,~,~,qvd,qvu,qpd,qpu,~,~] = LifetimeL(fac*z,n0,n,n1,fac*d0,fac*d,fac*d1);
    else
        tmpd = 0.5;
         z = 1:d_var-1;
        for l = 1:numel(z)
            
            [lvd(l),lvu(l),lpd(l),lpu(l),qvd(l),qvu(l),qpd(l),qpu(l)] = LifetimeL(fac*tmpd,n0,n,n1,fac*([d0 z(l)]),fac*0,fac*d1);
        end
       
    end
    
    Sv = qvu+qvd;
    Sp = qpu+qpd;
    
    if isnan(al_res) % use average oritentation (if you don't have angle information)
        life=tau_free./(1-qyield+qyield*(Sv.'+2.*Sp.')/(4*n));
        peak = find(diff(life)<0,1); % starting point of oscillations: 1st point after 1st maximum
        if ~isempty(peak) % do we even simulate up to the oscillations?
            if maxCurve
                life(peak:end) = NaN; % keep values up to starting point of oscillations
            else
                limit_LT = min(life(peak:end)); % lifetime up to which values are unique
                limit_height = find(life>limit_LT,1); % 1st non-unique lifetime value
                life(limit_height:end) = NaN;
            end
        end
    else            % calculate one curve for each angle theta
        if numel(al_res)==1
            theta = (90:-al_res:0).*pi./180;
        else
            theta = al_res;
        end
        life=zeros(numel(z),numel(theta));
        for i=1:numel(theta)
            Sr=Sv.*(cos(theta(i)).^2)+Sp.*(sin(theta(i)).^2);
            life(:,i)=tau_free./((1-qyield) + Sr/(4/3*n)*qyield);
        end
        for i=1:numel(theta)
            peak = find(diff(life(:,i))<0,1); % starting point of oscillations: 1st point after 1st maximum
            if ~isempty(peak) % do we even simulate up to the oscillations?
                if maxCurve
                    life(peak:end,i)=NaN; % keep values up to starting point of oscillations
                else
                    limit_LT = min(life(peak:end,i)); % lifetime up to which values are unique
                    [limit_height,~] = find(life>limit_LT,1); % 1st non-unique lifetime value
                    life(limit_height:end,i) = NaN;
                end
            end
        end
    end
    life_avg = life_avg + spectrum(wavelCount,2)*life; % weight lieftime-result with relative no. of photons of this wavelength
end

if fig && any(isnan(al_res))
    figure
    set (gcf,'name','MIET Calibration Curve for Random Orientation','NumberTitle','off')
    plot(z,life,'LineWidth',2);
    xlabel('distance from surface (nm)')
    ylabel('lifetime (ns)')
elseif fig
    figure
    set (gcf,'name','MIET Calibration Curves for Various Orientations','NumberTitle','off')
    c = hsv(numel(theta));
    for i=1:numel(theta)
        a(i) = plot(z,life(:,i),'color',c(i,:),'LineWidth',2);
        l{i} = ['polar angle = ',num2str(theta(i).*180./pi),'^{\circ}'];
        hold on
    end
    xlabel('distance from surface (nm)')
    ylabel('lifetime (ns)')
    axis tight
    legend(a,l,'Location','SouthEast')
    hold off
end


% nested function: find refractive index of "material" at wavelength lamem
    function value = refrIndex(material)
        switch material
            case 10
                value=metals.silver(metals.wavelength==lamem);
            case 20
                value=metals.gold(metals.wavelength==lamem);
            case 30
                value=metals.platinum(metals.wavelength==lamem);
            case 40
                value=metals.palladium(metals.wavelength==lamem);
            case 50
                value=metals.copper(metals.wavelength==lamem);
            case 60
                value=metals.aluminum(metals.wavelength==lamem);
            case 70
                value=metals.chromium(metals.wavelength==lamem);
            case 80
                value=metals.titan(metals.wavelength==lamem);
            case 90
                value=metals.tungsten(metals.wavelength==lamem);
            case 100
                value=metals.nickel(metals.wavelength==lamem);
            case 110
                value=metals.beryllium(metals.wavelength==lamem);
            case 120
                value=metals.ito(metals.wavelength==lamem);
            otherwise
                value=[];
                errormsg('The chosen material does not exist.');
        end
    end

end