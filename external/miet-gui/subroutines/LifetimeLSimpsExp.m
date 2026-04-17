function [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeLSimpsExp(z,n0,n,n1,d0,d,d1)

% [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeLSimpsExp(z,n0,n,n1,d0,d,d1) 
% calculates the lifetime and quenching of an oscillating dipole at a 
% distance z from an interface within a layer using Simpson's rule with
% exponential integration points for the integration of evanescent waves.
% Normalization: free dipole emission rate = 4/3*n
%
% INPUT VARIABLES:
% z  - molecule's distance from the bottom of its layer -> row vector
% n0 - vector of refractive indices of the stack below the molecule's layer
% n  - refractive index of the molecule's layer, must be a real number
% n1 - vector of refractive indices of the stack above the molecule's layer
% d0 - vector of layer thickness values of the stack below the molecule's layer ( length(d0)=length(n0)-1 )
% d  - thickness of molecule's layer
% d1 - vector of layer thickness values of the stack above the molecule's layer ( length(d1)=length(n1)-1 )
%
% OUTPUT VARIABLES:
% lvd - emission power of vertical dipole into infinite lower halfspace
% lvu - emission power of vertical dipole into infinite upper halfspace
% lpd - emission power of parallel dipole into infinite lower halfspace
% lpu - emission power of parallel dipole into infinite upper halfspace
% qvd - total emission power of vert. dipole towards lower stack of layers
% qvu - total emission power of vert. dipole towards upper stack of layers
% qpd - total emission power of par. dipole towards lower stack of layers
% qpu - total emission power of par. dipole towards upper stack of layers
% q and l values differ if there is absorption in a layer
%
% (c) Joerg Enderlein and Daja Ruhlandt, 2014


%% PART 1: total emission of modes that are real in the molecule's layer

z = z(:)';              % row vector -> column vector
col = ones(size(z));    % column vector, same length as z, filled with ones
dw = 1e-3;              % number of integration intervals = 1/dw
% integrate using midpoint rule: col. vec. w (has to be col.vec. for Fresnel()),
% contains the midpoints of the integration intervals
if all(imag([n0 n1])==0) && max([n0 n1])<n 
    % if there are no metal layers and if no evanescent modes from the 
    % molecule's layer can propagate anywhere else: 
    % lower integration boundary chosen such that there is one layer (the 
    % one with the highest index of refraction apart from the molecule's layer) 
    % such that q=n_{that layer}, in all the other layers w will be imaginary
    % => evanescent wave, but at least one layer where it can propagate 
    % as a normal wave
    % upper integration boundary = n (=> q=0), in 1/dw steps
    w = sqrt(n^2-max([n0 n1])^2) + (0.5:1/dw)'*dw*(n-sqrt(n^2-max([n0 n1])^2));
else
    % intervals from 0 to n in 1/dw steps
    w = (0.5:1/dw)'*dw*n; 
end
dw = w(2)-w(1);     % step width for integration (not equal to number in line 16!)

% get effective refl. & transm. coefficients for interaction with upper and
% lower stack of media for different directions of propagation (via "w")
[rpu, rsu, tpu, tsu] = Fresnel(w,[n n1],d1);                     % upper
[rpd, rsd, tpd, tsd] = Fresnel(w,[n n0(end:-1:1)],d0(end:-1:1)); % lower
rpu = rpu.'; rpd = rpd.'; rsu = rsu.'; rsd = rsd.'; % row vec. -> col. vec.
tpu = tpu.'; tpd = tpd.'; tsu = tsu.'; tsd = tsd.'; % row vec. -> col. vec.

ez0 = exp(2*1i*w*z);    % phase shift betw. 2 upwards trav. waves where one was initially trav. downwards & the other upwards
ez = exp(2*1i*w*(d-z)); % phase shift betw. 2 downwards trav. waves where one was initially trav. downwards & the other upwards
ed = exp(2*1i*w*d);     % phase shift from traversing molecule's layer twice

% take multiple reflections within molecule's layer into account;
% we calculate the amplitudes on the molecule's side of the interfaces
tp = (1./(1-rpu.*rpd.*ed));     % p-wave
ts = (1./(1-rsu.*rsd.*ed));     % s-wave

% calculate amplitude of el. field at the bottom of the molecule's layer
% for vertical (v) or parallel (pc,p-wave; ps,s-wave) oriented dipole; 
% (1) is wave travelling towards interface and (2) is wave that has just
% been reflected at interface
v1 = tp*col + ((rpu.*tp)*col).*ez;              
v2 = (rpd*col).*v1;  
pc1 = tp*col - ((rpu.*tp)*col).*ez;             
pc2 = (rpd*col).*pc1; 
ps1 = ts*col + ((rsu.*ts)*col).*ez;             
ps2 = (rsd*col).*ps1; 

% energy flux towards bottom = flux from downwards wave - flux from upwards
% wave; in the formulas (skript!) you find factors sin^2(theta) for vert. d., 
% 0.5 for parallel dipole and w^2/n^2 for p-wave of parallel dipole
st = sqrt(1-(w/n).^2);  % = sin(theta)
qvd = dw*sum((st.^2*col).*(abs(v1).^2-abs(v2).^2));
qpd = dw*sum((w.^2/n^2*col).*(abs(pc1).^2-abs(pc2).^2) + abs(ps1).^2-abs(ps2).^2)/2;

% calculate amplitude of el. field at the top of the molecule's layer
% for vertical (v) or parallel (pc,p-wave; ps,s-wave) oriented dipole; 
% (1) is wave travelling towards interface and (2) is wave that has just
% been reflected at interface
v1 = tp*col + ((rpd.*tp)*col).*ez0;
v2 = (rpu*col).*v1;
pc1 = tp*col - ((rpd.*tp)*col).*ez0;
pc2 = (rpu*col).*pc1;
ps1 = ts*col + ((rsd.*ts)*col).*ez0;
ps2 = (rsu*col).*ps1;

% energy flux towards top = flux from upwards wave - flux from downwards
% wave; in the formulas (skript!) you find factors sin^2(theta) for vert. d., 
% 0.5 for parallel dipole and w^2/n^2 for p-wave of parallel dipole
qvu = dw*sum((st.^2*col).*(abs(v1).^2-abs(v2).^2));
qpu = dw*sum((w.^2/n^2*col).*(abs(pc1).^2-abs(pc2).^2) + abs(ps1).^2-abs(ps2).^2)/2;


%% PART 2: emission of modes that are real in the molecule's layer to infinite halfspaces

if imag(n0(1))==0   % if lower halfspace is metal: emission damped to 0
    tp = (tpd./(1-rpu.*rpd.*ed));
    ts = (tsd./(1-rsu.*rsd.*ed));
    tp1 = tp.*(rpu.*ed);
    ts1 = ts.*(rsu.*ed);
    ez = exp(1i*w*z);
    % ww=0 if w_{bottom layer} is imaginary (emission would be damped to 0)
    % ww=w_{bottom layer}/w_{molecule's layer}
    ww = real(sqrt((n0(1)^2-n^2)./w.^2+1))*col;
    % calculate amplitudes at upper edge of bottommost medium
    v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
    pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
    ps = (ts*col).*ez + (ts1*col)./ez;
    % sum over those w from part 1 that are real in the lower halfspace
    lvd = dw*sum(ww.*abs(v).^2);
    lpd = dw*sum(ww.*(abs(pc).^2 + abs(ps).^2))/2;
end

if imag(n1(end))==0 % if upper halfspace is metal: emission damped to 0
    tp = (tpu./(1-rpu.*rpd.*ed));
    ts = (tsu./(1-rsu.*rsd.*ed));
    tp1 = tp.*(rpd.*ed);
    ts1 = ts.*(rsd.*ed);
    ez = exp(1i*w*(d-z));
    % ww=0 if w_{bottom layer} is imaginary (emission would be damped to 0)
    % ww=w_{bottom layer}/w_{molecule's layer}    
    ww = real(sqrt((n1(end)^2-n^2)./w.^2+1))*col;
    % calculate amplitudes at lower edge of topmost medium
    v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
    pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
    ps = (ts*col).*ez + (ts1*col)./ez;
    % sum over those w from part 1 that are real in the upper halfspace
    lvu = dw*sum(ww.*abs(v).^2);     % shouldn't it be ww.^2? and where is the n?
    lpu = dw*sum(ww.*(abs(pc).^2 + abs(ps).^2))/2;
end


%% PART 3: waveguide modes 
% (only occur when no metals are present and the indices of refraction of 
%   the top and bottom halfspace are both smaller than the largest index of
%   refraction in any of the other layers)

qv = 0;
qp = 0;

if all(imag([n0 n1])==0) && max([n0(2:end) n n1(1:end-1)])>max([n0(1) n1(end)]) % waveguide modes
    dw = 1i*1e-7;
    [ppos, spos] = HashWaveguideMode([n0 n n1],[d0 d d1]);        
    if ~isempty(ppos)
        wp = sqrt(n^2-ppos.^2);
        for j=1:length(ppos)
            ww = wp(j)+(-0.5:0.5)'*dw;
            rpd = Fresnel(ww, [n, n0(end:-1:1)], d0(end:-1:1));        
            rpu = Fresnel(ww, [n, n1], d1);
            rpu = rpu.'*col; rpd = rpd.'*col; 
            fp = diff((1-rpd.*rpu.*exp(2*1i*ww*d*col))./(rpd.*exp(2*1i*ww*z)+rpu.*exp(2*1i*ww*(d-z))+2*rpu.*rpd.*exp(2*1i*d*ww*col)))/dw;
            qv = qv + 4*pi*ppos(j).^2.*imag(1./fp)/n^2;
            qp = qp + 2*pi*imag(abs(wp(j)).^2./fp)/n^2;
        end
    end
    if ~isempty(spos)
        ws = sqrt(n^2-spos.^2);
        for j=1:length(spos)
            ww = ws(j)+(-0.5:0.5)'*dw;
            [~, rsd] = Fresnel(ww, [n, n0(end:-1:1)], d0(end:-1:1));
            [~, rsu] = Fresnel(ww, [n, n1], d1);
            rsu = rsu.'*col; rsd = rsd.'*col;
            fs = diff((1-rsd.*rsu.*exp(2*1i*ww*d*col))./(rsd.*exp(2*1i*ww*z)+rsu.*exp(2*1i*ww*(d-z))+2*rsu.*rsd.*exp(2*1i*d*ww*col)))/dw;
            qp = qp + 2*pi*imag(1./fs);
        end
    end
end


%% PART 4: modes that are imaginary in the molecule's layer

% CAREFUL: from now on, w denotes the imaginary part of the vert. component of k

% modes that are evanescent in the molecule's layer may be able to propagate/be absorbed
% in both the upper and the lower stack of layers or only in one of them or in none of them:
% treat emission into lower and upper halfspace separately

dIsMeaningful = (numel(n1)>1 || n1~=n); % if n1==n, i.e. if there is no "upper medium", the user may type any value for d

% lower halfspace ---------------------------------------------------------
if all(imag(n0)==0)
    % if n>=max[n0] & no metals are present, wmax=0 and we are finished (see next if-statement)
    % if n<max[n0] & no metals are present, imag(vert. comp. of k) goes
    % from exp(-6) to wmax=sqrt(n_max^2-n^2) (evanescent waves that can 
    % propagate in at least one of the other media)
    wmax = sqrt(max([n0 n]).^2-n^2);
    % w = exp(-6 + (6+log(wmax))*(0.5:1e3)'/1e3);
    % dw = (6+log(wmax))/1e3;
    if wmax>0
        w_temp = exp(-6 + (6+log(wmax))*(0:499)'/499); % edges of intervals exponentially distributed 
        w = zeros(2*numel(w_temp)-1,1);                      % all w-values: ...
        w(1:2:end)=w_temp(:);                           % ... edges of intervals
        w(2:2:end-1)=(w_temp(2:end)+w_temp(1:end-1))/2; % ... middle points of intervals
        w_diff=(w(3:2:end)-w(1:2:end-2))*ones(1,numel(z)); % matrix of interval widths
    end
else
    % if metals are present, imag(vert. component of k) goes from exp(-6) to a 
    % very high value; the rest of the integral to infinity is done analytically
    z_temp=z; z_temp(z>=d)=[];  % if there is a z==d, wmax will be inf; if there are z>d, wmax will be <0
    if dIsMeaningful 
        wmax = abs(log(1e-3))/min([d d0 d1 d-z_temp]);
    else % if n1==n, i.e. if there is no "upper medium", the user may type any value for d -> d and d-z are meaningless
        wmax = abs(log(1e-3))/min([d0 d1]);
    end
    w_temp = exp(-6 + (6+log(wmax))*(0:499)'/499); % edges of intervals exponentially distributed
    w = zeros(2*numel(w_temp)-1,1);                 % all w-values: ... 
    w(1:2:end)=w_temp(:);                           % ... edges of intervals
    w(2:2:end-1)=(w_temp(2:end)+w_temp(1:end-1))/2; % ... middle points of intervals
    w_diff=(w(3:2:end)-w(1:2:end-2))*ones(1,numel(z)); % matrix of interval widths
end
if wmax>0
    st = sqrt(1+(w/n).^2);  % q/n
    [rpu, rsu, ~, ~]     = Fresnel(1i*w,[n n1],d1);
    [rpd, rsd, tpd, tsd] = Fresnel(1i*w,[n n0(end:-1:1)],d0(end:-1:1));
    
    rpu = rpu.'; rpd = rpd.'; rsu = rsu.'; rsd = rsd.';
                 tpd = tpd.';              tsd = tsd.';
    ez0 = exp(-2*w*z);    % phase shift betw. 2 upwards trav. waves where one was initially trav. downwards & the other upwards
    if dIsMeaningful
        ez = exp(-2*w*(d-z)); % phase shift betw. 2 downwards trav. waves where one was initially trav. downwards & the other upwards
        ed = exp(-2*w*d);     % phase shift from traversing molecule's layer twice
    else
        ez = exp(-w*col);     % will be multiplied by zero anyway, so value doesn't matter
        ed = exp(-w);         %                            "
    end

    tp = (1./(1-rpu.*rpd.*ed));
    ts = (1./(1-rsu.*rsd.*ed));

    v1 = tp*col + ((rpu.*tp)*col).*ez;
    v2 = (rpd*col).*v1;
    pc1 = tp*col - ((rpu.*tp)*col).*ez;
    pc2 = (rpd*col).*pc1;
    ps1 = ts*col + ((rsu.*ts)*col).*ez;
    ps2 = (rsd*col).*ps1;

    summand = 2*(st.^2*col).*ez0.*imag(conj(v1).*v2);
    qvd = qvd + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );

    summand = ez0.*( (w.^2/n^2*col).*imag(conj(pc1).*pc2) + imag(conj(ps1).*ps2) );
    qpd = qpd + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
    
    if any(imag(n0)~=0) % ALL evanescent waves have to be taken into account: analytical result of integral tail
        % rpd(imag(w)=infty) = (n0(end)/n - n/n0(end)) / (n0(end)/n + n/n0(end))
        % rsd(imag(w)=infty) = 0
        qvd = qvd + 1./(z*n^2).* (n^2 + wmax^2 + wmax./z + 1/2./z.^2) .*exp(-2*wmax*z)...
            .*imag( abs(v1(end,:)).^2  *(n0(end)/n-n/n0(end))/(n0(end)/n+n/n0(end)) );
        qpd = qpd + 1./(2*z*n^2).* (wmax^2 + wmax./z + 1/2./z.^2) .*exp(-2*wmax*z)...
            .*imag( abs(pc1(end,:)).^2 *(n0(end)/n-n/n0(end))/(n0(end)/n+n/n0(end)) );
    end
    
    if imag(n0(1))==0 && n0(1)>n % emission into far-field only if modes are real in far-field
        if dIsMeaningful
            ez = exp(-w*z);       % phase shift from traveling from molecule to top of bottom stack
            ez0= exp(-w*(2*d-z)); % phase shift from traveling from molecule to upper stack and back to bottom stack
        else
            ez = exp(-w*z);       % phase shift from traveling from molecule to top of bottom stack
            ez0= exp(-w*z);       % will be multiplied by zero anyway, so value doesn't matter
        end
        tp = (tpd./(1-rpu.*rpd.*ed));
        ts = (tsd./(1-rsu.*rsd.*ed));
        tp1 = tp.*rpu;
        ts1 = ts.*rsu;

        v = ((st.*tp)*col).*ez + ((st.*tp1)*col).*ez0;
        pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col).*ez0;
        ps = (ts*col).*ez + (ts1*col).*ez0;
        ww = real(sqrt((n0(1)^2-n^2)./w.^2-1))*col;

        summand = ww.*abs(v).^2; % only sum over modes that are real in far-field (otherwise ww=0)
        lvd = lvd + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
        summand = ww.*(abs(pc).^2 + abs(ps).^2);
        lpd = lpd + 1/12*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
    end
end


% upper halfspace ---------------------------------------------------------
if all(imag(n1)==0)
    % if n>=max[n1] & no metals are present, wmax=0 and we are finished (see next if-statement)
    % if n<max[n1] & no metals are present, imag(vert. comp. of k) goes
    % from exp(-6) to wmax=sqrt(n_max^2-n^2) (evanescent waves that can 
    % propagate in at least one of the other media)
    wmax = sqrt(max([n1 n]).^2-n^2);
    % w = exp(-6 + (6+log(wmax))*(0.5:1e3)'/1e3);
    % dw = (6+log(wmax))/1e3;
    if wmax>0
        w_temp = exp(-6 + (6+log(wmax))*(0:499)'/499); % edges of intervals exponentially distributed 
        w = zeros(2*numel(w_temp)-1,1);                      % all w-values: ...
        w(1:2:end)=w_temp(:);                           % ... edges of intervals
        w(2:2:end-1)=(w_temp(2:end)+w_temp(1:end-1))/2; % ... middle points of intervals
        w_diff=(w(3:2:end)-w(1:2:end-2))*ones(1,numel(z)); % matrix of interval widths
    end
else
    % if metals are present, imag(vert. component of k) goes from exp(-6) to a 
    % very high value; the rest of the integral to infinity is done analytically
    z_temp=z; z_temp(z>=d)=[];  % if there is a z==d, wmax will be inf; if there are z>d, wmax will be <0
    if dIsMeaningful
        wmax = abs(log(1e-3))/min([d d0 d1 d-z_temp]);
    else % if n1==n, i.e. if there is no "upper medium", the user may type any value for d -> d and d-z are meaningless
        wmax = abs(log(1e-3))/min([d0 d1]);
    end
    w_temp = exp(-6 + (6+log(wmax))*(0:499)'/499); % edges of intervals exponentially distributed
    w = zeros(2*numel(w_temp)-1,1);                 % all w-values: ... 
    w(1:2:end)=w_temp(:);                           % ... edges of intervals
    w(2:2:end-1)=(w_temp(2:end)+w_temp(1:end-1))/2; % ... middle points of intervals
    w_diff=(w(3:2:end)-w(1:2:end-2))*ones(1,numel(z)); % matrix of interval widths
end
if wmax>0
    st = sqrt(1+(w/n).^2);  % q/n
    [rpu, rsu, tpu, tsu] = Fresnel(1i*w,[n n1],d1);
    [rpd, rsd, ~, ~]     = Fresnel(1i*w,[n n0(end:-1:1)],d0(end:-1:1));
    
    rpu = rpu.'; rpd = rpd.'; rsu = rsu.'; rsd = rsd.';
    tpu = tpu.';              tsu = tsu.';
    ez  = exp(-2*w*z);    % phase shift betw. 2 upwards trav. waves where one was initially trav. downwards & the other upwards
    if dIsMeaningful
        ez0 = exp(-2*w*(d-z));% square of phase shift from travelling from molecule to bottom of upper stack
        ed  = exp(-2*w*d);    % phase shift from traversing molecule's layer twice
    else
        ez0 = exp(-w*col);    % will be multiplied by zero anyway, so value doesn't matter
        ed  = exp(-w);        %                            "
    end
    
    tp = (1./(1-rpu.*rpd.*ed));
    ts = (1./(1-rsu.*rsd.*ed));

    v1 = tp*col + ((rpd.*tp)*col).*ez;
    v2 = (rpu*col).*v1;
    pc1 = tp*col - ((rpd.*tp)*col).*ez;
    pc2 = (rpu*col).*pc1;
    ps1 = ts*col + ((rsd.*ts)*col).*ez;
    ps2 = (rsu*col).*ps1;

    summand = 2*(st.^2*col).*ez0.*imag(conj(v1).*v2);
    qvu = qvu + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );

    summand = ez0.*( (w.^2/n^2*col).*imag(conj(pc1).*pc2) + imag(conj(ps1).*ps2) );
    qpu = qpu + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
    
    if any(imag(n1)~=0) % ALL evanescent waves have to be taken into account: analytical result of integral tail
        % rpu(imag(w)=infty) = (n1(1)/n - n/n1(1)) / (n1(1)/n + n/n1(1))
        % rsu(imag(w)=infty) = 0
        qvu = qvu + 1./((d-z)*n^2).* (n^2 + wmax^2 + wmax./(d-z) + 1/2./(d-z).^2) .*exp(-2*wmax*z)...
            .*imag( abs(v1(end,:)).^2  *(n1(1)/n-n/n1(1))/(n1(1)/n+n/n1(1)) );    
        qpu = qpu + 1./(2*(d-z)*n^2).* (wmax^2 + wmax./(d-z) + 1/2./(d-z).^2) .*exp(-2*wmax*z)...
            .*imag( abs(pc1(end,:)).^2 *(n1(1)/n-n/n1(1))/(n1(1)/n+n/n1(1)) );
    end
    
    if imag(n1(end))==0 && n1(1)>n % emission into far-field only if modes are real in far-field
        if dIsMeaningful
            ez = exp(-w*(d-z));   % phase shift from traveling from molecule to bottom of upper stack
            ez0= exp(-w*(d+z));   % phase shift from traveling from molecule to bottom stack and back to upper stack
        else
            ez = exp(-w*z);       % will be multiplied by zero anyway, so value doesn't matter
            ez0= exp(-w*z);       %                            "
        end
        tp = (tpu./(1-rpu.*rpd.*ed));
        ts = (tsu./(1-rsu.*rsd.*ed));
        tp1 = tp.*rpd;
        ts1 = ts.*rsd;

        v = ((st.*tp)*col).*ez + ((st.*tp1)*col).*ez0;
        pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col).*ez0;
        ps = (ts*col).*ez + (ts1*col).*ez0;
        ww = real(sqrt((n1(end)^2-n^2)./w.^2-1))*col;

        summand = ww.*abs(v).^2; % only sum over modes that are real in far-field (otherwise ww=0)
        lvu = lvu + 1/6*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
        summand = ww.*(abs(pc).^2 + abs(ps).^2);
        lpu = lpu + 1/12*sum( (summand(1:2:end-2,:)+4*summand(2:2:end-1,:)+summand(3:2:end,:)).*w_diff );
    end
end

%% Part 5: no emission into far-field if material is absorbing
if imag(n0(1))~=0
    lvd=0; lpd=0;
end
if imag(n1(end)) ~=0
    lvu=0; lpu=0;
end

end
