function [cx, k, offset, z, t, err] = DistTailfit(y, dt, flag, bild, N, taumax, T)
% The function DistTailfit performs a fit of a distributed decay curve.
% It is called by: 
% [cx, tau, offset, z, t, err] = DistTailfit(y, dt, flag, bild, N, taumax)
% The function arguments are:
% y         = 	fluorescence decay data
% dt        = 	time width of one TCSPC channel (in nanoseconds)
% flag      =   if flag==1 only lifetimes with significant contributions are
%               returned and clusters are averaged before they are returned
% bild      =   if bild==1, two figures are produced: one contains the raw 
%               data and the fitted curve, the other contains the weighted 
%               residuals and the distribution of decay times 
% N         =   number of lifetimes used for the fit, default is 100
% taumax    =   maximum lifetime expected
%
% The return parameters are:
% cx	    =	lifetime distribution (for exponential curves that have been 
%               normalized by division with sum(exp(-tau*t)) with t=(0:numel(y)-1)*dt)
% k         =   used decay rates
% offset    =	constant offset multiplied by numel(y)
% z 	    =	fitted fluorecence curve
% t         =   time axis
% err       =   chi2 value
% 
% (c) 2014 J?rg Enderlein

y = y(:);
n = length(y); 

if nargin<7 || isempty(T)
    longFlag=false;
else
    longFlag=true;
end

if nargin<6 || isempty(taumax)
    taumax=dt*n;
end

if nargin<5 || isempty(N)
    N = 200;
end

if nargin<4 || isempty(bild)
    bild = 0;
end

if nargin<3 || isempty(flag)
    flag = 0;
end

tp = dt*(0:n-1)';
t = (0:n-1)';

k = linspace(1/taumax,1/8/dt,N); % distribution of inverse decay times
if longFlag % take into account that emission might happen in next/third/fourth/... pulse
    M0 = [ones(size(t)) exp(-tp*k)./(1-ones(size(tp))*exp(-T*k))];
else
    M0 = [ones(size(t)) exp(-tp*k)]; % constant function and all exponential functions
end
M0 = M0./(ones(n,1)*sum(M0)); % each curve normalized to area=1

cx = PIRLSnonneg(M0,y); % find amplitudes of normalized exponential functions
z = M0*cx;            % fitted curve
err = sum((z-y).^2./abs(z))/n;

if bild
    figure
    t = dt*t;
    semilogy(t,y,'ob','linewidth',1);
    hold on
    semilogy(t,z,'r','linewidth',2);
    hold off
    
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    axis(v);
    xlabel('time [ns]');
    ylabel('lg count');
    figure;
    subplot(2,1,1);
    plot(t,(y-z)./sqrt(z)); 
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    axis(v);
    xlabel('time [ns]');
    ylabel('weighted residual');
    
    ind=1:length(cx)-2;
    len = length(ind);
    tau = 1./k;
    fac = sqrt(tau(1:end-1)/tau(2:end));
    tau1 = fac*tau(ind);
    subplot(2,1,2)
    semilogx(reshape([fac*tau(ind);fac*tau(ind);tau(ind)/fac;tau(ind)/fac],4*len,1),reshape([0*tau(ind);cx(ind+1)';cx(ind+1)';0*tau(ind)],4*len,1));
    patch(reshape([fac*tau(ind);fac*tau(ind);tau(ind)/fac;tau(ind)/fac],4*len,1),reshape([0*tau(ind);cx(ind+1)';cx(ind+1)';0*tau(ind)],4*len,1),'b');

    xlabel('decay time [ns]');
    ylabel('distribution');
end
k = k';
offset = cx(1).*M0(1);
cx(1) = [];

if flag>0
    cx = cx';
    tmp = cx>0; % find entries with significant contribution
    t = 1:length(tmp);
    t1 = t(tmp(2:end)>tmp(1:end-1)) + 1; % index of "start of cluster"
    t2 = t(tmp(1:end-1)>tmp(2:end));     % index of "end of cluser"
    if isempty(t1)|| isempty(t2) % t1 is empty if only entry is 1st, t2 is empty if only entry is last
        [cx, k, offset, z, t, err] = DistTailfit(y, dt, flag, bild, N, 2*taumax);
    elseif numel(t1)==1 && numel(t2)==1 % only one entry
        k = k(t1);
        cx = cx(t1);
    else
        if t1(1)>t2(1) % happens if smallest lifetime contributes significantly
            t2(1)=[];  % throw out first cluster
        end
        if t1(end)>t2(end) % happens if largest lifetime contributes significantly
            t1(end)=[];    % throw out last cluster
        end
        if length(t1)==length(t2)+1
            t1(end)=[];
        end
        if length(t2)==length(t1)+1
            t2(1)=[];
        end
        tmp_k = zeros(1,numel(t1)); tmp_cx = tmp_k;
        for j=1:length(t1)
            tmp_k(j) = cx(t1(j):t2(j))*k(t1(j):t2(j))/sum(cx(t1(j):t2(j))); % average decay rate of a cluster
            tmp_cx(j)  = sum( cx(t1(j):t2(j)) ); % total weight of a cluster
        end
        k = tmp_k;
        cx = tmp_cx; 
    end
end

