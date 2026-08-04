function [c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Fluofit(irf, y, p, dt, tau, lim, init, fitMode, plotFlag)
% The function FLUOFIT performs a fit of a multi-exponential decay curve.
% It is called by: 
% [c, offset, A, tau, dc, doffset, dtau, irs, z, t, chi] = fluofit(irf, y, p, dt, tau, limits, init).
% The function arguments are:
% irf 	= 	Instrumental Response Function
% y 	= 	Fluorescence decay data
% p 	= 	Time between laser exciation pulses (in nanoseconds)
% dt 	= 	Time width of one TCSPC channel (in nanoseconds)
% tau 	= 	Initial guess times. Use [] with init > 0 to let
%           DistFluofit choose the component seeds.
% lim   = 	limits for the lifetimes guess times
% init	=	Whether to use a initial guess routine or not
% fitMode = 'mle', 'ls', or 'pirls' (optional, default 'mle')
% plotFlag = draw diagnostic figure (optional, default true)
%
% The return parameters are:
% c	=	Color Shift (time shift of the IRF with respect to the fluorescence curve)
% offset	=	Offset
% A	    =   Amplitudes of the different decay components
% tau	=	Decay times of the different decay components
% dc	=	Color shift error
% doffset	= 	Offset error
% dtau	=	Decay times error
% irs	=	IRF, shifted by the value of the colorshift
% zz	    Fitted fluorecence component curves
% t     =   time axis
% chi   =   chi2 value
% 
% The program needs the following m-files: simplex.m, lsfit.m, mlfit.m, and convol.m.
% (c) 1996 J�rg Enderlein

if nargin < 8 || isempty(fitMode)
    fitMode = 'mle';
end
fitMode = lower(char(fitMode));
fitMode = strrep(fitMode, '-', '_');

switch fitMode
    case {'mle', 'ml', 'maximum_likelihood'}
        fitfun = 'mlfit';
        ampMode = 'mle';
    case {'ls', 'least_squares'}
        fitfun = 'lsfit';
        ampMode = 'ls';
    case {'pirls', 'pirlsnonneg', 'poisson_irls'}
        fitfun = 'pirlsfit';
        ampMode = 'pirls';
    otherwise
        error('fitMode must be ''mle'', ''ls'', or ''pirls''.');
end

if nargin < 9 || isempty(plotFlag)
    plotFlag = true;
end

if plotFlag
    close all
end
irf = irf(:);
offset = 0;
y = y(:);
n = length(irf); 
if nargin > 4
    tauInput = tau;
else
    tauInput = [];
end
if nargin>6
    if isempty(init)
        init = 1;
    end
elseif nargin>4 
    init = 0;
else
    init = 1;
end
if isempty(tauInput) && init == 0
    init = 1;
end

if init>0 
    try
        [~, tauRates, ~, c] = DistFluofit(irf, y, p, dt, [], 1);
        tau = distRatesToLifetimeSeeds(tauRates);
    catch
        [distAmps, tauRates, ~, c] = DistFluofit(irf, y, p, dt, [], 0);
        tau = distDistributionToLifetimeSeeds(distAmps, tauRates);
    end
    if isempty(tau)
        tau = validLifetimeSeeds(tauInput);
    end
    if isempty(tau)
        tau = fallbackLifetimeSeed(y, p, dt);
    end
else
    c = 0;
end

if (nargin<6)||isempty(lim)
    lim = [zeros(1,length(tau)) 100.*ones(1,length(tau))];
else
    lim = double(lim(:)).';
    if numel(lim) ~= 2*numel(tau)
        error('lim must contain lower and upper lifetime bounds for each fitted component.');
    end
end;

p = p/dt;
tp = (1:p)';
tau = tau(:)'/dt; 
lim_min = lim(1:numel(tau))./dt;
lim_max = lim(numel(tau)+1:end)./dt;
t = 1:length(y);
m = length(tau);
x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = Convol(irs, x);
z = [ones(size(z,1),1) z];
basis = z;
[A, z] = solveFluofitAmplitudes(basis, y, ampMode);
zz = basis.*(ones(size(basis,1),1)*A');
if plotFlag
    close all
end
if init<2
%     disp('Fit =                Parameters =');
    param = [c; tau'];
    % Decay times and Offset are assumed to be positive.
    paramin = [-1/dt lim_min];
    paramax = [ 1/dt lim_max];
    [param, dparam] = Simplex(fitfun, param, paramin, paramax, 1e-12, Inf, irf(:), y(:), p);
    c = param(1);
    dc = dparam(1);
    tau = param(2:length(param))';
    dtau = dparam(2:length(param));
    x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
    irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
    z = Convol(irs, x);
    z = [ones(size(z,1),1) z];
    zsum = sum(z);
    zsum(~isfinite(zsum) | zsum == 0) = 1;
    z = z./(ones(n,1)*zsum);
%     A = z\y;
    [A, zfit] = solveFluofitAmplitudes(z, y, ampMode);
    zz = z.*(ones(size(z,1),1)*A');
    z = zfit;
    dtau = dtau;
    dc = dt*dc;
else
    dtau = 0;
    dc = 0;
end
chi = sum((y-z).^2./max(abs(z), eps))/(n-m);
t = dt*t;
tau = dt*tau';
c = dt*c;
offset = zz(1,1); 
A(1) = [];
if plotFlag
	hold off
    subplot('position',[0.1 0.4 0.8 0.5])
	plot(t,log10(y),t,log10(irs),t,log10(z));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);
	axis(v);
	xlabel('Time in ns');
	ylabel('Log Count');
	s = sprintf('COF = %3.3f   %3.3f', c, offset);
	text(max(t)/2,v(4)-0.05*(v(4)-v(3)),s);
	s = ['AMP = '];
	for i=1:length(A)
		s = [s sprintf('%1.3f',A(i)/sum(A)) '   '];
	end
	text(max(t)/2,v(4)-0.12*(v(4)-v(3)),s);
	s = ['TAU = '];
	for i=1:length(tau)
		s = [s sprintf('%3.3f',tau(i)) '   '];
	end
	text(max(t)/2,v(4)-0.19*(v(4)-v(3)),s);
    subplot('position',[0.1 0.1 0.8 0.2])
	plot(t,(y-z)./sqrt(abs(z)));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);

    axis(v);
	xlabel('Time in ns');
	ylabel('Residue');
	s = sprintf('%3.3f', chi);
	text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
    set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])
end
end

function [A, zfit] = solveFluofitAmplitudes(M, y, ampMode)
M = real(double(M));
y = real(double(y(:)));

switch lower(ampMode)
    case 'ls'
        A = M \ y;
    case 'pirls'
        if exist('PIRLSnonneg', 'file') == 2
            A = PIRLSnonneg(M, y, 10);
        else
            warning('PIRLSnonneg.m not found; falling back to lsqnonneg.');
            A = lsqnonneg(M, y);
        end
    otherwise
        A = lsqnonneg(M, y);
end

A = real(double(A));
A(~isfinite(A)) = 0;
if ~strcmpi(ampMode, 'ls')
    A = max(A, 0);
end
zfit = M * A;
end

function tauNs = distRatesToLifetimeSeeds(tauRates)
tauRates = real(double(tauRates(:).'));
tauRates = tauRates(isfinite(tauRates) & tauRates > 0);
tauNs = 1 ./ tauRates;
tauNs = tauNs(isfinite(tauNs) & tauNs > 0);
tauNs = unique(tauNs, 'stable');
end

function tauNs = distDistributionToLifetimeSeeds(distAmps, tauRates)
distAmps = real(double(distAmps(:).'));
tauRates = real(double(tauRates(:).'));
valid = isfinite(distAmps) & distAmps > 0 & isfinite(tauRates) & tauRates > 0;
if ~any(valid)
    tauNs = [];
    return;
end

distAmps(~valid) = 0;
threshold = 0.2 * max(distAmps(valid));
active = valid & distAmps >= threshold;
if ~any(active)
    [~, idx] = max(distAmps);
    active(idx) = true;
end

edges = diff([false active false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
tauNs = zeros(1, numel(starts));
for jj = 1:numel(starts)
    idx = starts(jj):stops(jj);
    weights = distAmps(idx);
    totalWeight = sum(weights);
    if totalWeight <= 0
        continue;
    end
    lifetimeNs = 1 ./ tauRates(idx);
    tauNs(jj) = sum(weights .* lifetimeNs) / totalWeight;
end
tauNs = tauNs(isfinite(tauNs) & tauNs > 0);
tauNs = unique(tauNs, 'stable');
end

function tauNs = validLifetimeSeeds(tauInput)
tauNs = real(double(tauInput(:).'));
tauNs = tauNs(isfinite(tauNs) & tauNs > 0);
tauNs = unique(tauNs, 'stable');
end

function tauNs = fallbackLifetimeSeed(y, p, dt)
y = real(double(y(:)));
y(~isfinite(y) | y < 0) = 0;
t = (0:numel(y)-1).' * dt;
s = sum(y);
if s > 0
    tauNs = sum(t .* y) / s;
else
    tauNs = p / 4;
end
if ~isfinite(tauNs) || tauNs <= 0
    tauNs = p / 4;
end
tauNs = min(max(tauNs, dt), p);
end
