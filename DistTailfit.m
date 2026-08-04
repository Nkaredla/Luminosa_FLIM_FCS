function [cx, k, offset, z, t, err] = DistTailfit(y, dt, flag, bild, N, taumax, T)
%DISTTAILFIT Distributed exponential tail fit without IRF reconvolution.
%
% [cx, k, offset, z, t, err] = DistTailfit(y, dt, flag, bild, N, taumax, T)
%
% y      fluorescence decay tail data
% dt     TCSPC channel width in ns
% flag   if true, return clustered significant decay-rate components
% bild   if true, draw fit/residual/distribution figures
% N      number of distributed rates, default 200
% taumax maximum lifetime in ns, default dt*numel(y)
% T      optional pulse period in ns for wrap-around correction

y = real(double(y(:)));
y(~isfinite(y)) = 0;
y = max(y, 0);
n = length(y);
if n < 3
    error('DistTailfit requires at least three tail samples.');
end

if nargin < 7 || isempty(T)
    longFlag = false;
else
    T = double(T);
    longFlag = isfinite(T) && T > 0;
end

if nargin < 6 || isempty(taumax)
    taumax = dt*n;
end
taumax = double(taumax);
if ~isfinite(taumax) || taumax <= 0
    taumax = dt*n;
end

if nargin < 5 || isempty(N)
    N = 200;
end
N = max(3, round(double(N)));

if nargin < 4 || isempty(bild)
    bild = 0;
end
if nargin < 3 || isempty(flag)
    flag = 0;
end

dt = double(dt);
if ~isfinite(dt) || dt <= 0
    error('DistTailfit requires a positive dt in ns.');
end

tp = dt*(0:n-1)';
t = (0:n-1)';

k = linspace(1/taumax, 1/max(8*dt, eps), N);
if longFlag
    denom = 1 - ones(size(tp))*exp(-T*k);
    denom = max(denom, eps);
    M0 = [ones(size(t)) exp(-tp*k)./denom];
else
    M0 = [ones(size(t)) exp(-tp*k)];
end
M0 = M0./max(ones(n,1)*sum(M0), eps);

cx = PIRLSnonneg(M0, y);
z = M0*cx;
err = sum((z-y).^2./max(abs(z), eps))/n;

if bild
    tNs = dt*t;
    figure;
    semilogy(tNs, y, 'ob', 'linewidth', 1);
    hold on
    semilogy(tNs, z, 'r', 'linewidth', 2);
    hold off
    axis tight
    xlabel('time [ns]');
    ylabel('count');

    figure;
    subplot(2,1,1);
    plot(tNs, (y-z)./sqrt(max(abs(z), eps)));
    axis tight
    xlabel('time [ns]');
    ylabel('weighted residual');

    ind = 1:length(cx)-2;
    if ~isempty(ind)
        tau = 1./k;
        fac = sqrt(tau(1:end-1)./tau(2:end));
        len = length(ind);
        subplot(2,1,2);
        xPatch = reshape([fac(ind).*tau(ind); fac(ind).*tau(ind); ...
            tau(ind)./fac(ind); tau(ind)./fac(ind)], 4*len, 1);
        yPatch = reshape([0*tau(ind); cx(ind+1)'; cx(ind+1)'; ...
            0*tau(ind)], 4*len, 1);
        semilogx(xPatch, yPatch);
        patch(xPatch, yPatch, 'b');
        xlabel('decay time [ns]');
        ylabel('distribution');
    end
end

k = k(:);
offset = cx(1).*M0(1);
cx(1) = [];

if flag > 0
    [cx, k] = clusterDistTailComponents(cx, k, y, dt, bild, N, taumax, T, longFlag);
end
end

function [cx, k] = clusterDistTailComponents(cx, k, y, dt, bild, N, taumax, T, longFlag)
cx = cx(:).';
tmp = cx > 0;
idx = 1:length(tmp);
t1 = idx(tmp(2:end) > tmp(1:end-1)) + 1;
t2 = idx(tmp(1:end-1) > tmp(2:end));

if isempty(t1) || isempty(t2)
    if taumax < 1e6*dt*length(y)
        if longFlag
            [cx, k] = DistTailfit(y, dt, 1, bild, N, 2*taumax, T);
        else
            [cx, k] = DistTailfit(y, dt, 1, bild, N, 2*taumax);
        end
    else
        k = k(:).';
    end
    return;
end

if numel(t1) == 1 && numel(t2) == 1
    k = k(t1);
    cx = cx(t1);
    return;
end

if t1(1) > t2(1)
    t2(1) = [];
end
if ~isempty(t1) && ~isempty(t2) && t1(end) > t2(end)
    t1(end) = [];
end
if length(t1) == length(t2)+1
    t1(end) = [];
end
if length(t2) == length(t1)+1
    t2(1) = [];
end

tmpK = zeros(1, numel(t1));
tmpCx = tmpK;
for jj = 1:numel(t1)
    weights = cx(t1(jj):t2(jj));
    rates = k(t1(jj):t2(jj)).';
    tmpK(jj) = weights*rates/sum(weights);
    tmpCx(jj) = sum(weights);
end
k = tmpK;
cx = tmpCx;
end
