function MIET_lifedist(t,tcspc,tau_min, tau_max)

t = t(:);
N = 200;
n = numel(t);
k = linspace(1/tau_max,1/tau_min,N); % distribution of inverse decay times

M0 = [ones(size(t)) exp(-t*k)];

M0 = M0./(ones(n,1)*sum(M0)); % each curve normalized to area=1

cx = lsqnonneg(M0,tcspc); % find amplitudes of normalized exponential functions
z = M0*cx;            % fitted curve
err = sum((z-tcspc).^2./abs(z))/n;