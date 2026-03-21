function err = mlfit(param, irf, y, p)
%	MLFIT(param, irf, y, p) returns the ML error between the data y
%	and the computed values.
%	MLFIT assumes a function of the form:
%
%	  y =  yoffset + A(1)*convol(irf,exp(-t/tau(1)/(1-exp(-p/tau(1)))) + ...
%
%	param(1) is the color shift value between irf and y.
%	param(2) is the irf offset.
%	param(3:...) are the decay times.
%	irf is the measured Instrumental Response Function.
%	y is the measured fluorescence decay curve.
%	p is the time between to laser excitations (in number of TCSPC channels).


% global bild Plothandle
bild = 0;
n = length(irf);
t = 1:n;
tp = (1:p)';
c = param(1);
tau = param(2:length(param)); tau = tau(:)';
x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = Convol(irs, x);
z = [ones(size(z,1),1) z];
A = lsqnonneg(z,y);
z = z*A;
if bild
    % 	set(Plothandle,'ydata',z)
    figure(Plothandle)
    hold all
    plot(t,log10(z));
    hold off
    drawnow
end
ind = y>0;
err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind))/(n-length(tau))



