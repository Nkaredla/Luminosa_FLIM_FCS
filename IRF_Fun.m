function [z] = IRF_Fun(p, t, pic)

%  Computes the following model IRF 
%
%  IRF(t_0, w, a, dt, T1, T2, t) = 
%           1/Q [ exp(- t'^2/(2 w^2) +
%                 a P H(t'') exp(- t''/T1) (1 - exp(-t''/T2)) ] 
%
%   t'  = t - t_0
%   t'' = t - t_0 - dt;
%   Q   = sqrt(2 pi w^2) + a T1(1+T1/T2)^{T2/T1}
%   P   = (1-T2/T1)(1+T1/T2)^{T2/T1}
%
%  (see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)

t = t(:);
p = p'; 

t_0 = p(1);
w1  = p(2);
T1  = p(3);
T2  = p(4);
a   = p(5);
b   = p(6);
dt  = p(7);

t1 = t-t_0;
t2 = t-t_0-dt;
t3 = t+t_0;

H  = ones(size(t));
H(t<(t_0+dt)) = 0;
G  = zeros(size(t));
G(t<(t_0-dt)) = 1;

% IRF =  [exp(- t1.^2/(2*w1)) G.*(exp(t1./T1)) H.*(exp(-t2./T1).*(1-exp(-t2./T2)))];
IRF =  [exp(- t1.^2/(2*w1)) G.*exp(t3./T1) H.*(exp(-t2./T1).*(1-exp(-t2./T2)))];
IRF =  ones(numel(t),1)*[1 b a].*IRF./(ones(numel(t),1)*sum(IRF));

tm = 0.5*max(IRF(:,1));
ind = IRF(:,2)>tm; 
IRF(ind,2) = tm;
ind = IRF(:,3)>tm; 
IRF(ind,3) = tm;
IRF(isnan(IRF)) = 0;
IRF =  sum(IRF,2)./sum(sum(IRF));
IRF(IRF<0) = 0;


[tmp, t0] = max(IRF);
tmp = IRF(1:t0);
tmp = diff(tmp);
tmp(tmp<0) = 0;
tmp = [0; cumsum(tmp)];
IRF(1:t0) = tmp;
z =  sum(IRF,2)./sum(sum(IRF));

if nargin>2 && ~isempty(pic)
    if pic==1
        plot(t, z, 'r'); drawnow
    else
        semilogy(t, z); drawnow
    end
end

