function [err, c, zz, z] = TCSPC_Fun(p, t, y, para)

%  Fits the following model function to TCSPC-Data
%
%  TCSPC(N,{q_1},{tau_i},t_0, w, a, dt, T1, T2, bg, t) = 
%   bg + N/Q Sum[ q_i (sqrt(pi/2) w/tau_i exp(w^2/(2 tau_i^2) - t'/tau_i) Erfc((w^2- t' tau_i)/(sqrt(2)w tau_t)) +
%                 a P H(t')) (tau_i/(tau_i/T12 -1)/(tau_1/T1 -1)/T2 exp(-t''\tau_i) + 
%                                 1/(tau_1/T12 -1) exp(-t''/T12) - 
%                                 1/(tau_i/T1  -1) exp(-t''/T1) ) ]
%   t'  = t - t_0
%   t'' = t - t_0 - dt;
%   Q   = sqrt(2 pi w^2) + a T1(1+T1/T2)^{T2/T1}
%   P   = (1-T2/T1)(1+T1/T2)^{T2/T1}
%
%  (see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)
%
%   The program needs the following m-file: Convol.m.

p = p(:);

if (nargin>3)&&(~isempty(para))
    para = para(:);    
    n = numel(para);
    if n > 6 
        p = [para; p];
    else
        p = [p(1:7-n); para; p(8-n:end)];
    end
end


if length(t)<length(y)
    
    c = t;
    t = y(:);

    nex = numel(p) - 7;
    
    tau = p(8:end);
    
    IRF = IRF_Fun(p(1:7), t);
  
    zz = zeros(numel(t),nex+1);
    zz(:,1) = ones(size(t));
    for i = 1:nex
      tmp = Convol(IRF, exp(- t1./tau(i))./tau(i));
      zz(:,1+i) =  tmp(1:numel(t));
    end;

%     if nargin>4 && ~isempty(nrm) 
%         zz = zz./(ones(length(t),1)*sum(zz));
%     end
    
    for j=1:size(c,2)
        err(:,j) = zz*c(:,j);
    end

else

    t = t(:);

    [m, n] = size(y);
    if m<n 
        y   = y'; 
        [m, n] = size(y);
    end
    
	t = t(isfinite(sum(y,2)));
    y = y(isfinite(sum(y,2)),:);
        
    nex = numel(p) - 7;
          
    IRF = IRF_Fun(p(1:7), t);

    tau = p(8:end);
    t1  = t-p(1);
    
    zz = zeros(numel(t),nex+1);
    zz(:,1) = ones(size(t));
    for i = 1:nex
      tmp = Convol(IRF, exp(- t1./tau(i))./tau(i));
      zz(:,1+i) =  tmp(1:numel(t));
    end;
        
%     if nargin>5 && ~isempty(nrm) 
%         zz = zz./(ones(m,1)*sum(zz));
%     end

    for j=1:n
        c(:,j) = lsqnonneg(zz,y(:,j));
%         c(:,j) = zz\y(:,j);
        z(:,j) = zz*c(:,j);
    end
  
    err = sum(sum((y-z).^2./(10+abs(z))));
end

