function [tau,c,z,zz,chi2,p_init,count] = MLEfitTail(y,t,tol,steps)
%[tau,c,z,zz,chi2,p_init] = MLEfitTail(t,y) Summary of this function goes here
%   Detailed explanation goes here
% y         TCSPC data
% t         time
% 
% Sebastian Isbaner, 07/2018

if nargin<2 || isempty(t)
    t = 1:length(y);
end
if nargin<3 || isempty(tol)
    tol=[];
end
if nargin<4 || isempty(steps)
    steps=[];
end

% inittial guess
dt = t(2)-t(1);
[cx, k, offset] = DistTailfit(y, dt,1,0,5,5);
p_init = [offset*numel(y) sum(cx)./sum(cx.*k) sum(cx) ];
count=0;
while true && count<10
        [p_fit,~,stps] = Simplex(@ExpFunMLE_full,p_init,[],[],tol,steps,t,y);
        if stps<steps; break; else count = count+1; end
end
        

[err, c, z, zz, residues] = ExpFunMLE_full(p_fit, t, y);

tau = p_fit(2);
chi2 = nansum(residues.^2)./(numel(y)-numel(p_init));

end

