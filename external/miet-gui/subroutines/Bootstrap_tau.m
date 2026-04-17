function [tcspc, tau, dtau] = Bootstrap_tau(bin,microtime,shift,num,n)

photons = numel(microtime);
if nargin<5||isempty(n)
    n = 100; % number of bootstrap calculations
end
if nargin<4||isempty(num)
    num = 10000; % 1000 photons per bin
end
if nargin<3||isempty(shift)
    shift = 10; % cutoff from the maximum of the peak
end
if photons<num
    disp('Too few photons for bootstrapping')
    tcspc = []; tau = [];
    return
end

Ngate = numel(bin);
tcspc = mHist(double(microtime),1:Ngate);
[~,pos] = max(tcspc);

tmptau = zeros(n,1);
for i = 1:n
    idx = randi([1 photons],[num,1]);
    tmptcspc = mHist(double(microtime(idx)),1:Ngate);
    p = [1 2];
    p = Simplex('ExpFun',p,[],[],1e-10,Inf,bin(pos+shift:end),tmptcspc(pos+shift:end),[],[],1);
    [~, c, ~, ~] = ExpFun(p,bin(pos+shift:end),tmptcspc(pos+shift:end),[],[],1);
    tmptau(i) = sum(c(2:end))./sum(1./p.*c(2:end));
end
tau = mean(tmptau); dtau = std(tmptau);
