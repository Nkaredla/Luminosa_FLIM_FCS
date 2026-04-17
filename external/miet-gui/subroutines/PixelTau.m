function [tau, tau_c] = PixelTau(tcspc,sync,t0,deadtime,apd,dt,tsync)

if nargin<7 || isempty(tsync)
    tsync = 50; % 50 ns
end
if nargin<6 || isempty(dt)
    dt = 0.064; % 32 ps
end
if nargin<5 || isempty(apd)
    apd = 33; % 33ns;
end
if nargin<4 || isempty(deadtime)
    deadtime = 80; % 80 ns
end

apd = ceil(apd./dt); deadtime = ceil(deadtime./dt); % in number of bins

t = sync*ceil(tsync./dt)+tcspc; % in tcspc bins
meas = mHist(tcspc,1:ceil(tsync/dt)); % the tcspc of the data
ipd = mHist(diff(t),1:10*ceil(tsync/dt)+apd+deadtime); % bins of tcspc bins
N = 1:10;
for i = 1:10
    range = (i-1)*ceil(tsync./dt)+1:i*(tsync./dt); range = range+apd+deadtime;
    mN(i) = sum(ipd(range));
end

epsilon = Simplex('ExpFun',1e3,[],[],[],[],N(mN>0),mN(mN>0)./tsync,1,[],1);
epsilon = 1./epsilon;

%  c = polyfit(N(mN>0),log(nphot(mN>0)./tsync),1);
% c = -c(1);
% plot(N,exp(polyval(c,N)))

meas = meas./sum(meas);
back = meas;
cons = fft(gallery('triw',numel(meas),1,numel(meas))); cons2 = gallery('circul',circshift(meas,deadtime+1))';
for jj=1:10
    mm = real(ifft(repmat(fft(back),[1,numel(meas)]).*cons));
    if apd==0
        ww = exp(-[zeros(numel(meas),1) mm(:,1:end-1)])/(1-exp(-epsilon));
    else
        ww = [exp(-repmat(mm(:,apd),1,apd-1)) + exp(-mm(:,1:apd-1))/(exp(epsilon)-1) exp(-mm(:,apd:end))/(1-exp(-epsilon))];
    end
    ww = circshift(ww,1);
    ww = sum(cons2.*ww,2);
    back = meas./(0.5*ww+0.5*circshift(ww,-1));
    back = back/sum(back)*epsilon;
end

% p = [1 2];
% p = Simplex('ExpFun',p,[],[],[],[],[t0:ceil(tsync./dt)].*dt,meas(t0:end),[],[],1);
[c, p, ~, ~, ~, ~] = DistTailfit(meas(t0:end), dt,1);
% [~,c] = ExpFun(p,[t0:ceil(tsync./dt)].*dt,meas(t0:end),[],[],1);
tau = sum(c)/sum(c.*p);

% p = Simplex('ExpFun',p,[],[],[],[],[t0:ceil(tsync./dt)].*dt,back(t0:end),[],[],1);
[c, p, ~, ~, ~, ~] = DistTailfit(back(t0:end), dt);
% [~,c] = ExpFun(p,[t0:ceil(tsync./dt)].*dt,back(t0:end),[],[],1);
tau_c = sum(c)/sum(c.*p);
