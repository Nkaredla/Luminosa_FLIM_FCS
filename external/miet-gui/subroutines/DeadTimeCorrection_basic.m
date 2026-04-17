function back = DeadTimeCorrection(meas,deadtime,apd,eps)

% program for recovering an unbiased decay-curve form a measured TCSPC
% curve

% input parameters
% meas - measured decay curve
% deadtime - electroncis dead-time in time.units of TCSPC channels
% apd - detector dead-time in time.units of TCSPC channels
% eps - average number of hitting photons per excitation cycle

back = meas;
cons = fft(gallery('triw',numel(meas),1,numel(meas)));
for jj=1:10
    mm = real(ifft(repmat(fft(back),[1,numel(meas)]).*cons));
    if apd==0
        ww = exp(-[zeros(numel(meas),1) mm(:,1:end-1)])/(1-exp(-eps));
    else
        ww = [exp(-repmat(mm(:,apd),1,apd-1)) + exp(-mm(:,1:apd-1))/(exp(eps)-1) exp(-mm(:,apd:end))/(1-exp(-eps))];
    end
    ww = circshift(ww,1);
    ww = sum(gallery('circul',circshift(meas,deadtime+1))'.*ww,2);
    back = meas./(0.5*ww+0.5*circshift(ww,-1));
    back = back/sum(back)*eps;
end

