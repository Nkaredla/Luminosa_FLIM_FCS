function back = DeadTimeCorrectionC(meas,deadtime,apd,epsilon)
% program for recovering an unbiased decay-curve form a measured TCSPC
% curve

% input parameters
% meas - measured decay curve
% deadtime - electroncis dead-time in time.units of TCSPC channels
% apd - detector dead-time in time.units of TCSPC channels
% eps - average number of hitting photons per excitation cycle

% If the maximum relative error between two iterations drops below 1-e3, we stop iterating.
REL_ERR_THRESH = 1e-3; 
% Maximum number of iterations to perform
MAX_IT = 10; 

circulMat = gallery('circul',circshift(meas,deadtime+1))'; % Temporary matrix
shiftidx_down = [numel(meas); (1:numel(meas)-1).']; % indexing Nx1 vector with shiftindex_down identical to circshift(vec,1)
shiftidx_up   = [2:numel(meas), 1]; % indexing Nx1 vector with shiftindex_down identical to circshift(vec,-1)

back = meas(:); % Take input as row vector
old_back = back;

for iteration=1:MAX_IT
%     -- Equivalent MATLAB code for the following MEX call: --
%     boxcarSet_fft = fft(gallery('triw',numel(back),1,numel(back)));
%     mm = real(ifft(repmat(fft(back),[1,numel(back)]).*boxcarSet_fft));
    mm = circConv_BoxcarSet(back);
    mm = mm(shiftidx_down,:);
    
    if apd==0
        ww = [ones(numel(meas),1), exp(-mm(:,1:end-1))]/(1-exp(-epsilon));
    else
%     -- Equivalent MATLAB code for the following MEX call: --
%         ww = [repmat(exp(-mm(:,apd)),1,apd-1) + exp(-mm(:,1:apd-1))/(exp(epsilon)-1), exp(-mm(:,apd:end))/(1-exp(-epsilon))];
        ww = computeWW(mm,apd,epsilon);
    end
    
    ww = sum(circulMat.*ww,2);
    back = meas./(ww+ww(shiftidx_up))*2;
    back = back/sum(back)*epsilon;
    
    % Convergence check
    % If the maximum relative error drops below 1-e3, we stop iterating.
    relative_error = max(abs(back-old_back)./old_back);
    if(relative_error< REL_ERR_THRESH)
        return
    end
    old_back = back;
end

end



