function [meanbg, minchan, bgpix, smtcspc]= BG_val(tcspc, cutoff, tim, nx, ny)
% Inputs:
% tcspc = combined tcspc curve of the whole image
% tim = number of bins (we need to smooth the tcspc curve)
% nx, ny = number of pixels in x- and y-directio, respectively
%
% Outputs:
% meanbg = estimated background signal of the combined tcspc curve
% minchan = number of the last channel before the signal rises again at the
%           end of the tcspc curve (and therefore last channel used in the
%           calculation of the average arrival time)
% bgpix = estimated number of background counts per pixel (needs nx and ny)
% smtcspc = smoothed tcspc curve

if nargin<3 || isempty(tim)
    tim = 10; 
end

if nargin<4 || isempty(nx)
    flagbg = 1;
elseif nargin<5 || isempty(ny)
    ny = nx;
else
    flagbg = 0;
end

n = ndims(tcspc); % dimensions of tcspc: time channels, detectors, pulses
dim = zeros(n,1);
for i = 1:n
    dim(i) = size(tcspc,i);
end

% Shift dimensions such that the first dimension of tcspc is the time channels
% (assuming that the number of time channels is greater than the number of
% detectors or pulses per sync).
[NChannels,maxdim] = max(dim); 
tcspc = shiftdim(tcspc,maxdim-1);

% Make smoothed tcspc curve "smtcspc" with bin edges "window", starting at peak.
for i = 1:n
    dim(i) = size(tcspc,i);
end
smtcspc = zeros((tim),dim(2:end));

[~,pos_of_max]=max(tcspc);
window = pos_of_max:round((NChannels-pos_of_max+1)/(tim)):NChannels;
if (NChannels-window(end))/round(NChannels/(tim-1))>0.4
    window = [window,NChannels];
end

for k = 1:numel(window)-1
    for i = 1:dim(2)
        if n ==3
            for j = 1:dim(3)
                smtcspc(k,i,j) = mean(tcspc(window(k):window(k+1),i,j));
            end
        else
            smtcspc(k,i) = mean(tcspc(window(k):window(k+1),i));
        end
    end
end

% For each pair of (detector,pulse), find minchan, meanbg and bgpix.
if n ==3
    for i = 1:dim(2)
        for j = 1:dim(3)
            minchan(i,j) = firstmin(smtcspc(:,i,j));                      %#ok<*AGROW>
            meanbg(i,j) = smtcspc(minchan(i,j),i,j);
            minchan(i,j) = window(minchan(i,j));
            if ~flagbg
                bgpix(i,j) = meanbg(i,j)/nx/ny * (minchan(i,j)-cutoff); % background photons per bin * number of bins
            else
                bgpix = NaN;
            end
            
        end
    end
else
    bgpix=zeros(dim(2),1);
    for i = 1:dim(2)
        minchan(i) = firstmin(smtcspc(:,i));   
        meanbg(i)  = smtcspc(minchan(i),i);
        minchan(i) = window(minchan(i));
        if ~flagbg
            bgpix(i)   = meanbg(i)/nx/ny * (minchan(i)-cutoff); % background photons per timebin * number of bins
        else
            bgpix = NaN;
        end
    end
end


