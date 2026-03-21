function [auto, autotime] = tttr2xfcsSym(y, num, Ncasc, Nsub)
% -------------------------------------------------------------------------
% Symmetrically normalised auto/cross FCS:
%   g_s(k) =  G(k) / ( n̄_left · n̄_right )  – 1
%
%  – y      :  photon arrival times (integer clock ticks)
%  – num    :  logical/weighting matrix – which photon belongs to which stream
%  – Ncasc  :  number of cascades     (bin width doubles each cascade)
%  – Nsub   :  channels per cascade   (e.g. 8 → “multiple-tau” spacing)
% -------------------------------------------------------------------------

y = round(y(:));
if size(num,1) ~= numel(y)
    num = num.';   % force rows = photons
end
nch = size(num,2);

autotime = zeros(Ncasc*Nsub,1);
auto     = nan(Ncasc*Nsub,nch,nch);

shift = 0;
delta = 1;

for j = 1:Ncasc

    % --- correct binning: sum counts for identical bin indices ---
    [y,k] = unique(y,'stable');
    tmp   = cumsum(num);
    num   = diff([zeros(1,size(num,2)); tmp(k,:)]);


    ymin = min(y);
    ymax = max(y);
    nbins_full = ymax - ymin + 1;

    for ks = 1:Nsub
        idx = ks + (j-1)*Nsub;

        shift = shift + delta;
        lag   = round(shift/delta);

        autotime(idx) = lag * delta;

        if lag <= 0 || lag >= nbins_full
            continue;
        end

        npair_tot = nbins_full - lag;

        % means over the full overlap windows (include zeros via division)
        maskL = (y >= ymin)       & (y <= ymax - lag);
        maskR = (y >= ymin + lag) & (y <= ymax);

        nbarL = sum(num(maskL,:),1) / npair_tot;
        nbarR = sum(num(maskR,:),1) / npair_tot;

        den = (nbarL.' * nbarR);
        den(den==0) = NaN;  % avoid Inf when a channel is empty

        % numerator: only bins where both sides are populated contribute
        [~,iL,iR] = intersect(y, y+lag);
        if isempty(iL)
            G = zeros(nch,nch);
        else
            kL = num(iL,:);
            kR = num(iR,:);
            G  = (kL.' * kR) / npair_tot;
        end

        auto(idx,:,:) = G ./ den;   % ELEMENT-WISE normalisation
    end

    % next cascade: bin width doubles
    y     = floor(0.5 * y);
    delta = 2 * delta;
end
end
