function [lt_trace Amp Tau] = LT_trace(tcspc_roi_frame,resolution)

if isstruct(tcspc_roi_frame)
    tmp = tcspc_roi_frame;
    tcspc_roi_frame = tmp.tcspc_roi_frame;
    if isfield(tmp, 'Resolution')
        resolution = tmp.Resolution;
    else
        resolution = 0.256; % 256 ps resolution
    end
end

if nargin<2||isempty(resolution)
    resolution = 0.256;
end

[nFrames nChannel nMol]= size(tcspc_roi_frame);
int_trace = squeeze(sum(tcspc_roi_frame,2));
tcspc_roi = squeeze(sum(tcspc_roi_frame));
tcspc_tot = sum(tcspc_roi,2);
% calc IRF
tmp.head.Resolution = resolution;
tmp.head.SyncRate = tmp.head.TTResult_SyncRate;
tcspcIRF = Calc_mIRF(tmp.head, tcspc_tot.');

% [cx, tau, offset, csh, z, t, err] = DistFluofit_extension(tcspcIRF, squeeze(sum(sum(tcspc_roi_frame(1:10,:,:),3),1)), nChannel*resolution, resolution, [-10 10],1,1);

% tresh =

% DistFloutfit calculations that are repetitive

tcspcIRF = tcspcIRF(:);
n = length(tcspcIRF);
tp = resolution*(1:nChannel)';
t = (1:n)';
N = 100; % number of lifetime decay patterns to be calculated
sh_min = -10;
sh_max = 10;

tau = (1/resolution)./exp((0:N)/N*log(nChannel)); % distribution of decay times
M0 = [ones(size(t)) Convol(tcspcIRF,exp(-tp*tau))];
M0 = M0./(ones(n,1)*sum(M0));
% err = [];

lt_trace = zeros(nFrames,nMol);
Amp = cell(nFrames,nMol);
Tau = Amp;

for i = 1:nFrames
    disp(i)
    for j = 1:nMol
        y = squeeze(tcspc_roi_frame(i,:,j));
        if sh_max-sh_min>0
            for c=sh_min:sh_max
                M = (1-c+floor(c))*M0(rem(rem(t-floor(c)-1, n)+n,n)+1,:) + (c-floor(c))*M0(rem(rem(t-ceil(c)-1, n)+n,n)+1,:);
                ind = max([1,1+c]):min([n,n+c]);
                cx = lsqnonneg(M(ind,:),y(ind));
                z = M*cx;
                err = [err sum((z-y).^2./abs(z))/n];
            end
            
            shv = sh_min:0.1:sh_max;
            tmp = interp1(sh_min:sh_max, err, shv);
            [~, pos] = min(tmp);
            csh = shv(pos);
        else
            csh = sh_min;
        end
        
        M = (1-csh+floor(csh))*M0(rem(rem(t-floor(csh)-1, n)+n,n)+1,:) + (csh-floor(csh))*M0(rem(rem(t-ceil(csh)-1, n)+n,n)+1,:);
        c = ceil(abs(csh))*sign(csh);
        ind = max([1,1+c]):min([n,n+c]);
        cx = lsqnonneg(M(ind,:),y(ind));
%         z = M*cx; % Fitted curve
%         err = sum((z-y).^2./abs(z))/n;
        tau = tau';
        cx(1) = [];
        
        if flag>0
            cx = cx';
            tmp = cx>0.2*max(cx);
            t = 1:length(tmp);
            t1 = t(tmp(2:end)>tmp(1:end-1)) + 1;
            t2 = t(tmp(1:end-1)>tmp(2:end));
            try
                if t1(1)>t2(1)
                    t2(1)=[];
                end
                if t1(end)>t2(end)
                    t1(end)=[];
                end
                if length(t1)==length(t2)+1
                    t1(end)=[];
                end
                if length(t2)==length(t1)+1
                    t2(1)=[];
                end
                tmp = []; bla = [];
                for j=1:length(t1)
                    tmp = [tmp cx(t1(j):t2(j))*tau(t1(j):t2(j))/sum(cx(t1(j):t2(j)))];
                    bla = [bla sum(cx(t1(j):t2(j)))];
                end
                cx = bla./tmp;
                cx = cx/sum(cx);
                tau = tmp;
            catch
                tau = 0;
                cx = NaN;
            end
        end
        Amp{i,j} = cx;
        Tau{i,j} = 1./tau;
        lt_trace(i,j) =  sum(cx)./sum(cx.*tau);
    end
end


end
