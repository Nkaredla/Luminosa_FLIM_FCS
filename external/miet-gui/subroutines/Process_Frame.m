function [tag, tau, tcspc_pix, time] = Process_Frame(im_sync,im_col,im_line,im_chan,im_tcspc,head)

Resolution = max([head.MeasDesc_Resolution*1e9 0.256]); % resolution of 1 ns to calculate average lifetimes
chDiv      = ceil(1e-9*Resolution./head.MeasDesc_Resolution);
SyncRate   = 1./head.MeasDesc_GlobalResolution;
nx      = double(max(im_col));
ny      = double(max(im_line));
dind    = double(unique(im_chan));
Ngate   = double(max(im_tcspc./chDiv));
maxch_n = numel(dind);

tcspc_pix = zeros(nx,ny,Ngate,maxch_n);
time = zeros(numel(im_sync),maxch_n);
tag = zeros(nx,ny,maxch_n);
tau = tag;
bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1])*Resolution; % 3D time axis
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    tcspc_pix(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)./chDiv),1:nx,1:ny,1:Ngate); % tcspc histograms for all the pixels at once!
    time(:,ch) = im_sync(ind)*ceil(1/SyncRate/Resolution/1e-9)+double(im_tcspc(ind)); % in tcspc bins
    tag(:,:,ch) = sum(tcspc_pix,3);
    tau(:,:,ch) = real(sqrt((sum(bin.^2.*tcspc_pix,3)./tag)-(sum(bin.*tcspc_pix,3)./tag).^2));   
end
                           