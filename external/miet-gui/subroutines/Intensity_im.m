function [tag tcspc maxch] = Intensity_im(im_chan, im_line, im_col, im_tcspc, cutoff)

nx = double(max(im_col)); ny = double(max(im_line)); 
resolution = 32e-3; % resolution in nanoseconds
shift = ceil(cutoff./resolution);
Ngate = double(max(im_tcspc));

maxch = unique(im_chan);

for ch = 1:numel(maxch)
    tcspc(ch,:) = mHist(double(im_tcspc(im_chan == maxch(ch))),1:Ngate);
    indel(ch) = sum(tcspc(ch,:))<500; %#ok<*AGROW> % deleting image planes containing less than 500 photons
end
tcspc(indel,:) = [];
maxch(indel)    = [];
maxch_n        = numel(maxch);

[~,pos] = max(sum(tcspc,1));
ind = im_tcspc<(pos+shift);
im_chan(ind) =[];
im_line(ind) =[];
im_col(ind)  =[];

im_pixel=zeros(size(im_line));
for i=1:maxch_n
    ind = im_chan==maxch(i);
    im_pixel(ind)=double(im_line(ind))+(double(im_col(ind))-1)*ny + nx*ny*(i-1);
end

tag=zeros(ny*nx*maxch_n,1);
tag=reshape(mHist(im_pixel,1:numel(tag)),[ny,nx,maxch_n]);
 