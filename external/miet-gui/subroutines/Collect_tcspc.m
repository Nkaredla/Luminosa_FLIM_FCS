function [tcspc_sm, photons_sm] = Collect_tcspc(tcspc_pix,field)

num = size(field,3);
Nchan = size(tcspc_pix,3);
if size(tcspc_pix,4)>1
    tcspc_pix=sum(tcspc_pix,4);
end
tcspc_sm = zeros(Nchan,num);

tcspc_pix = reshape(tcspc_pix,[size(tcspc_pix,1)*size(tcspc_pix,2),Nchan]);
for i = 1:num
    pix_ind = reshape(field(:,:,i),[size(tcspc_pix,1),1]);
    tcspc_sm(:,i) = sum(tcspc_pix(logical(pix_ind),:),1);
end
photons_sm = sum(tcspc_sm,1);
