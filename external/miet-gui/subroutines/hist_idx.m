function [N X idx] = hist_idx(y,nbins)

y = y(:);
if nargin<2|| isempty(nbins)
    nbins=10; % default number of bins
end
[N X] = hist(y,nbins);

idx = false([length(y),numel(X)]);
Xtmp = [-inf X inf];
for i = 1:numel(X)
    [idx(:,i)] = y<=Xtmp(i+1)& y>Xtmp(i);
end
    