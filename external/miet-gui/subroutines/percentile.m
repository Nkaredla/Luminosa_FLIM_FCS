function [out, out2, bin] = percentile(in,val)

bin = unique(val); % bins of theta
out = zeros(numel(bin),3); out2 = zeros(numel(bin),2);
for i = 1:numel(bin)
    idx = val == bin(i);
    sort_in = sort(in(idx));
    tmp1 = ceil(numel(sort_in)./4); % 25 percent
    tmp2 = ceil(numel(sort_in)*3/4); % 75 percent
    tmp3 = ceil(numel(sort_in)/2); % 50 percent
    out(i,:)  = [sort_in(tmp1) sort_in(tmp2) sort_in(tmp3)];
    out2(i,:) = [mean(sort_in) std(sort_in)];
end