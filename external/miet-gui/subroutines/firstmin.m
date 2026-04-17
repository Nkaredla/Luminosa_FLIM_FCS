function idx = firstmin(array1D)
dip = min(find(diff(array1D)>0)); %#ok<MXFND>

if isempty(dip)
    idx = numel(array1D);
else
    idx = dip;
end