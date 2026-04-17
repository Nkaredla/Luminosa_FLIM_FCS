function [err, c, zz, z] = ExpFunMLE(varargin)
%[err, c, zz, z] = ExpFun(p, t, y, pic, nrm, pos)
% p   - parameters (lifetimes times)
% t   - time (x) axis
% y   - function values (y-axis)
% pic - plot toggle: 0, no plot; 1 plot; 2, semilogy; 3, semilogx
% nrm - normelise the exponetial componets
% pos - toggle for amplitute fitting: 
%           0 - unconstrained fit (by matrix left division); 
%           1 - nonnegative least square fit of componets; 
%           cell array of constrains, min or [min max]
%
% err - error: total relative square
% c   - compontens amplitudes
% zz  - compontens function value
% z   - total function value

[err, c, zz, z] = ExpFun(varargin{:});

if numel(err)==1 % Normal mode, not ExpFun(p,c,t) to calculate curve
    y = varargin{3};
    err = sum(sum(-y.*log(z+1e-30)+z));
end

end

