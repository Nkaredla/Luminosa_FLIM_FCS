function s = mnum2str(x,m,n)
% mnum2str(x,m,n) converts the number x into a string with m digits before and n digits after the period, 
% patching with blanks and zeros if necessary

if nargin==1
    m=-1;
    n=-1;
end
if nargin==2
    n=floor(log10(x)+1);
    n(n<0)=1;
    n=m-n;
    m=m-n;
end
if nargin==3
    s = num2str(round(x*10^n)/10^n);
else
    s = num2str(x);
end
pos = findstr(s,'.');
if isempty(pos)
    l2 = 0;
    l1 = length(s);
    s = [s '.'];
else
    l2 = length(s)-pos;
    l1 = length(s)-l2-1;
end
if m>0
    for j=1:m-l1
        s = [' ',s];
    end
end
if n>0
    for j=1:n-l2
        s = [s,'0'];
    end
end

