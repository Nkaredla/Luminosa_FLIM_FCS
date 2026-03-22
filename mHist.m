function [z, xv] = mHist(x,xv,weight)

if size(x,1)==1 || size(x,2)==1
    x = x(:);
    if nargin>2
        weight = weight(:);
        weight(~isfinite(x)) = [];
    end
    x(~isfinite(x)) = [];
end
if nargin>1 && ~isempty(xv)
    if size(xv)==1
        n = xv;
    else
        n = [];
        xmin = xv(1);
        xmax = xv(end);
        ind = x>xmax | x<xmin;
        x(ind) = [];
        if nargin>2
            weight(ind) = [];
        end
        if length(unique(diff(xv)))==1
            dx = xv(2)-xv(1);
            x = round((x-xmin)/dx)+1;
            xmax = round((xmax-xmin)/dx)+1;
        else
            x = round(interp1(xv,1:length(xv),x));
            xmax = round(interp1(xv,1:length(xv),xmax));
        end
    end
else
    n = 100;
end
if ~isempty(n)
    xmin = min(min(x));
    xmax = max(max(x));
    dx = (xmax-xmin)/n;
    xv = xmin:dx:xmax;
    x = round((x-xmin)/dx)+1;
    xmax = round((xmax-xmin)/dx)+1;
end   
z = zeros(length(xv)*size(x,2),1);
if nargin<3
    num = sort(x+xmax*ones(size(x,1),1)*(0:size(x,2)-1));
    num = num(:);
    z(num) = 1;
    tmp = diff(diff([0; num; 0])==0);
    ind = (1:length(num))';
    z(num(tmp==1)) = z(num(tmp==1))-ind(tmp==1)+ind(tmp==-1);
else
    num = x+xmax*ones(size(x,1),1)*(0:size(x,2)-1);
    num = num(:);
    [num, ord] = sort(num);
    weight = weight(:);
    weight = weight(ord);
    tmp = diff([0; num])>0;
    z(num(tmp)) = weight(tmp);    
    tmp = diff(diff([0; num; 0])==0);
    ind = cumsum(weight);
    z(num(tmp==1)) = z(num(tmp==1)) + ind(tmp==-1)-ind(tmp==1); 
end    
z = reshape(z,length(xv),size(x,2));
if nargout==0
    bar(xv,z,1,'b'); 
    clear z xv
end
