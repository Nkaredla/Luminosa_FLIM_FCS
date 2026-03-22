function [z, xv, yv] = mHist2(x,y,xv,yv,weight)

x = x(:);
y = y(:);
ind = ~isfinite(x) | ~isfinite(y);
x(ind) = []; 
y(ind) = [];
if nargin==3
   if size(xv)==1
      nx = xv; ny = xv;
   else
      nx = xv(1); ny = xv(2);
   end
 else
   nx = 100; ny = 100;
end
if nargin<4 || (isempty(xv) && isempty(yv))
   xmin = min(x);
   xmax = max(x);
   dx = (xmax-xmin)/nx;
   ymin = min(y);
   ymax = max(y);
   dy = (ymax-ymin)/ny;
	xv = xmin:dx:xmax;
	yv = ymin:dy:ymax;
	x = round((x-xmin)/dx)+1;
   y = round((y-ymin)/dy)+1;
   xmax = round((xmax-xmin)/dx)+1;
else
   xmin = xv(1);
   xmax = xv(end);
   ymin = yv(1);
   ymax = yv(end);
   ind = x>xmax | x<xmin | y>ymax | y<ymin;
   x(ind) = [];
   y(ind) = [];
   if sum(diff(diff(xv)))==0
      dx = xv(2)-xv(1);
      x = round((x-xmin)/dx)+1;
      xmax = round((xmax-xmin)/dx)+1;
   else
      x = round(interp1(xv,1:length(xv),x));
      xmax = round(interp1(xv,1:length(xv),xmax));
   end
   if sum(diff(diff(yv)))==0
      dy = yv(2)-yv(1);
      y = round((y-ymin)/dy)+1;
   else
      y = round(interp1(yv,1:length(yv),y));
   end
end
z = zeros(length(xv)*length(yv),1);
if nargin<5
    num = sort(x + xmax*(y-1));
    z(num) = 1;
    tmp = diff(diff([0; num; 0])==0);
    ind = (1:length(num))';
    z(num(tmp==1)) = z(num(tmp==1))-ind(tmp==1)+ind(tmp==-1);
else
    num = x + xmax*(y-1);
    [num, ord] = sort(num);
    weight = weight(:);
    weight = weight(ord);
    tmp = diff([0; num])>0;
    z(num(tmp)) = weight(tmp);    
    tmp = diff(diff([0; num; 0])==0);
    ind = cumsum(weight);
    z(num(tmp==1)) = ind(tmp==-1)-ind(tmp==1)+weight(tmp==1);
end
z = reshape(z,length(xv),length(yv)).';
if nargout==0
   pcolor(xv,yv,log(z)); 
   axis image;
	set(get(gca,'children'),'edgecolor',[0.5 0.5 0.5]);
   clear z xv yv
end
