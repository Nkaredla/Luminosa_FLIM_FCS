function [h, xv, yv, zv] = mHist3(x,y,z,xv,yv,zv)

x = x(:);
y = y(:);
z = z(:);
ind = ~isfinite(x) | ~isfinite(y) | ~isfinite(z);
x(ind) = []; 
y(ind) = [];
z(ind) = [];
if nargin==4
    if size(xv)==1
        nx = xv; ny = xv; nz = zv;
    else
        nx = xv(1); ny = xv(2); nz = xv(3);
    end
else
    nx = 100; ny = 100; nz = 100;
end
if nargin<6
    xmin = min(x);
    xmax = max(x);
    dx = (xmax-xmin)/nx;
    ymin = min(y);
    ymax = max(y);
    dy = (ymax-ymin)/ny;
    zmin = min(z);
    zmax = max(z);
    dz = (zmax-zmin)/nz;
    xv = xmin:dx:xmax;
    yv = ymin:dy:ymax;
    zv = zmin:dz:zmax;
    x = round((x-xmin)/dx)+1;
    y = round((y-ymin)/dy)+1;
    z = round((z-zmin)/dz)+1;
    xmax = round((xmax-xmin)/dx)+1;
    ymax = round((ymax-ymin)/dy)+1;
else
    xmin = xv(1);
    xmax = xv(end);
    ymin = yv(1);
    ymax = yv(end);
    zmin = zv(1);
    zmax = zv(end);
    x(x>xmax) = xmax;
    x(x<xmin) = xmin;
    y(y>ymax) = ymax;
    y(y<ymin) = ymin;
    z(z>zmax) = zmax;
    z(z<zmin) = zmin;
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
        ymax = round((ymax-ymin)/dy)+1;
    else
        y = round(interp1(yv,1:length(yv),y));
        ymax = round(interp1(yv,1:length(yv),ymax));
    end
    if sum(diff(diff(zv)))==0
        dz = zv(2)-zv(1);
        z = round((z-zmin)/dz)+1;
    else
        z = round(interp1(zv,1:length(zv),z));
    end
end
h = zeros(length(xv)*length(yv)*length(zv),1);
num = sort(x + xmax*(y-1)+xmax*ymax*(z-1));
h(num) = 1;
tmp = diff(diff([0; num; 0])==0);
ind = (1:length(num))';
h(num(tmp==1)) = h(num(tmp==1))-ind(tmp==1)+ind(tmp==-1);
h = reshape(h,length(xv),length(yv),length(zv));
