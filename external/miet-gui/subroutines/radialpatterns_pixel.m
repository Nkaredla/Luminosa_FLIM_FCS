function [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c, mask, molim] = radialpatterns_pixel(tag, pattern, pixel, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic)
% radialpatterns_pixel computes patterns for the givel pixel size and the
% excitation wavelength and then identifies molecules by calling
% FindPatternRapid.m
if nargin<14 || isempty(pic)
    pic=0;
end

if nargin<13 || isempty(focpos)
    focpos=0;
end
pixel = pixel.*(0.9:0.2/5:1.1);
close all
nn = round(0.4./pixel);
resolution = [lamex/0.02 lamex/0.001];
rhofield = [-lamex/resolution(1)/2 nn(1)*pixel(1)*1.1];
zfield = [0 0.01];
fd = 3e3;

over = inf;
maxm = 3;
atf  = [];
ring = 'cos(psi).*rad';



[fxc, fxs, fyc, fys, fzc, fzs, rho, ~] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);            %#ok<*AGROW>
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;

mask = zeros(2*max(nn)+1,2*max(nn)+1,numel(theta)*numel(nn));
c=0;
for i = 1:numel(pixel)
    [xx,yy] = meshgrid(-nn(i):nn(i),-nn(i):nn(i));
    rr = pixel(i)*sqrt(xx.^2+yy.^2);
    psi = angle(xx + 1i*yy);
    
    for k=1:numel(theta)
        c=c+1;
        mask(max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,c) = interp1(rho(:,1),(fxc(:,1,1)*sin(theta(k))*cos(phi(k))+fyc(:,1,1)*sin(theta(k))*sin(phi(k))+fzc(:,1,1)*cos(theta(k))) ,rr,'cubic',0);
        for j=1:maxm
            mask(max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,c) = mask(max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,max(nn)-nn(i)+1:max(nn)-nn(i)+2*nn(i)+1,c) + interp1(rho(:,1),(fxc(:,1,j+1)*sin(theta(k))*cos(phi(k))+fyc(:,1,j+1)*sin(theta(k))*sin(phi(k))+fzc(:,1,j+1)*cos(theta(k))) ,rr,'cubic',0).*cos(j*psi) + ...
                interp1(rho(:,1),(fxs(:,1,j)*sin(theta(k))*cos(phi(k))+fys(:,1,j)*sin(theta(k))*sin(phi(k))+fzs(:,1,j)*cos(theta(k))) ,rr,'cubic',0).*sin(j*psi);
        end
    end
end
mask = abs(mask).^2;

if strcmpi(pattern,'radial')
    tresh = 1;
else
    tresh = 0.8;
end

[err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, ~] = FindPatternRapid(tag,mask,mask(:,:,1)>0,mask(:,:,1)>0,[],tresh,[],0);

scmod  = ceil(sc/numel(theta));
modesc = mode(scmod);
delind = scmod<modesc-1 | scmod>modesc+1; 
xc(delind)  =[];
yc(delind)  =[];
bc(delind)  =[];
cc(delind)  =[];
sc(delind)  =[];
len(delind) =[];
imm(:,:,delind) =[];
imm_c = sum(imm,3);



if pic==1
    figure
    CombineImages(cat(3,tag,imm_c),1,2)
end

if ~isempty(len)
    nn = max(nn);
    nx = 2*round(sqrt(numel(len)/2));
    ny = ceil(numel(len)/nx);
    molim = zeros(2*(2*nn+1)*ny,(2*nn+1)*nx);
    tag_ext=zeros(size(tag,1)+2.*nn+1,size(tag,2)+2.*nn+1);
    tag_ext(nn+1:end-nn-1,nn+1:end-nn-1)=tag;
    imm_c_ext=zeros(size(imm_c,1)+2.*nn+1,size(imm_c,2)+2.*nn+1);
    imm_c_ext(nn+1:end-nn-1,nn+1:end-nn-1)=imm_c;
    for jx=1:nx
        for jy=1:ny
            if (jy-1)*nx+jx<=numel(len)
                j = (jy-1)*nx+jx;
                molim((2*jy-2)*(2*nn+1)+(1:(2*nn+1)),(2*nn+1)*(jx-1)+1:(2*nn+1)*jx) = tag_ext(yc(j)+(0:2*nn),xc(j)+(0:2*nn));
                molim((2*jy-1)*(2*nn+1)+(1:(2*nn+1)),(2*nn+1)*(jx-1)+1:(2*nn+1)*jx) = imm_c_ext(yc(j)+(0:2*nn),xc(j)+(0:2*nn));
            end
        end
    end
else
    molim=zeros(size(tag));
end

end