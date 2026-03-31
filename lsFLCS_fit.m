function [err, c, z, zz] = lsFLCS_fit(p,t,x,y,v,pixel,flag)

w0 = p(end);

if strcmpi(flag,'2d_stick') % 2 diffusion coeff (one global) + sticking model
    tauD = p(1:end-3);
    kp = p(end-2);
    km = p(end-1);
    
    
    ytop    = y(:,:,1);
    ybottom = y(:,:,2);
    [xixi, tt] = meshgrid(x*pixel*1e-6,t);
    
    tmp = Stick2DFCSFun(x*pixel,t,w0/1e-6,((w0/1e-6)^2)./(4*tauD(end)),kp,km).';
    ytmp = tmp(:)./max(tmp(:));
    ytmp(:,end+1)=1;
    c1 = lsqnonneg(ytmp,ybottom(:));
    zz1 = cat(3,ones(size(ybottom(:,:,1))), tmp./max(tmp(:)));
    z1 = c1(2)*zz1(:,:,1)+c1(1)*zz1(:,:,2);
    
    % err1 = sum(sum((ytop-z1).^2./(ybottom.^2)));
    err1 = sum(sum((ybottom-z1).^2./abs(z1)));
    
    
    zz = []; ytmp = [];
    for j = 1:numel(tauD)
        tmp = exp(-(xixi/w0).^2./(1+ (tt/tauD(j) + xixi/(v*tauD(j)))))./(1 + tt/tauD(j) + xixi/(v*tauD(j)));
        zz = cat(3,zz,tmp);
        ytmp(:,j)=tmp(:)./max(tmp(:));
    end
    zz = cat(3,ones(size(ybottom)), zz);
    zz = zz./repmat(max(max(zz)),size(tt,1),size(tt,2));
    ytmp(:,end+1)=1;
    c2 = lsqnonneg(ytmp,ytop(:));
    % c = y(:)\ytmp;
    % weight = repmat(ytop(1,:),size(t,1),1);
    
    z2 = c2(3)*zz(:,:,1)+c2(1)*zz(:,:,2)+c2(2)*zz(:,:,3);
    % z = c(1)*zz(:,:,1)+c(2)*zz(:,:,2);
    
    % err2 = sum(sum(((ytop-z2)).^2));
    err2 = sum(sum((ytop-z2).^2./abs(z2)));
    
    err = err1+err2;
    disp([err1 err2 err])
    z = cat(3,z1,z2);
    
    c = [c1;c2];
    
elseif strcmpi(flag,'diff_mult') % multiple diffusion (all global) model
    tauD = p(1:end-1);
    ytop    = y(:,:,1);
    ybottom = y(:,:,2);
    [xixi, tt] = meshgrid(x*pixel*1e-6,t);
    
    zz = []; ytmp = [];
    for j = 1:numel(tauD)
        tmp = exp(-(xixi/w0).^2./(1+ (tt/tauD(j) + xixi/(v*tauD(j)))))./(1 + tt/tauD(j) + xixi/(v*tauD(j)));
        zz = cat(3,zz,tmp);
        ytmp(:,j)=tmp(:)./max(tmp(:));
    end
    zz = cat(3, zz, ones(size(ytop)));
    zz = zz./repmat(max(max(zz)),size(t,1),size(t,2));
    ytmp(:,end+1)=1;
    c1 = lsqnonneg(ytmp,ytop(:));
    % c = y(:)\ytmp;
    % weight = repmat(ytop(1,:),size(t,1),1);
    z1 = 0;
    for j = 1:size(zz,3)
        z1 = z1+c1(j)*zz(:,:,j);
    end
    
    % z1 = c1(3)*zz(:,:,1)+c1(1)*zz(:,:,2)+c1(2)*zz(:,:,3);
    
    err1 = sum(sum((ytop-z1).^2./abs(z1)));
    
    
    zz = []; ytmp = [];
    for j = 1:numel(tauD)
        tmp = exp(-(xixi/w0).^2./(1+ (tt/tauD(j) + xixi/(v*tauD(j)))))./(1 + tt/tauD(j) + xixi/(v*tauD(j)));
        zz = cat(3,zz,tmp);
        ytmp(:,j)=tmp(:)./max(tmp(:));
    end
    zz = cat(3, zz, ones(size(ybottom)));
    zz = zz./repmat(max(max(zz)),size(t,1),size(t,2));
    ytmp(:,end+1)=1;
    c2 = lsqnonneg(ytmp,ybottom(:));
    % c = y(:)\ytmp;
    % weight = repmat(ytop(1,:),size(t,1),1);
    z2 = 0;
    for j = 1:size(zz,3)
        z2 = z2+c2(j)*zz(:,:,j);
    end
    
    % z2 = c2(3)*zz(:,:,1)+c2(1)*zz(:,:,2)+c2(2)*zz(:,:,3);
    % z = c(1)*zz(:,:,1)+c(2)*zz(:,:,2);
    
    % err2 = sum(sum(((ytop-z2)).^2));
    err2 = sum(sum((ybottom-z2).^2./abs(z2)));
    
    err = err1+err2;
    disp([err1 err2 err])
    z = cat(3,z1,z2);
    
    c = [c1;c2];
    
elseif strcmpi(flag,'diff_intermediate') % 3 diffusion (one intermediate / global) model
    tauD = p(1:end-1);
    ytop    = y(:,:,1);
    ybottom = y(:,:,2);
    [xixi, tt] = meshgrid(x*pixel*1e-6,t);
    
    zz = []; ytmp = [];
    for j = 1:numel(tauD)-1
        tmp = exp(-(xixi/w0).^2./(1+ (tt/tauD(j) + xixi/(v*tauD(j)))))./(1 + tt/tauD(j) + xixi/(v*tauD(j)));
        zz = cat(3,zz,tmp);
        ytmp(:,j)=tmp(:)./max(tmp(:));
     end
    zz = cat(3, zz, ones(size(ytop)));
    zz = zz./repmat(max(max(zz)),size(t,1),size(t,2));
    ytmp(:,end+1)=1;
    c1 = lsqnonneg(ytmp,ytop(:));
    % c = y(:)\ytmp;
    % weight = repmat(ytop(1,:),size(t,1),1);
    z1 = 0;
    for j = 1:size(zz,3)
        z1 = z1+c1(j)*zz(:,:,j);
    end
    
    % z1 = c1(3)*zz(:,:,1)+c1(1)*zz(:,:,2)+c1(2)*zz(:,:,3);
    
    err1 = sum(sum((ytop-z1).^2./abs(z1)));
    
    
    zz = []; ytmp = [];
    for j = 2:numel(tauD)
        tmp = exp(-(xixi/w0).^2./(1+ (tt/tauD(j) + xixi/(v*tauD(j)))))./(1 + tt/tauD(j) + xixi/(v*tauD(j)));
        zz = cat(3,zz,tmp);
        ytmp(:,j-1)=tmp(:)./max(tmp(:));
     end
    zz = cat(3, zz, ones(size(ybottom)));
    zz = zz./repmat(max(max(zz)),size(t,1),size(t,2));
    ytmp(:,end+1)=1;
    c2 = lsqnonneg(ytmp,ybottom(:));
    % c = y(:)\ytmp;
    % weight = repmat(ytop(1,:),size(t,1),1);
    z2 = 0;
    for j = 1:size(zz,3)
        z2 = z2+c2(j)*zz(:,:,j);
    end
    
    % z2 = c2(3)*zz(:,:,1)+c2(1)*zz(:,:,2)+c2(2)*zz(:,:,3);
    % z = c(1)*zz(:,:,1)+c(2)*zz(:,:,2);
    
    % err2 = sum(sum(((ytop-z2)).^2));
    err2 = sum(sum((ybottom-z2).^2./abs(z2)));
    
    err = err1+err2;
    disp([err1 err2 err])
    z = cat(3,z1,z2);
    
    c = [c1;c2];
else
    disp('flag unknown')
    err = NaN;
    c = NaN;
    z = NaN;
    zz = NaN;
end


if 1
    xx = x*pixel;
%     figure
    subplot(231)
    h=surf(xx,t,ytop);
    set(get(h,'Parent'),'YScale','Log');
    title('data')
    shading interp
    view(2)
    axis tight
    box on
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    clims = caxis;
    colorbar
    subplot(232)
    h=surf(xx,t,z2);
    caxis(clims)
    set(get(h,'Parent'),'YScale','Log')
    title('fit')
    shading interp
    view(2)
    axis tight
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    box on
    colorbar
    subplot(233)
    h=surf(xx,t,(ytop-z2)./(ytop));
    set(get(h,'Parent'),'YScale','Log')
    title('residues')
    shading interp
    view(2)
    axis tight
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    box on
    colorbar
    subplot(234)
    h=surf(xx,t,ybottom);
    set(get(h,'Parent'),'YScale','Log');
    title('data')
    shading interp
    view(2)
    axis tight
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    box on
    colorbar
    clims = caxis;
    subplot(235)
    h=surf(xx,t,z1);
    caxis(clims)
    set(get(h,'Parent'),'YScale','Log')
    title('fit')
    shading interp
    view(2)
    axis tight
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    box on
    colorbar
    subplot(236)
    h=surf(xx,t,(ybottom-z1)./(ybottom));
    set(get(h,'Parent'),'YScale','Log')
    title('residues')
    shading interp
    view(2)
    axis tight
    ylabel('lag time / s')
    xlabel('lag space / \mum')
    
    box on
    colorbar
    drawnow
    
end
