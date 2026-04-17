function [field,imm,imm_c, xc, yc, mask, al_res, be_res, theta_mol, phi_mol,err]=pattern_detect_auto(pixel,tag,pattern, NA, n0, n, n1, d0, d, d1, lamex, focpos, flag, pic)
% pattern_detect initiates pattern matching program and collects the
% image of identified patterns. The identified image is filtered against
% the background in the image.
%
% Output parameters:
% field is a 3D matrix where each plane contains the identified pixels of
% only one molecule pattern from the intensity image.
% size(field,1)=size(tag,1) size(field,2)=size(tag,2) and size(field,3) =
% number of molecules.
% imm_c is the image containing the patterns of all the molecules
% identified.
% xc is the vector containing the central x-coordinates of all the
% molecule patterns identified.
% yc is the vector containing the central x-coordinates of all the
% molecule patterns identified.
% mask is the matrix containing the patterns calculated by the program.
% al_res is the resolution of the out-of-plane angle used for computing the
% patterns.
% be_res is the resolution of the in-of-plane angle used for computing the
% patterns.
% theta_mol is the vector of the out-of-plane angles of all the molecules.
% phi_mol is the vector of the in-of-plane angles of all the molecules.

err=0;
if nargin<14 || isempty(pic)
    pic=[];
end

if nargin<13 || isempty(flag)
    flag='lifetime';
end

if nargin<12 || isempty(focpos)
    focpos=0;
end

if nargin<3 || isempty(pattern)
    pattern='linear';
end


if strcmpi(flag,'lifetime')|| strcmpi(flag,'MIET') || strcmpi(flag,'pattern match')|| strcmpi(flag,'combined')
    if strcmpi(pattern,'linear')
        al_res=NaN;
        be_res=NaN;
        disp('choose an examplary single molecule for identifying the rest')
        [a,b] = getroi(tag);
        [mx,my,wx,wy] = Gauss2D(tag(a(1):a(2),b(1):b(2)));
        ww = mean([wx,wy]);
        [x,y] = meshgrid(-10:10,-10:10);
        mask = exp(-2*(x.^2+y.^2)/ww^2);
        [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm, imm_c] = FindPatternRapid(tag,mask,mask>0,mask>0,[],0.5,[],0);
        if pic==1
            figure
            mim(cat(3,tag,imm_c));
            
            figure
%             mim(tag(11:end-10,11:end-10));
            mim(tag)
            phi = 0:pi/50:2*pi;
            hold on
            for j=1:length(xc)
%                 plot(xc(j)-10+ww*cos(phi),yc(j)-10+ww*sin(phi),'c');
                 plot(xc(j)+ww*cos(phi),yc(j)+ww*sin(phi),'c');
            end
            hold off
        end
        
        %     mm=round(size(mask,1)./2);
        
    elseif strcmpi(pattern,'azimuthal')
        
        phi= 0:pi/50:pi; % 50 patterns in total, can be changed
        theta=ones(size(phi))*pi/2;
        al_res=NaN;
        be_res=phi(2);
        %                  h.waitbar=waitbar(0,'Detecting Molecules and Identifying patterns....');
        %                  for i=-5:5
        %                     [Pattern(i+6).err,Pattern(i+6).bim, Pattern(i+6).cim, Pattern(i+6).sim, Pattern(i+6).xc, Pattern(i+6).yc, Pattern(i+6).bc,...
        %                      Pattern(i+6).cc, Pattern(i+6).sc, Pattern(i+6).len, Pattern(i+6).imm, Pattern(i+6).imm_c, Pattern(i+6).mask, Pattern(i+6).molim]...
        %                      = radialpatterns_focus(tag, pattern, pixel+0.001*i, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);
        %                     imtest{i+6,1} = Pattern(i+6).molim;
        %                     h.waitbar=waitbar((i+6)/11);
        %                     drawnow
        %                  end
        %                  delete(h.waitbar);
        
        [err, ~, ~, ~, xc, yc, ~, ~, sc, ~, imm, imm_c, mask, ~] = radialpatterns_pixel(tag, pattern, pixel, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);
    elseif strcmpi(pattern,'radial')
        be_res = 5; % minimum resolution of in-plane angle
        al_res = 5; % minimum resolution of out-of-plane angle
        cnt = 1;
        for k=90:-al_res:0
            al = k/180*pi;
            if k==90
                jj = round(180/be_res);
                dbe = pi/jj;
            elseif k==0
                jj = 1;
                dbe = 0;
            else
                jj = round(sin(al)*360/be_res);
                dbe = 2*pi/jj;
            end
            for j=1:jj
                theta(cnt) = al;
                phi(cnt) = dbe*(j-1);
                cnt = cnt+1;
            end
        end
        %                 h.waitbar=waitbar(0,'Detecting Molecules and Identifying patterns....');
        %                 for i=-5:5
        %                      [Pattern(i+6).err,Pattern(i+6).bim, Pattern(i+6).cim, Pattern(i+6).sim, Pattern(i+6).xc, Pattern(i+6).yc, Pattern(i+6).bc,...
        %                      Pattern(i+6).cc, Pattern(i+6).sc, Pattern(i+6).len, Pattern(i+6).imm, Pattern(i+6).imm_c, Pattern(i+6).mask, Pattern(i+6).molim]...
        %                      = radialpatterns_focus(tag, pattern, pixel+0.001*i, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);
        %                      imtest{i+6,1} = Pattern(i+6).molim;
        %                      h.waitbar=waitbar((i+6)/11);
        %                      drawnow;
        %                 end
        %                 delete(h.waitbar);
        [err, ~, ~, ~, xc, yc, ~, ~, sc, ~, imm, imm_c, mask, ~] = radialpatterns_pixel(tag, pattern, pixel, theta, phi, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);
    end
    close all
%     if ~flagauto
%         if strcmpi(pattern,'azimuthal')|| strcmpi(pattern,'radial')
%             disp('Click leftarrow/right arrow to scroll images and find the best fitted patterns for the molecules identified.')
%             disp('Click Enter while viewing the best fitted figure')
%             
%             pind  = SelectImage(imtest);
%             close all
%             imm   = Pattern(pind).imm;
%             imm_c = Pattern(pind).imm_c;
%             xc    = Pattern(pind).xc;
%             yc    = Pattern(pind).yc;
%             sc    = Pattern(pind).sc;
%             mask  = Pattern(pind).mask;
%         end
%     else
%         nmol=zeros(numel(Pattern,1));
%         for i=1:numel(Pattern)
%             nmol(i)=numel(Pattern(i).len);
%         end
%         [~,pind]=max(nmol);
%         imm   = Pattern(pind).imm;
%         imm_c = Pattern(pind).imm_c;
%         xc    = Pattern(pind).xc;
%         yc    = Pattern(pind).yc;
%         sc    = Pattern(pind).sc;
%         mask  = Pattern(pind).mask;
%     end
%     
    
elseif strcmpi(flag,'random')
    al_res=NaN;
    be_res=NaN;
    [err, ~, ~, ~, xc, yc, ~, ~, sc, ~, imm, imm_c, mask] = radialpatterns_pixel(tag, pattern, pixel, NA, n0, n, n1, d0, d, d1, lamex, focpos, pic);
end

bcgind=imm_c==0;
if size(imm_c)==size(tag)
    bcgimm=bcgind.*tag;
else
    err=1;
    field=[];
    imm=[];imm_c=[]; xc=[]; yc=[]; mask=[]; al_res=[]; be_res=[]; theta_mol=[]; phi_mol=[];
    return
end
bcgval=bcgimm(:);
bcgind=bcgval==0;
bcgval(bcgind)=[];
% bcgind=bcgval>=mean(bcgval)+2.*std(bcgval);
% bcgval(bcgind)=[];
% bg=mean(bcgval);
bg=mode(bcgval);
posimm=zeros(size(imm_c));
for i=1:numel(xc)
    posimm(yc(i),xc(i))=1;
end

field  = zeros(size(imm,1), size(imm,2),numel(sc));
bright = zeros(numel(sc),1);
if numel(xc)>0
    for i=1:numel(sc)
        ind = imm(:,:,i)>max(max(imm(:,:,i))).*0.04;
        field(:,:,i)=ind.*tag.*(mConv2(posimm,Disk(ceil(lamex/2./pixel)))>1e-5);
        %              bright(i)=max(max(field(:,:,i)))./bg;
        pix(i)=sum(sum(ind.*(mConv2(posimm,Disk(ceil(lamex/2./pixel)))>1e-5)));
        bright(i)=sum(sum(field(:,:,i)))./(pix(i).*bg);
        if bright(i)<1.5 || bright(i)==1.5
            field(:,:,i)=0;
        else
            ind=field(:,:,i)>0.05.*(max(max(field(:,:,i)))); % 5%
            field(:,:,i)=ind.*1;
        end
        empt(i)=sum(sum(field(:,:,i)));
        sum_phot(i)=sum(sum(field(:,:,i).*tag));
    end
    
    ind=sum_phot<100; % cutoff 500 photons
    field(:,:,ind)=[];
    empt(ind)=[];
    sc(ind)=[];
    xc(ind)=[];
    yc(ind)=[];
    
    ind=empt<=ceil(0.1.*max(empt));
    field(:,:,ind)=[];
    sc(ind)=[];
    xc(ind)=[];
    yc(ind)=[];
    if strcmpi(flag,'lifetime')|| strcmpi(flag,'MIET') || strcmpi(flag,'pattern match')|| strcmpi(flag,'combined')
        if strcmpi(pattern,'radial')  % || strcmpi(pattern,'azimuthal')
            theta     = repmat(theta,[1 6]);
            phi       = repmat(phi, [1 6]);
            theta_mol = theta(sc);
            phi_mol   = phi(sc);
        elseif strcmpi(pattern,'azimuthal')
            phi        = repmat(phi, [1 6]);
            theta_mol  = NaN;
            phi_mol    = phi(sc);
        else
            theta_mol = NaN;
            phi_mol   = NaN;
        end
    elseif strcmpi(flag,'random')
        theta_mol = NaN;
        phi_mol   = NaN;
    end
    
    
    disp(['Molecules detected...',num2str(numel(xc))]);
    
else
    disp('No single molecules found!! Check model patterns')
end
end

