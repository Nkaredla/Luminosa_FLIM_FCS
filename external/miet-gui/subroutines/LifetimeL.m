function [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(z,n0,n,n1,d0,d,d1)

% [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(z,n0,n,n1,d0,d,d1) calculates the lifetime and quenching 
% of an oscillating dipole at a distance z from an interface within a layer 
% normalization: free dipole emission rate = 4/3*n
% z  - molecule's distance from the bottom of its layer
% n0 - vector of refractive indices of the stack below the molecule's layer
% n  - refractive index of the molecule's layer
% n1 - vector of refractive indices of the stack above the molecule's layer
% d0 - vector of layer thickness values of the stack below the molecule's layer ( length(d0)=length(n0)-1 )
% d  - thickness of molecule's layer
% d1 - vector of layer thickness values of the stack above the molecule's layer ( length(d1)=length(n1)-1 )

z = z(:)';
col = ones(size(z));
dw = 1e-3;
if all(imag([n0 n1])==0) && max([n0 n1])<n
    w = sqrt(n^2-max([n0 n1])^2) + (0.5:1/dw)'*dw*(n-sqrt(n^2-max([n0 n1])^2));
else
    w = (0.5:1/dw)'*dw*n;
end
dw = w(2)-w(1);
st = sqrt(1-(w/n).^2);
[rpu, rsu, tpu, tsu] = Fresnel(w,[n n1],d1);
[rpd, rsd, tpd, tsd] = Fresnel(w,[n n0(end:-1:1)],d0(end:-1:1));    

% plot(w,abs(rpu),w,abs(rpd));%,w,(angle(rpu)+angle(rpd)+2*w'*d)/2/pi)

rpu = rpu.'; rpd = rpd.'; rsu = rsu.'; rsd = rsd.'; 
tpu = tpu.'; tpd = tpd.'; tsu = tsu.'; tsd = tsd.'; 
ez0 = exp(2*1i*w*z);
ez = exp(2*1i*w*(d-z));
ed = exp(2*1i*w*d);
tp = (1./(1-rpu.*rpd.*ed));
ts = (1./(1-rsu.*rsd.*ed));
v1 = tp*col + ((rpu.*tp)*col).*ez;
v2 = (rpd.*tp)*col + ((rpd.*rpu.*tp)*col).*ez;
pc1 = tp*col - ((rpu.*tp)*col).*ez;
pc2 = (rpd.*tp)*col - ((rpd.*rpu.*tp)*col).*ez;
ps1 = ts*col + ((rsu.*ts)*col).*ez;
ps2 = (rsd.*ts)*col + ((rsu.*rsd.*ts)*col).*ez;

qvd = dw*sum((st.^2*col).*(abs(v1).^2-abs(v2).^2));
qpd = dw*sum((w.^2/n^2*col).*(abs(pc1).^2-abs(pc2).^2) + abs(ps1).^2-abs(ps2).^2)/2;

v1 = tp*col + ((rpd.*tp)*col).*ez0;
v2 = (rpu.*tp)*col + ((rpd.*rpu.*tp)*col).*ez0;
pc1 = tp*col - ((rpd.*tp)*col).*ez0;
pc2 = (rpu.*tp)*col - ((rpd.*rpu.*tp)*col).*ez0;
ps1 = ts*col + ((rsd.*ts)*col).*ez0;
ps2 = (rsu.*ts)*col + ((rsu.*rsd.*ts)*col).*ez0;

qvu = dw*sum((st.^2*col).*(abs(v1).^2-abs(v2).^2));
qpu = dw*sum((w.^2/n^2*col).*(abs(pc1).^2-abs(pc2).^2) + abs(ps1).^2-abs(ps2).^2)/2;

if imag(n0(1))==0
    tp = (tpd./(1-rpu.*rpd.*ed));
    ts = (tsd./(1-rsu.*rsd.*ed));
    tp1 = tp.*(rpu.*ed);
    ts1 = ts.*(rsu.*ed);
    ez = exp(1i*w*z);
    ww = real(sqrt((n0(1)^2-n^2)./w.^2+1))*col;
    
    v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
    pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
    ps = (ts*col).*ez + (ts1*col)./ez;
    
    lvd = dw*sum(ww.*abs(v).^2);
    lpd = dw*sum(ww.*(abs(pc).^2 + abs(ps).^2))/2;
end

if imag(n1(end))==0
    tp = (tpu./(1-rpu.*rpd.*exp(2*1i*w*d)));
    ts = (tsu./(1-rsu.*rsd.*exp(2*1i*w*d)));
    tp1 = tp.*(rpd.*exp(2*1i*w*d));
    ts1 = ts.*(rsd.*exp(2*1i*w*d));
    ez = exp(1i*w*(d-z));
    
    v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
    pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
    ps = (ts*col).*ez + (ts1*col)./ez;
    ww = real(sqrt((n1(end)^2-n^2)./w.^2+1))*col;
    
    lvu = dw*sum(ww.*abs(v).^2);
    lpu = dw*sum(ww.*(abs(pc).^2 + abs(ps).^2))/2;
end

qv = 0;
qp = 0;

if all(imag([n0 n1])==0) && max([n0(2:end) n n1(1:end-1)])>max([n0(1) n1(end)]) % waveguide modes
    dw = 1i*1e-7;
    [ppos, spos] = HashWaveguideMode([n0 n n1],[d0 d d1]);        
    if ~isempty(ppos)
        wp = sqrt(n^2-ppos.^2);
        for j=1:length(ppos)
            ww = wp(j)+(-0.5:0.5)'*dw;
            rpd = Fresnel(ww, [n, n0(end:-1:1)], d0(end:-1:1));        
            rpu = Fresnel(ww, [n, n1], d1);
            rpu = rpu.'*col; rpd = rpd.'*col; 
            fp = diff((1-rpd.*rpu.*exp(2*1i*ww*d*col))./(rpd.*exp(2*1i*ww*z)+rpu.*exp(2*1i*ww*(d-z))+2*rpu.*rpd.*exp(2*1i*d*ww*col)))/dw;
            qv = qv + 4*pi*ppos(j).^2.*imag(1./fp)/n^2;
            qp = qp + 2*pi*imag(abs(wp(j)).^2./fp)/n^2;
        end
    end
    if ~isempty(spos)
        ws = sqrt(n^2-spos.^2);
        for j=1:length(spos)
            ww = ws(j)+(-0.5:0.5)'*dw;
            [~, rsd] = Fresnel(ww, [n, n0(end:-1:1)], d0(end:-1:1));
            [~, rsu] = Fresnel(ww, [n, n1], d1);
            rsu = rsu.'*col; rsd = rsd.'*col;
            fs = diff((1-rsd.*rsu.*exp(2*1i*ww*d*col))./(rsd.*exp(2*1i*ww*z)+rsu.*exp(2*1i*ww*(d-z))+2*rsu.*rsd.*exp(2*1i*d*ww*col)))/dw;
            qp = qp + 2*pi*imag(1./fs);
        end
    end
end

if all(imag([n0 n1])==0) % this part is most shaky!! find a better way for stable integration!!
    wmax = sqrt(1.01*max([n0 n n1]).^2-n^2);
    w = exp(-6 + (6+log(wmax))*(0.5:1e3)'/1e3);
    dw = (6+log(wmax))/1e3;
else
    wmax = sqrt(1e3*max([real(n0) n real(n1)]).^2-n^2);
    w = exp(-6 + (6+log(wmax))*(0.5:1e4)'/1e4);
    dw = (6+log(wmax))/1e4;
    %     wmax = exp(10);
    %     w = exp(-6 + 14*(0.5:1e6)'/1e6);
    %     dw = 14/1e6;
end
if wmax>0
    st = sqrt(1+(w/n).^2);
    [rpu, rsu, tpu, tsu] = Fresnel(1i*w,[n n1],d1);
    [rpd, rsd, tpd, tsd] = Fresnel(1i*w,[n n0(end:-1:1)],d0(end:-1:1));

    if sum(isfinite(rpu))>0
        %         rpu(~isfinite(rpu)) = rpu(sum(isfinite(rpu)));
        %         rsu(~isfinite(rsu)) = rsu(sum(isfinite(rsu)));
        %         tpu(~isfinite(tpu)) = tpu(sum(isfinite(tpu)));
        %         tsu(~isfinite(tsu)) = tsu(sum(isfinite(tsu)));
        %         rpd(~isfinite(rpd)) = rpd(sum(isfinite(rpd)));
        %         rsd(~isfinite(rsd)) = rsd(sum(isfinite(rsd)));
        %         tpd(~isfinite(tpd)) = tpd(sum(isfinite(tpd)));
        %         tsd(~isfinite(tsd)) = tsd(sum(isfinite(tsd)));

        rpu = rpu.'; rpd = rpd.'; rsu = rsu.'; rsd = rsd.';
        tpu = tpu.'; tpd = tpd.'; tsu = tsu.'; tsd = tsd.';
        ed = exp(-2*w*d);

        tp = (1./(1-rpu.*rpd.*ed));
        ts = (1./(1-rsu.*rsd.*ed));

        %         tp(~isfinite(tp)) = tp(sum(isfinite(tp)));
        %         ts(~isfinite(ts)) = ts(sum(isfinite(ts)));

        ez0 = exp(-2*w*z);
        ez = exp(-2*w*(d-z));

        v1 = tp*col + ((rpu.*tp)*col).*ez;
        v2 = (rpd.*tp)*col + ((rpd.*rpu.*tp)*col).*ez;
        pc1 = tp*col - ((rpu.*tp)*col).*ez;
        pc2 = (rpd.*tp)*col - ((rpd.*rpu.*tp)*col).*ez;
        ps1 = ts*col + ((rsu.*ts)*col).*ez;
        ps2 = (rsd.*ts)*col + ((rsu.*rsd.*ts)*col).*ez;

        qvd = qvd + (1./z + 1/2/n^2./z.^3).*imag(conj(v1(end,:)).*v2(end,:));
        qvd = qvd + 2*dw*sum((w.*st.^2*col).*ez0.*(imag(conj(v1).*v2) - ones(numel(w),1)*imag(conj(v1(end,:)).*v2(end,:))));
        qpd = qpd + 1/4/n^2./z.^3.*imag(conj(pc1(end,:)).*pc2(end,:)) + 1/2./z.*imag(conj(ps1(end,:)).*ps2(end,:));
        qpd = qpd + dw*sum(ez0.*((w.^3/n^2*col).*(imag(conj(pc1).*pc2) - ones(numel(w),1)*imag(conj(pc1(end,:)).*pc2(end,:))) + (w*col).*(imag(conj(ps1).*ps2) - ones(numel(w),1)*imag(conj(ps1(end,:)).*ps2(end,:)))));

        flp = ez; ez = ez0; ez0 = flp;
        v1 = tp*col + ((rpd.*tp)*col).*ez;
        v2 = (rpu.*tp)*col + ((rpd.*rpu.*tp)*col).*ez;
        pc1 = tp*col - ((rpd.*tp)*col).*ez;
        pc2 = (rpu.*tp)*col - ((rpd.*rpu.*tp)*col).*ez;
        ps1 = ts*col + ((rsd.*ts)*col).*ez;
        ps2 = (rsu.*ts)*col + ((rsu.*rsd.*ts)*col).*ez;

        qvu = qvu + (1./(d-z) + 1/2/n^2./(d-z).^3).*imag(conj(v1(end,:)).*v2(end,:));
        qvu = qvu + 2*dw*sum((w.*st.^2*col).*ez0.*(imag(conj(v1).*v2) - ones(numel(w),1)*imag(conj(v1(end,:)).*v2(end,:))));
        qpu = qpu + 1/4/n^2./(d-z).^3.*imag(conj(pc1(end,:)).*pc2(end,:)) + 1/2./(d-z).*imag(conj(ps1(end,:)).*ps2(end,:));
        qpu = qpu + dw*sum(ez0.*((w.^3/n^2*col).*(imag(conj(pc1).*pc2) - ones(numel(w),1)*imag(conj(pc1(end,:)).*pc2(end,:))) + (w*col).*(imag(conj(ps1).*ps2) - ones(numel(w),1)*imag(conj(ps1(end,:)).*ps2(end,:)))));
        
        %         qvu = qvu + 2*dw*sum((w.*st.^2*col).*ez0.*imag(conj(v1).*v2));
        %         qpu = qpu + dw*sum(ez0.*((w.^3/n^2*col).*imag(conj(pc1).*pc2) + (w*col).*imag(conj(ps1).*ps2)));

        if imag(n0(1))==0
            tp = (tpd./(1-rpu.*rpd.*ed));
            ts = (tsd./(1-rsu.*rsd.*ed));
            tp1 = tp.*(rpu.*ed);
            ts1 = ts.*(rsu.*ed);
            
            %         tp(~isfinite(tp)) = tp(sum(isfinite(tp)));
            %         ts(~isfinite(ts)) = ts(sum(isfinite(ts)));
            %         tp1(~isfinite(tp1)) = tp1(sum(isfinite(tp1)));
            %         ts1(~isfinite(ts1)) = ts1(sum(isfinite(ts1)));
            
            ez = exp(-w*z);
            ez(w*z>100) = 1e-20;
            v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
            pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
            ps = (ts*col).*ez + (ts1*col)./ez;
            ww = real(sqrt((n0(1)^2-n^2)./w.^2-1))*col;
            
            lvd = lvd + dw*sum((w*col).*ww.*abs(v).^2);
            lpd = lpd + dw*sum((w*col).*ww.*(abs(pc).^2 + abs(ps).^2))/2;
        else
            lvd = qvd;
            lpd = qpd;
        end
        
        if imag(n1(end))==0
            tp = (tpu./(1-rpu.*rpd.*ed));
            ts = (tsu./(1-rsu.*rsd.*ed));
            tp1 = tp.*(rpd.*ed);
            ts1 = ts.*(rsd.*ed);
            
            %             tp(~isfinite(tp)) = tp(sum(isfinite(tp)));
            %             ts(~isfinite(ts)) = ts(sum(isfinite(ts)));
            %             tp1(~isfinite(tp1)) = tp1(sum(isfinite(tp1)));
            %             ts1(~isfinite(ts1)) = ts1(sum(isfinite(ts1)));
            
            ez = exp(-w*(d-z));
            ez(w*(d-z)>100) = 1e-20;
            v = ((st.*tp)*col).*ez + ((st.*tp1)*col)./ez;
            pc = ((w/n.*tp)*col).*ez - ((w/n.*tp1)*col)./ez;
            ps = (ts*col).*ez + (ts1*col)./ez;
            ww = real(sqrt((n1(end)^2-n^2)./w.^2-1))*col;
            
            lvu = lvu + dw*sum((w*col).*ww.*abs(v).^2);
            lpu = lpu + dw*sum((w*col).*ww.*(abs(pc).^2 + abs(ps).^2))/2;
        else
            lvu = qvu;
            lpu = qpu;
        end
           
    else
        qvu = inf;
        qpu = inf;
        qvd = inf;
        qpd = inf;
        lvu = inf;
        lpu = inf;
        lvd = inf;
        lpd = inf;
    end
end

