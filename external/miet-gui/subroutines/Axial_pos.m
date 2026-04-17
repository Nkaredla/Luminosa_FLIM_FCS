function [ax_pos,err] = Axial_pos(pol_ang, lifetime, z, lifecurve)
% This program calculates the axial location of single molecules, given
% their polar angle, quantum yield and lifetime.

% (c) Narain Karedla, 2018
    z=z(:); pol_ang=pol_ang(:);
    pol_ang=round(pol_ang.*180/pi);
    ind = isnan(pol_ang);
    pol_ang(ind)=90;
    tmp=size(lifecurve,2);
    theta= (90:-(90/(tmp-1)):0);
    if size(lifecurve,1)<size(lifecurve,2)
        lifecurve=lifecurve.';
    end
    d_life=zeros(size(lifecurve));
    for i=1:size(lifecurve,2)
        d_life(:,i)=gradient(lifecurve(:,i));
    end
    ax_pos=zeros(numel(lifetime),1);
    err=zeros(numel(lifetime),1);
    for i=1:numel(lifetime)
        hdif=lifetime(:,i)-lifecurve(:,pol_ang(i,:));
        if min(abs(hdif))<max(d_life(:,pol_ang(i)==theta))|| min(abs(hdif))==max(d_life(:,pol_ang(i)==theta))
            [~,pos]=min(abs(hdif));
            ax_pos(i)=z(pos);
            err(i)=hdif(pos)*d_life(pos,pol_ang(i)==theta);
        else
            ax_pos(i)=NaN; err=NaN;
        end
        clear pos; clear hdif
    end
end
