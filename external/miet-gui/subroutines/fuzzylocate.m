function [loc idx] = fuzzylocate(num, coord)
% This program gives the location and the index for data points on a list
% of axic values. However, this program will fail if one of the points lies
% equidistant from two axis points. loc is the axis value and idx is the
% index such that coord(idx)=loc.

% (c) Narain Karedla, 2014
if ~isempty(coord)&& ~isempty(num)
    if size(coord,1)>size(coord,2)
        coord=coord.';
    end
    [m,n]=size(coord);
    if min(size(num))>1
        if size(num,1)>size(num,2)
            num=num.';
        end
        p=size(num,2);
        if size(num,1)~=size(coord,1)
            disp('Dimension mismatch')
            loc=[]; idx=[];
            return
        else
            loc=zeros(p,m); idx=loc;
            for i=1:m
                tmpnum=num(i,:);
                dloc=repmat(coord(i,:),[p,1])-repmat(tmpnum.',[1,n]);
                [~,idx(:,i)]=min(abs(dloc.'));
                loc(:,i)=coord(i,idx(:,i));
            end
        end
    else
        idx=zeros(m,1);
        loc=idx;
        num = num(:);
        p = numel(num);
        for i=1:p
            dloc=coord(:)-num(i);
            [~,idx(i)]=min(abs(dloc.'));
            loc(i)=coord(idx(i));
        end
    end
end