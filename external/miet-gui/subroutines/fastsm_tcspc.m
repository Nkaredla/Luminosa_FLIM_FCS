function[field, fieldcell, tcspc_sm, tcspcirf]=fastsm_tcspc(bcgind, field, Ngate, im_line, im_col, im_tcspc, im_chan)
% fastsm_tcspc.m collects the photons from the pixels identified for each
% single molecule separately.

% (c) Narain Karedla, 2014
disp('Collecting photons from identified molecules...')

    maxch=size(field,4);
    tcspc_sm=zeros(size(field,3), Ngate, maxch);
    fieldcell= cell(size(field,1),size(field,2),maxch);
    z=size(field,3);

    if ~isempty(bcgind)
        tcspcirf=zeros(Ngate,maxch);
        disp('Calculating IRF...')
    else
        tcspcirf=[];
    end

    for ch=1:size(maxch)
        for p=1:z
            field(:,:,p,ch)=p.*field(:,:,p,ch);
        end
        for x=1:size(field,1)
            for y=1:size(field,2)
                fieldcell{x,y,ch}=unique(field(x,y,:,ch));
                tcspc_pix = zeros(Ngate,1);
                if isempty(bcgind);
                    if sum(ismember(fieldcell{x,y,ch},1:z))>0
                       zind = fieldcell{x,y,ch}(:);
                       ind = zind==0;
                       zind(ind) = [];
                       ind1 = (im_line == y);
                       ind2 = (im_col == x);
                       tmp1 = double(im_tcspc).*(ind1.*ind2);
                       tmp2 = double(im_chan).*(ind1.*ind2);
                       tcspc_pix=mHist(double(tmp1(tmp2 == ch)),1:Ngate);
                       for i=1:numel(zind)
                           tcspc_sm(zind(i),:,ch)=tcspc_sm(zind(i),:,ch)+tcspc_pix';
                       end
                    end
                else
                    if sum(ismember(fieldcell{x,y,ch},1:z))>0
                       zind = fieldcell{x,y,ch}(:);
                       ind = zind==0;
                       zind(ind) = [];
                       ind1 = (im_line == y);
                       ind2 = (im_col == x);
                       tmp1 = double(im_tcspc).*(ind1.*ind2);
                       tmp2 = double(im_chan).*(ind1.*ind2);
                       tcspc_pix=mHist(double(tmp1(tmp2 == ch)),1:Ngate);
                       for i=1:numel(zind)
                           tcspc_sm(zind(i),:,ch)=tcspc_sm(zind(i),:,ch)+tcspc_pix';
                       end
                    elseif bcgind(y,x)==1
                        ind1 = (im_line == y);
                        ind2 = (im_col == x);
                        tmp1 = double(im_tcspc).*(ind1.*ind2);
                        tmp2 = double(im_chan).*(ind1.*ind2);
                        tcspcirf=tcspcirf+mHist(double(tmp1(tmp2 == ch)),1:Ngate);
                    end

                end
            end
        end
    end

    
end
               

