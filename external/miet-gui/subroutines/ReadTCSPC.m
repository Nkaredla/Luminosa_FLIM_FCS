function [bin, tcspcdata, head] = ReadTCSPC(name)

% [bin, tcspcdata, head] = ReadTCSPC(name)
% name : file name

head      = [];
bin       = [];
tcspcdata = [];

if strcmp(name(end-2:end),'ht3')
    head = HT3_Read(name);
end

if ~isempty(head)
    Resolution = max([0.032 head.Resolution]);
    chDiv      = Resolution/head.Resolution;

    NGate     = ceil(1e9/head.SyncRate/Resolution)+1;
    NCounts   = head.Records;
    resfile   = [name(1:end-4) '.res'];
else
    disp('Not a valid file!');
    return
end


if ( exist(resfile, 'file') == 2) % check if tcspc data was saved
    load(resfile, '-mat');
    if (isfield(res,'bin'))&&(isfield(res,'tcspcdata'))
        bin = res.bin;
        tcspcdata = res.tcspcdata;
    end
end

if isempty(bin)
    photons = 5e6; % number of photons processed one at a time

    dind    = [];
    num     = 1;
    cnt     = 0;
    bin     = 1:NGate;
    
    h = waitbar(0,'Reading ht3-data, please wait...');
    
    while num>0
        [tmpy, tmpx, chan, markers, num] = HT3_Read(name, [cnt+1 photons head.length]);
                       
        % build TCSPC-Hisogramm

        if (num>0)

            cnt = cnt + num;

            tmpx(markers~=0) = [];
            tmpx = round(0.5+tmpx./chDiv);
            chan(markers~=0) = [];
            
            if isempty(dind)
                dind = unique(chan);                
                tcspcdata = zeros(NGate, max(dind)+1);
            end

            for ch = 1:numel(dind)
                tcspcdata(:,dind(ch)+1) = tcspcdata(:,dind(ch)+1) + mHist(tmpx(chan==(dind(ch))), bin);
            end;
        end
        waitbar(cnt/NCounts,h)
        drawnow;
    end

    close(h);
                       
    ind = sum(tcspcdata,2)==0;
    bin(ind) = [];
    tcspcdata(ind,:) = [];
    
    res.bin = bin';
    res.tcspcdata = tcspcdata;
    head.Resolution = Resolution;
    
    save(resfile, 'res', 'head','-v7.3');
end