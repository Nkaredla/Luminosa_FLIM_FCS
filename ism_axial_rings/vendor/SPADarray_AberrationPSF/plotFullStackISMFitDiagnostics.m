function files = plotFullStackISMFitDiagnostics(resultInput,varargin)
%PLOTFULLSTACKISMFITDIAGNOSTICS Visualize measured and fitted ISM intensities.
%
%   plotFullStackISMFitDiagnostics(result)
%   plotFullStackISMFitDiagnostics('full_stack_ism_wavefront_fit.mat')
%
%   Detector values are displayed on the regular 23-pixel honeycomb. The
%   measured/calibrated detector coordinates remain in the forward model and
%   are not replaced by this display geometry.

    p=inputParser;
    addParameter(p,'outputDir','');
    addParameter(p,'visible','off');
    addParameter(p,'planesPerPage',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(p,varargin{:});
    opts=p.Results;

    result=loadResult(resultInput);
    if isempty(opts.outputDir)
        if isfield(result,'outputDir') && ~isempty(result.outputDir)
            outDir=result.outputDir;
        else
            outDir=pwd;
        end
    else
        outDir=char(opts.outputDir);
    end
    if exist(outDir,'dir')~=7, mkdir(outDir); end

    [metrics,displayXY]=detectorMetrics(result);
    csvFile=fullfile(outDir,'selected_plane_detector_fit.csv');
    writetable(metrics.table,csvFile);

    quality=planeFitQuality(result,metrics);
    qualityFile=fullfile(outDir,'selected_plane_fit_quality.csv');
    writetable(quality.table,qualityFile);

    overviewPositions=selectOverviewPlanes(result,metrics);
    hexFile=fullfile(outDir,'selected_plane_detector_hex_fit.png');
    writeDetectorHexFigure(result,metrics,quality,displayXY, ...
        overviewPositions,hexFile,opts.visible,'overview');

    spatialFile=fullfile(outDir,'selected_plane_spatial_fit.png');
    writeSpatialFigure(result,metrics,quality,overviewPositions, ...
        spatialFile,opts.visible,'overview');

    pageDir=fullfile(outDir,'selected_plane_fit_pages');
    if exist(pageDir,'dir')~=7, mkdir(pageDir); end
    [hexPages,spatialPages]=writeSelectedPlanePages(result,metrics,quality, ...
        displayXY,pageDir,opts.visible,round(opts.planesPerPage));

    files=struct('detectorTable',csvFile,'qualityTable',qualityFile, ...
        'detectorHexFigure',hexFile,'spatialFigure',spatialFile, ...
        'detectorHexPages',{hexPages},'spatialPages',{spatialPages});
end

function result=loadResult(inputValue)
    if isstruct(inputValue)
        result=inputValue;
        return;
    end
    S=load(char(inputValue));
    if isfield(S,'result') && isstruct(S.result)
        result=S.result;
        return;
    end
    names=fieldnames(S);
    for k=1:numel(names)
        if isstruct(S.(names{k})) && isfield(S.(names{k}),'fit')
            result=S.(names{k});
            return;
        end
    end
    error('plotFullStackISMFitDiagnostics:MissingResult', ...
        'No full-stack result structure was found.');
end

function [metrics,displayXY]=detectorMetrics(result)
    data=result.data;
    fit=result.fit;
    planeIdx=result.selectedPlaneIndices(:).';
    nPlane=numel(planeIdx);
    nCh=size(data.rawCounts,3);

    if nCh==23
        [displayXY,~,~]=detectorLayout('honeycomb23',1);
    else
        theta=linspace(0,2*pi,nCh+1).';
        displayXY=[cos(theta(1:end-1)) sin(theta(1:end-1))];
    end

    measuredSignal=zeros(nCh,nPlane);
    fittedSignal=zeros(nCh,nPlane);
    pearson=zeros(nCh,nPlane);
    rawTotal=zeros(nCh,nPlane);
    backgroundTotal=zeros(nCh,nPlane);

    nPixel=size(data.rawCounts,1)*size(data.rawCounts,2);
    for ip=1:nPlane
        iz=planeIdx(ip);
        raw=double(data.rawCounts(:,:,:,iz));
        background=repmat(double(data.backgroundPerPixel(:,:,:,iz)), ...
            size(raw,1),size(raw,2),1);
        fitted=fit.globalPhotonScale*double(fit.model(:,:,:,iz));
        expected=max(fitted+background,realmin);

        rawTotal(:,ip)=squeeze(sum(sum(raw,1),2));
        backgroundTotal(:,ip)=squeeze(sum(sum(background,1),2));
        measuredSignal(:,ip)=squeeze(sum(sum(max(raw-background,0),1),2));
        fittedSignal(:,ip)=squeeze(sum(sum(fitted,1),2));
        pearson(:,ip)=squeeze(sum(sum(raw-expected,1),2))./ ...
            sqrt(max(squeeze(sum(sum(expected,1),2)),1));
    end

    detectorIndex=repmat((1:nCh).',nPlane,1);
    channelID=repmat(data.channelIDs(:),nPlane,1);
    planeIndex=kron(planeIdx(:),ones(nCh,1));
    stageZUm=kron(data.stageZUm(planeIdx(:)).',ones(nCh,1));
    T=table(planeIndex,stageZUm,detectorIndex,channelID, ...
        rawTotal(:),backgroundTotal(:),measuredSignal(:), ...
        fittedSignal(:),(measuredSignal(:)-fittedSignal(:)),pearson(:), ...
        'VariableNames',{'planeIndex','stageZUm','detectorIndex','channelID', ...
        'rawCounts','backgroundCounts','measuredSignalCounts', ...
        'fittedSignalCounts','signalResidualCounts','pearsonResidual'});

    metrics=struct('planeIdx',planeIdx,'measuredSignal',measuredSignal, ...
        'fittedSignal',fittedSignal,'pearson',pearson,'table',T, ...
        'nPixel',nPixel);
end

function quality=planeFitQuality(result,M)
    nPlane=numel(M.planeIdx);
    detectorCorrelation=nan(nPlane,1);
    detectorNRMSE=nan(nPlane,1);
    detectorPearsonRMS=nan(nPlane,1);
    spatialCorrelation=nan(nPlane,1);
    spatialNRMSE=nan(nPlane,1);
    spatialPearsonRMS=nan(nPlane,1);

    for ip=1:nPlane
        detectorCorrelation(ip)=vectorCorrelation( ...
            M.measuredSignal(:,ip),M.fittedSignal(:,ip));
        detectorNRMSE(ip)=relativeRmse( ...
            M.measuredSignal(:,ip),M.fittedSignal(:,ip));
        detectorPearsonRMS(ip)=finiteRms(M.pearson(:,ip));

        [measured,fitted,pearson]=spatialPlaneData(result,M.planeIdx(ip));
        spatialCorrelation(ip)=vectorCorrelation(measured,fitted);
        spatialNRMSE(ip)=relativeRmse(measured,fitted);
        spatialPearsonRMS(ip)=finiteRms(pearson);
    end

    planeIndex=M.planeIdx(:);
    stageZUm=reshape(result.data.stageZUm(planeIndex),[],1);
    T=table(planeIndex,stageZUm,detectorCorrelation,detectorNRMSE, ...
        detectorPearsonRMS,spatialCorrelation,spatialNRMSE, ...
        spatialPearsonRMS);
    quality=struct('detectorCorrelation',detectorCorrelation, ...
        'detectorNRMSE',detectorNRMSE, ...
        'detectorPearsonRMS',detectorPearsonRMS, ...
        'spatialCorrelation',spatialCorrelation, ...
        'spatialNRMSE',spatialNRMSE, ...
        'spatialPearsonRMS',spatialPearsonRMS,'table',T);
end

function positions=selectOverviewPlanes(result,M)
    z=result.data.stageZUm(M.planeIdx);
    if isfield(result,'fit') && isfield(result.fit,'estZ0Um') && ...
            isfinite(result.fit.estZ0Um)
        z0=result.fit.estZ0Um;
    else
        [~,centerPosition]=max(sum(M.measuredSignal,1));
        z0=z(centerPosition);
    end
    targets=z0+[-0.5 0 0.5];
    positions=zeros(1,numel(targets));
    available=1:numel(z);
    for k=1:numel(targets)
        [~,j]=min(abs(z(available)-targets(k)));
        positions(k)=available(j);
        available(j)=[];
        if isempty(available), break; end
    end
    positions=sort(unique(positions(positions>0)));
end

function [hexPages,spatialPages]=writeSelectedPlanePages( ...
        result,M,Q,displayXY,pageDir,visible,planesPerPage)
    nPlane=numel(M.planeIdx);
    nPage=ceil(nPlane/planesPerPage);
    hexPages=cell(nPage,1);
    spatialPages=cell(nPage,1);
    for page=1:nPage
        first=(page-1)*planesPerPage+1;
        positions=first:min(first+planesPerPage-1,nPlane);
        pageLabel=sprintf('page %d of %d',page,nPage);
        hexPages{page}=fullfile(pageDir,sprintf( ...
            'detector_hex_fit_page_%02d.png',page));
        spatialPages{page}=fullfile(pageDir,sprintf( ...
            'spatial_fit_page_%02d.png',page));
        writeDetectorHexFigure(result,M,Q,displayXY,positions, ...
            hexPages{page},visible,pageLabel);
        writeSpatialFigure(result,M,Q,positions,spatialPages{page}, ...
            visible,pageLabel);
    end
end

function writeDetectorHexFigure(result,M,Q,displayXY,positions, ...
        outFile,visible,pageLabel)
    nPlane=numel(positions);
    fig=figure('Visible',visible,'Color','w', ...
        'Position',[40 40 1500 max(520,360*nPlane)]);
    tl=tiledlayout(fig,nPlane,3,'Padding','compact','TileSpacing','compact');
    for row=1:nPlane
        ip=positions(row);
        measured=M.measuredSignal(:,ip);
        fitted=M.fittedSignal(:,ip);
        residual=M.pearson(:,ip);
        signalLimit=max([measured;fitted],[],'omitnan');
        if ~isfinite(signalLimit) || signalLimit<=0, signalLimit=1; end
        residualLimit=max(abs(residual),[],'omitnan');
        if ~isfinite(residualLimit) || residualLimit<=0, residualLimit=1; end

        ax=nexttile(tl);
        plotDetectorHexMap(displayXY,measured,'Parent',ax,'CLim',[0 signalLimit]);
        colormap(ax,'turbo');
        colorbar(ax);
        title(ax,sprintf('Measured signal, z=%.3f um', ...
            result.data.stageZUm(M.planeIdx(ip))));

        ax=nexttile(tl);
        plotDetectorHexMap(displayXY,fitted,'Parent',ax,'CLim',[0 signalLimit]);
        colormap(ax,'turbo');
        colorbar(ax);
        title(ax,sprintf('Fitted: r=%.3f, NRMSE=%.3f', ...
            Q.detectorCorrelation(ip),Q.detectorNRMSE(ip)));

        ax=nexttile(tl);
        plotDetectorHexMap(displayXY,residual,'Parent',ax, ...
            'CLim',[-residualLimit residualLimit]);
        colormap(ax,redBlueMap(256));
        colorbar(ax);
        title(ax,sprintf('Pearson residual: RMS=%.2f', ...
            Q.detectorPearsonRMS(ip)));
    end
    status=fitStatus(result);
    title(tl,sprintf('%s: detector-integrated fit (%s)', ...
        status.text,pageLabel), ...
        'Color',status.color,'FontWeight','bold');
    exportFigure(fig,outFile);
end

function writeSpatialFigure(result,M,Q,positions,outFile,visible,pageLabel)
    nPlane=numel(positions);
    fig=figure('Visible',visible,'Color','w', ...
        'Position',[40 40 1500 max(520,360*nPlane)]);
    tl=tiledlayout(fig,nPlane,3,'Padding','compact','TileSpacing','compact');
    for row=1:nPlane
        ip=positions(row);
        iz=M.planeIdx(ip);
        [measured,fitted,pearson]=spatialPlaneData(result,iz);

        signalLimit=max([measured(:);fitted(:)],[],'omitnan');
        if ~isfinite(signalLimit) || signalLimit<=0, signalLimit=1; end
        residualLimit=robustAbsLimit(pearson,0.995);
        if ~isfinite(residualLimit) || residualLimit<=0, residualLimit=1; end

        ax=nexttile(tl);
        imagesc(ax,measured,[0 signalLimit]);
        formatImageAxes(ax);
        colormap(ax,'hot');
        colorbar(ax); title(ax,sprintf('Measured XY, z=%.3f um', ...
            result.data.stageZUm(iz)));

        ax=nexttile(tl);
        imagesc(ax,fitted,[0 signalLimit]);
        formatImageAxes(ax);
        colormap(ax,'hot');
        colorbar(ax); title(ax,sprintf('Fitted: r=%.3f, NRMSE=%.3f', ...
            Q.spatialCorrelation(ip),Q.spatialNRMSE(ip)));

        ax=nexttile(tl);
        imagesc(ax,pearson,[-residualLimit residualLimit]);
        formatImageAxes(ax);
        colormap(ax,redBlueMap(256));
        colorbar(ax); title(ax,sprintf('Pearson residual: RMS=%.2f', ...
            Q.spatialPearsonRMS(ip)));
    end
    status=fitStatus(result);
    title(tl,sprintf('%s: detector-summed spatial fit (%s)', ...
        status.text,pageLabel), ...
        'Color',status.color,'FontWeight','bold');
    exportFigure(fig,outFile);
end

function [measured,fitted,pearson]=spatialPlaneData(result,iz)
    raw=sum(double(result.data.rawCounts(:,:,:,iz)),3);
    background=sum(repmat(double( ...
        result.data.backgroundPerPixel(:,:,:,iz)), ...
        size(raw,1),size(raw,2),1),3);
    measured=max(raw-background,0);
    fitted=result.fit.globalPhotonScale* ...
        sum(double(result.fit.model(:,:,:,iz)),3);
    expected=max(fitted+background,realmin);
    pearson=(raw-expected)./sqrt(expected);
end

function formatImageAxes(ax)
    axis(ax,'image');
    axis(ax,'tight');
    pbaspect(ax,[1 1 1]);
    xlabel(ax,'x pixel');
    ylabel(ax,'y pixel');
end

function value=vectorCorrelation(a,b)
    a=double(a(:));
    b=double(b(:));
    valid=isfinite(a)&isfinite(b);
    a=a(valid);
    b=b(valid);
    if numel(a)<2
        value=NaN;
        return;
    end
    a=a-mean(a);
    b=b-mean(b);
    denominator=norm(a)*norm(b);
    if denominator<=0
        value=NaN;
    else
        value=(a.'*b)/denominator;
    end
end

function value=relativeRmse(a,b)
    a=double(a(:));
    b=double(b(:));
    valid=isfinite(a)&isfinite(b);
    a=a(valid);
    b=b(valid);
    if isempty(a)
        value=NaN;
    else
        value=norm(a-b)/max(norm(a),eps);
    end
end

function value=finiteRms(x)
    x=double(x(:));
    x=x(isfinite(x));
    if isempty(x)
        value=NaN;
    else
        value=sqrt(mean(x.^2));
    end
end

function limit=robustAbsLimit(x,quantileValue)
    values=sort(abs(double(x(isfinite(x)))));
    if isempty(values)
        limit=NaN;
        return;
    end
    index=max(1,min(numel(values),ceil(quantileValue*numel(values))));
    limit=values(index);
end

function status=fitStatus(result)
    accepted=false;
    if isfield(result,'acceptance') && isfield(result.acceptance,'allAccepted')
        accepted=logical(result.acceptance.allAccepted);
    end
    if accepted
        status=struct('text','FIT ACCEPTED','color',[0 0.45 0]);
    else
        status=struct('text','FIT REJECTED','color',[0.75 0 0]);
    end
end

function map=redBlueMap(n)
    x=linspace(-1,1,n).';
    map=[min(1,1+x),1-abs(x),min(1,1-x)];
    map=max(0,min(1,map));
end

function exportFigure(fig,fileName)
    try
        exportgraphics(fig,fileName,'Resolution',180);
    catch
        set(fig,'PaperPositionMode','auto');
        print(fig,fileName,'-dpng','-r180');
    end
    close(fig);
end
