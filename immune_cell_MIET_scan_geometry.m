function info = immune_cell_MIET_scan_geometry(headOrFile)
%IMMUNE_CELL_MIET_SCAN_GEOMETRY Conservatively classify a Luminosa PTU scan plane.
%   INFO = IMMUNE_CELL_MIET_SCAN_GEOMETRY(HEAD) examines a header returned by
%   PTU_Read_Head, or the .head member of a PTU reader result.
%
%   INFO = IMMUNE_CELL_MIET_SCAN_GEOMETRY(PTUFILE) reads the PTU header first.
%
%   INFO contains these stable fields:
%       plane           - 'XY', 'XZ', 'YZ', 'XZ/YZ', or 'unknown'
%       confidence      - 'high', 'medium', or 'low'
%       reasons         - cell array of explanatory character vectors
%       relevantFields  - raw header values used for audit/reporting
%
%   Luminosa writes ImgHdr_ScanDirection as a numeric code.  PTU_Read_Head
%   exposes the code but does not define its enumeration.  For files created
%   by Luminosa, this helper uses the dataset-derived project convention
%   0 = XY, 1 = YZ, and 2 = XZ.  This is not claimed to be an official PTU
%   enumeration, so a result based only on that code has medium confidence.
%
%   ImgHdr_Dimensions is retained for reporting but is not used as a plane
%   label.  In observed Luminosa 2-D PTU files it is 3 for XY and orthogonal
%   scans alike, so interpreting it as an XYZ volume would be unsafe.

    narginchk(1, 1);
    [head, sourceFile] = resolveHeader(headOrFile);

    info = struct();
    info.plane = 'unknown';
    info.confidence = 'low';
    info.reasons = {};
    info.relevantFields = collectRelevantFields(head);
    info.sourceFile = sourceFile;
    info.convention = [ ...
        'Dataset-derived Luminosa convention: ' ...
        'ImgHdr_ScanDirection 0=XY, 1=YZ, 2=XZ (not an official PTU enum).'];

    [textPlane, textEvidence, textConflict] = planeFromTextMetadata(head);
    [directionPlane, directionReason, rawDirection] = ...
        planeFromLuminosaDirection(head);
    info.rawScanDirection = rawDirection;

    if textConflict
        info.reasons = {textEvidence, ...
            'Conflicting textual plane labels prevent a safe classification.'};
        return;
    end

    if ~strcmp(textPlane, 'unknown')
        if planesConflict(textPlane, directionPlane)
            info.reasons = {textEvidence, directionReason, ...
                'Text metadata and the Luminosa direction category disagree.'};
            return;
        end

        info.plane = textPlane;
        info.confidence = 'high';
        info.reasons = {textEvidence};
        if ~strcmp(directionPlane, 'unknown')
            info.reasons = [info.reasons, {directionReason}];
        end
        return;
    end

    if ~strcmp(directionPlane, 'unknown')
        info.plane = directionPlane;
        info.confidence = 'medium';
        info.reasons = {directionReason};
        return;
    end

    if ~isempty(directionReason)
        info.reasons = {directionReason};
    else
        info.reasons = { ...
            'No usable scan-plane or scan-direction metadata was found.'};
    end

    if isfield(head, 'ImgHdr_Dimensions')
        value = scalarFiniteDouble(head.ImgHdr_Dimensions);
        if ~isempty(value)
            info.reasons = [info.reasons, {sprintf( ...
                ['ImgHdr_Dimensions=%g was not treated as a plane label; ' ...
                 'Luminosa writes 3 for multiple kinds of 2-D scan.'], value)}];
        end
    end
end

function [head, sourceFile] = resolveHeader(headOrFile)
    sourceFile = '';

    if ischar(headOrFile) || ...
            (isstring(headOrFile) && isscalar(headOrFile))
        sourceFile = char(headOrFile);
        if exist(sourceFile, 'file') ~= 2
            error('immune_cell_MIET_scan_geometry:FileNotFound', ...
                'PTU file does not exist: %s', sourceFile);
        end
        if exist('PTU_Read_Head', 'file') ~= 2
            error('immune_cell_MIET_scan_geometry:MissingReader', ...
                ['PTU_Read_Head.m must be on the MATLAB path when the input ' ...
                 'is a filename.']);
        end
        head = PTU_Read_Head(sourceFile);
    elseif isstruct(headOrFile) && isscalar(headOrFile)
        if isfield(headOrFile, 'head') && isstruct(headOrFile.head) && ...
                isscalar(headOrFile.head)
            head = headOrFile.head;
        else
            head = headOrFile;
        end
    else
        error('immune_cell_MIET_scan_geometry:InvalidInput', ...
            'Input must be a scalar header/PTU struct or a PTU filename.');
    end

    if ~isstruct(head) || ~isscalar(head) || isempty(fieldnames(head))
        error('immune_cell_MIET_scan_geometry:EmptyHeader', ...
            'The PTU header is empty or invalid.');
    end
end

function fields = collectRelevantFields(head)
    names = { ...
        'CreatorSW_Name', 'CreatorSW_Version', 'File_Comment', ...
        'Measurement_Mode', 'Measurement_SubMode', ...
        'ImgHdr_Ident', 'ImgHdr_Dimensions', 'ImgHdr_ScanDirection', ...
        'ImgHdr_XCalibration', 'ImgHdr_YCalibration', ...
        'ImgHdr_ZCalibration', 'ImgHdr_XCalibOffset', ...
        'ImgHdr_YCalibOffset', 'ImgHdr_ZCalibOffset', ...
        'ImgHdr_X0', 'ImgHdr_Y0', 'ImgHdr_Z0', ...
        'ImgHdr_PixX', 'ImgHdr_PixY', 'ImgHdr_PixResol', ...
        'ImgHdr_MaxFrames', 'ImgHdr_ObjectiveName'};

    fields = struct();
    for k = 1:numel(names)
        name = names{k};
        if isfield(head, name)
            fields.(name) = head.(name);
        end
    end

    % Preserve non-standard plane/axis tags for auditability.
    allNames = fieldnames(head);
    for k = 1:numel(allNames)
        name = allNames{k};
        if ~isempty(regexpi(name, '(plane|scanaxis|scan_axis|axisname)', 'once'))
            fields.(name) = head.(name);
        end
    end
end

function [plane, evidence, conflict] = planeFromTextMetadata(head)
    plane = 'unknown';
    evidence = '';
    conflict = false;
    hits = {};
    hitFields = {};

    names = fieldnames(head);
    for k = 1:numel(names)
        name = names{k};
        if ~isTextMetadataField(name)
            continue;
        end
        text = textValue(head.(name));
        if isempty(text)
            continue;
        end

        planes = planesNamedInText(text);
        hits = [hits, planes]; %#ok<AGROW>
        hitFields = [hitFields, repmat({name}, 1, numel(planes))]; %#ok<AGROW>
    end

    if isempty(hits)
        return;
    end

    uniqueHits = stableUnique(hits);
    if numel(uniqueHits) > 1
        conflict = true;
        evidence = sprintf('Header text contains multiple plane labels (%s).', ...
            strjoin(uniqueHits, ', '));
        return;
    end

    plane = uniqueHits{1};
    evidence = sprintf('Header field %s explicitly labels the scan plane as %s.', ...
        hitFields{1}, plane);
end

function tf = isTextMetadataField(name)
    fixed = {'File_Comment', 'ImgHdr_ScanDirection', 'ImgHdr_ScanPlane', ...
        'ImgHdr_Plane', 'ImgHdr_Axis', 'ImgHdr_ScanAxis', ...
        'ScanDirection', 'ScanPlane', 'ScanAxis'};
    tf = any(strcmpi(name, fixed)) || ...
        ~isempty(regexpi(name, '(plane|scan.?axis)', 'once'));
end

function text = textValue(value)
    text = '';
    if ischar(value)
        text = value(:).';
    elseif isstring(value) && isscalar(value) && ~ismissing(value)
        text = char(value);
    elseif iscell(value) && isscalar(value)
        text = textValue(value{1});
    end
end

function planes = planesNamedInText(text)
    planes = {};
    text = lower(text);
    if ~isempty(regexp(text, '(^|[^a-z0-9])x[\s_\-/]*y([^a-z0-9]|$)', 'once'))
        planes = [planes, {'XY'}];
    end
    if ~isempty(regexp(text, '(^|[^a-z0-9])x[\s_\-/]*z([^a-z0-9]|$)', 'once'))
        planes = [planes, {'XZ'}];
    end
    if ~isempty(regexp(text, '(^|[^a-z0-9])y[\s_\-/]*z([^a-z0-9]|$)', 'once'))
        planes = [planes, {'YZ'}];
    end
    if isempty(planes) && ...
            ~isempty(regexpi(text, '(cross[ -]?section|orthogonal scan)', 'once'))
        planes = {'XZ/YZ'};
    end
end

function [plane, reason, rawDirection] = planeFromLuminosaDirection(head)
    plane = 'unknown';
    rawDirection = [];

    if ~isfield(head, 'ImgHdr_ScanDirection')
        reason = 'ImgHdr_ScanDirection is absent.';
        return;
    end

    rawDirection = scalarFiniteDouble(head.ImgHdr_ScanDirection);
    if isempty(rawDirection) || rawDirection ~= round(rawDirection)
        reason = 'ImgHdr_ScanDirection is not a finite integer scalar.';
        return;
    end

    creator = '';
    if isfield(head, 'CreatorSW_Name')
        creator = textValue(head.CreatorSW_Name);
    end
    if isempty(regexpi(creator, 'luminosa', 'once'))
        reason = sprintf( ...
            ['ImgHdr_ScanDirection=%g is retained as raw metadata, but its ' ...
             'enumeration is not assumed for creator "%s".'], ...
            rawDirection, creator);
        return;
    end

    switch rawDirection
        case 0
            plane = 'XY';
            reason = [ ...
                'CreatorSW_Name is Luminosa and ImgHdr_ScanDirection=0; ' ...
                'the dataset-derived project convention classifies this as ' ...
                'XY (medium confidence; this is not an official PTU enum).'];
        case 1
            plane = 'YZ';
            reason = [ ...
                'CreatorSW_Name is Luminosa and ImgHdr_ScanDirection=1; ' ...
                'the dataset-derived project convention classifies this as ' ...
                'YZ (medium confidence; this is not an official PTU enum).'];
        case 2
            plane = 'XZ';
            reason = [ ...
                'CreatorSW_Name is Luminosa and ImgHdr_ScanDirection=2; ' ...
                'the dataset-derived project convention classifies this as ' ...
                'XZ (medium confidence; this is not an official PTU enum).'];
        otherwise
            reason = sprintf([ ...
                'Luminosa ImgHdr_ScanDirection=%g is outside the supported ' ...
                'project convention (0, 1, or 2).'], rawDirection);
    end
end

function value = scalarFiniteDouble(value)
    if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
        value = [];
        return;
    end
    value = double(value);
    if ~isfinite(value)
        value = [];
    end
end

function tf = planesConflict(a, b)
    tf = false;
    if strcmp(a, 'unknown') || strcmp(b, 'unknown')
        return;
    end
    if strcmp(a, b)
        return;
    end
    aCross = any(strcmp(a, {'XZ', 'YZ', 'XZ/YZ'}));
    bCross = any(strcmp(b, {'XZ', 'YZ', 'XZ/YZ'}));
    if aCross && bCross && (strcmp(a, 'XZ/YZ') || strcmp(b, 'XZ/YZ'))
        return;
    end
    tf = true;
end

function values = stableUnique(values)
    keep = true(size(values));
    for k = 2:numel(values)
        keep(k) = ~any(strcmp(values{k}, values(1:k - 1)));
    end
    values = values(keep);
end
