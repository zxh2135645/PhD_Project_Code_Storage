function cineData = read_all_cine_dicoms(rootFolder, recursive, seriesNameFilter)
% READ_ALL_CINE_DICOMS
% Read all cine DICOM series from a folder into MATLAB memory.
%
% Usage:
%   cineData = read_all_cine_dicoms(rootFolder);
%   cineData = read_all_cine_dicoms(rootFolder, true, 'cine_sax');
%
% Inputs:
%   rootFolder       : folder containing DICOMs
%   recursive        : true/false (default true)
%   seriesNameFilter : optional substring filter for SeriesDescription
%
% Output:
%   cineData : struct array with fields
%       .vol               [Ny Nx Nz Nt] single
%       .baseName          generated readable name
%       .seriesUID
%       .seriesNumber
%       .seriesDescription
%       .studyDescription
%       .patientName
%       .Nz
%       .Nt
%       .zVals
%       .tVals
%       .files
%       .hitMap
%
% Notes:
% - Multi-frame DICOMs are returned as [Ny Nx 1 Nt]
% - Single-frame cine series are assembled into [Ny Nx Nz Nt]
% - Intensity is normalized to [0,1]

if nargin < 3
    seriesNameFilter = '';
end
if nargin < 2 || isempty(recursive)
    recursive = true;
end

% ---- Collect files ----
dicomFiles = list_files(rootFolder, recursive);
if isempty(dicomFiles)
    error('No files found under: %s', rootFolder);
end

fprintf('Found %d files. Reading DICOM headers...\n', numel(dicomFiles));

% ---- Read headers ----
entries = struct( ...
    'file', {}, 'seriesUID', {}, 'seriesNumber', {}, 'studyUID', {}, ...
    'rows', {}, 'cols', {}, 'instance', {}, 'trigger', {}, 'acqTime', {}, ...
    'slicePos', {}, 'ipp', {}, 'modality', {}, 'imageType', {}, 'cineLikely', {} );

for i = 1:numel(dicomFiles)
    f = dicomFiles{i};

    try
        info = dicominfo(f);
    catch
        continue;
    end

    if ~isfield(info, 'Rows') || ~isfield(info, 'Columns')
        continue;
    end

    imageTypeStr = '';
    if isfield(info, 'ImageType')
        if iscell(info.ImageType)
            imageTypeStr = strjoin(info.ImageType, '\');
        else
            imageTypeStr = char(info.ImageType);
        end
    end

    nFrames = 1;
    if isfield(info, 'NumberOfFrames')
        nFrames = double(info.NumberOfFrames);
    end

    cineLikely = false;
    if contains(upper(imageTypeStr), 'CINE') || contains(upper(imageTypeStr), 'CARD')
        cineLikely = true;
    end
    if nFrames > 1
        cineLikely = true;
    end
    if isfield(info, 'TriggerTime')
        cineLikely = true;
    end

    e.file = f;
    e.seriesUID    = get_field(info, 'SeriesInstanceUID', '');
    e.seriesNumber = double(get_field(info, 'SeriesNumber', NaN));
    e.studyUID     = get_field(info, 'StudyInstanceUID', '');

    e.rows = double(info.Rows);
    e.cols = double(info.Columns);

    e.instance = double(get_field(info, 'InstanceNumber', NaN));
    e.trigger  = double(get_field(info, 'TriggerTime', NaN));
    e.acqTime  = get_field(info, 'AcquisitionTime', '');

    e.slicePos = double(get_field(info, 'SliceLocation', NaN));
    if isfield(info, 'ImagePositionPatient')
        e.ipp = double(info.ImagePositionPatient(:))';
    else
        e.ipp = [NaN NaN NaN];
    end

    e.modality   = get_field(info, 'Modality', '');
    e.imageType  = imageTypeStr;
    e.cineLikely = cineLikely;

    entries(end+1) = e; %#ok<AGROW>
end

if isempty(entries)
    error('No readable DICOM image files found in: %s', rootFolder);
end

fprintf('Kept %d DICOM image files.\n', numel(entries));

% ---- Group by series ----
uids = {entries.seriesUID};
[uniqUIDs, ~, idxUID] = unique(uids);

fprintf('Found %d series.\n', numel(uniqUIDs));

cineData = struct( ...
    'vol', {}, ...
    'baseName', {}, ...
    'seriesUID', {}, ...
    'seriesNumber', {}, ...
    'seriesDescription', {}, ...
    'studyDescription', {}, ...
    'patientName', {}, ...
    'Nz', {}, ...
    'Nt', {}, ...
    'zVals', {}, ...
    'tVals', {}, ...
    'files', {}, ...
    'hitMap', {} );

for s = 1:numel(uniqUIDs)
    seriesEntries = entries(idxUID == s);

    try
        info0 = dicominfo(seriesEntries(1).file);
    catch
        continue;
    end

    seriesDesc = safe_str(get_field(info0, 'SeriesDescription', 'Series'));

    if ~isempty(seriesNameFilter) && ~contains(lower(seriesDesc), lower(seriesNameFilter))
        fprintf('Skipping series: %s (does not contain "%s")\n', seriesDesc, seriesNameFilter);
        continue;
    end

    studyDesc = safe_str(get_field(info0, 'StudyDescription', 'Study'));

    patName = '';
    if isfield(info0, 'PatientName')
        try
            patName = safe_str(info0.PatientName.FamilyName);
        catch
            patName = '';
        end
    end

    uidShort = uniqUIDs{s};
    if numel(uidShort) > 8
        uidShort = uidShort(end-7:end);
    end

    baseName = sprintf('%s_Ser%g_%s_UID%s', ...
        patName, get_field(info0, 'SeriesNumber', NaN), seriesDesc, uidShort);
    baseName = regexprep(baseName, '[^\w\-]+', '_');
    if isempty(baseName)
        baseName = sprintf('Series_%d_UID%s', s, uidShort);
    end

    fprintf('\n[%d/%d] Reading: %s (n=%d)\n', ...
        s, numel(uniqUIDs), baseName, numel(seriesEntries));

    % ---- Multi-frame handling ----
    nFrames0 = double(get_field(info0, 'NumberOfFrames', 1));
    if nFrames0 > 1
        fprintf('  Detected multi-frame DICOM (NumberOfFrames=%d)\n', nFrames0);

        try
            V = dicomread(info0);
        catch ME
            warning('  Failed reading multi-frame: %s', ME.message);
            continue;
        end

        V = squeeze(single(V));
        if ndims(V) == 3
            V = normalize01(V);
            V = reshape(V, size(V,1), size(V,2), 1, size(V,3));

            cineData(end+1).vol               = V; %#ok<AGROW>
            cineData(end).baseName            = baseName;
            cineData(end).seriesUID           = uniqUIDs{s};
            cineData(end).seriesNumber        = get_field(info0, 'SeriesNumber', NaN);
            cineData(end).seriesDescription   = seriesDesc;
            cineData(end).studyDescription    = studyDesc;
            cineData(end).patientName         = patName;
            cineData(end).Nz                  = 1;
            cineData(end).Nt                  = size(V,4);
            cineData(end).zVals               = 1;
            cineData(end).tVals               = 1:size(V,4);
            cineData(end).files               = {seriesEntries.file};
            cineData(end).hitMap              = true(1, size(V,4));
        else
            warning('  Unexpected multi-frame dimensions: %s. Skipping.', mat2str(size(V)));
        end

        continue;
    end

    % ---- Single-frame handling ----
    sliceKey = nan(numel(seriesEntries), 1);
    for i = 1:numel(seriesEntries)
        if ~isnan(seriesEntries(i).ipp(3))
            sliceKey(i) = seriesEntries(i).ipp(3);
        elseif ~isnan(seriesEntries(i).slicePos)
            sliceKey(i) = seriesEntries(i).slicePos;
        else
            sliceKey(i) = NaN;
        end
    end

    timeKey = nan(numel(seriesEntries), 1);
    for i = 1:numel(seriesEntries)
        if ~isnan(seriesEntries(i).trigger)
            timeKey(i) = seriesEntries(i).trigger;
        else
            timeKey(i) = seriesEntries(i).instance;
        end
    end

    tolZ = 1e-3;
    [zVals, zIdx] = unique_tol(sliceKey, tolZ);
    [tVals, ~, tIdx] = unique(timeKey);

    Nz = numel(zVals);
    Nt = numel(tVals);

    if Nz < 2 || Nt < 2
        warning('  Not enough slices/phases to look like cine (Nz=%d, Nt=%d). Reading anyway...', Nz, Nt);
    end

    Ny = seriesEntries(1).rows;
    Nx = seriesEntries(1).cols;

    vol = zeros(Ny, Nx, Nz, Nt, 'single');
    hit = false(Nz, Nt);

    for i = 1:numel(seriesEntries)
        z = zIdx(i);
        t = tIdx(i);

        if hit(z, t)
            continue;
        end

        try
            img = dicomread(seriesEntries(i).file);
        catch
            continue;
        end

        img = single(squeeze(img));
        if ~isequal(size(img), [Ny Nx])
            continue;
        end

        vol(:, :, z, t) = img;
        hit(z, t) = true;
    end

    completeness = nnz(hit) / numel(hit);
    fprintf('  Filled %d/%d frames (%.1f%%)\n', nnz(hit), numel(hit), 100 * completeness);

    [~, zOrder] = sort(zVals, 'ascend');
    [~, tOrder] = sort(tVals, 'ascend');

    vol = vol(:, :, zOrder, tOrder);
    zVals = zVals(zOrder);
    tVals = tVals(tOrder);
    hit = hit(zOrder, tOrder);

    vol = normalize01(vol);

    cineData(end+1).vol             = vol; %#ok<AGROW>
    cineData(end).baseName          = baseName;
    cineData(end).seriesUID         = uniqUIDs{s};
    cineData(end).seriesNumber      = get_field(info0, 'SeriesNumber', NaN);
    cineData(end).seriesDescription = seriesDesc;
    cineData(end).studyDescription  = studyDesc;
    cineData(end).patientName       = patName;
    cineData(end).Nz                = Nz;
    cineData(end).Nt                = Nt;
    cineData(end).zVals             = zVals;
    cineData(end).tVals             = tVals;
    cineData(end).files             = {seriesEntries.file};
    cineData(end).hitMap            = hit;
end

fprintf('\nDone. Read %d cine series.\n', numel(cineData));

end

% ===================== Helpers =====================

function files = list_files(rootFolder, recursive)
if recursive
    d = dir(fullfile(rootFolder, '**', '*'));
else
    d = dir(fullfile(rootFolder, '*'));
end
d = d(~[d.isdir]);
files = fullfile({d.folder}, {d.name});
end

function v = get_field(s, field, defaultVal)
if isfield(s, field)
    v = s.(field);
else
    v = defaultVal;
end
end

function s = safe_str(x)
if isempty(x)
    s = '';
    return;
end
if isstring(x) || ischar(x)
    s = char(x);
else
    try
        s = char(string(x));
    catch
        s = '';
    end
end
s = strtrim(s);
end

function X = normalize01(X)
X = single(X);
mx = max(X(:));
if mx <= 0 || ~isfinite(mx)
    X = zeros(size(X), 'single');
else
    X = X ./ mx;
end
end

function [valsU, idx] = unique_tol(vals, tol)
vals = vals(:);
idx = zeros(size(vals));
valsU = [];

nanMask = isnan(vals);
if any(nanMask)
    valsU = [valsU; NaN];
    idx(nanMask) = numel(valsU);
end

finiteMask = ~nanMask;
vf = vals(finiteMask);
if isempty(vf)
    return;
end

[vs, ord] = sort(vf);
rep = vs(1);
valsU = [valsU; rep];
idx_f = zeros(size(vs));
idx_f(1) = numel(valsU);

for k = 2:numel(vs)
    if abs(vs(k) - rep) <= tol
        idx_f(k) = numel(valsU);
    else
        rep = vs(k);
        valsU = [valsU; rep];
        idx_f(k) = numel(valsU);
    end
end

idx_tmp = zeros(size(vf));
idx_tmp(ord) = idx_f;
idx(finiteMask) = idx_tmp;
end