function export_cine_dicoms_to_gifs(rootFolder, outFolder, recursive, seriesNameFilter, num_f)
% EXPORT_CINE_DICOMS_TO_GIFS
% Find CINE DICOMs in a folder, group by series, and export GIFs.
%
% Usage:
%   export_cine_dicoms_to_gifs('D:\data\cine', 'D:\data\gifs', true);
%
% Inputs:
%   rootFolder : folder containing DICOMs
%   outFolder  : where GIFs will be written
%   recursive  : true/false (default true)
%
% Notes:
% - Groups by SeriesInstanceUID.
% - Uses TriggerTime to order cardiac phases when available.
% - Uses SliceLocation/ImagePositionPatient to order slices when available.
% - If series is not 2D+time, it will still try but may warn/skip.
if nargin < 5         
    num_f = 1; 
    specialNotes = '';
else
    specialNotes = cat(2, '_Echoes', num2str(num_f));
end
if nargin < 4                      , seriesNameFilter = ''; end
if nargin < 3 || isempty(recursive), recursive = true; end
if nargin < 2 || isempty(outFolder), outFolder = fullfile(rootFolder, 'cine_gifs'); end
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

% ---- Params ----
delay_sec  = 0.07;
delay_sec  = delay_sec/2;
loopcount  = inf;
do_montage = true;
do_perSlice = false; % set true if you want one GIF per slice as well
maxFilesScan = inf;  % optionally cap for testing

% ---- Collect candidate DICOM files ----
dicomFiles = list_files(rootFolder, recursive);
if isempty(dicomFiles)
    error('No files found under: %s', rootFolder);
end
if isfinite(maxFilesScan)
    dicomFiles = dicomFiles(1:min(numel(dicomFiles), maxFilesScan));
end

fprintf('Found %d files. Reading DICOM headers...\n', numel(dicomFiles));

% ---- Read headers; keep only likely cine images ----
entries = struct( ...
    'file', {}, 'seriesUID', {}, 'seriesNumber', {}, 'studyUID', {}, ...
    'rows', {}, 'cols', {}, 'instance', {}, 'trigger', {}, 'acqTime', {}, ...
    'slicePos', {}, 'ipp', {}, 'modality', {}, 'imageType', {}, 'cineLikely', {} );

for i = 1:numel(dicomFiles)
    f = dicomFiles{i};

    try
        info = dicominfo(f);
    catch
        continue; % not a dicom
    end

    % Must have pixel data-ish
    if ~isfield(info,'Rows') || ~isfield(info,'Columns')
        continue;
    end

    % Heuristic: cine-likely if
    % - ImageType contains 'CINE' or 'CARDIAC' or 'M'
    % - OR NumberOfFrames exists and > 1
    % - OR TriggerTime exists (common in cardiac cine)
    imageTypeStr = '';
    if isfield(info,'ImageType')
        if iscell(info.ImageType), imageTypeStr = strjoin(info.ImageType, '\'); 
        else, imageTypeStr = char(info.ImageType);
        end
    end
    nFrames = 1;
    if isfield(info,'NumberOfFrames')
        nFrames = double(info.NumberOfFrames);
    end

    cineLikely = false;
    if contains(upper(imageTypeStr), 'CINE') || contains(upper(imageTypeStr), 'CARD')
        cineLikely = true;
    end
    if nFrames > 1
        cineLikely = true;
    end
    if isfield(info,'TriggerTime')
        cineLikely = true;
    end

    % You can tighten this if your folder contains lots of non-cine:
    % if ~cineLikely, continue; end

    e.file = f;

    e.seriesUID = get_field(info,'SeriesInstanceUID','');
    e.seriesNumber = double(get_field(info,'SeriesNumber', NaN));
    e.studyUID  = get_field(info,'StudyInstanceUID','');

    e.rows = double(info.Rows);
    e.cols = double(info.Columns);

    e.instance = double(get_field(info,'InstanceNumber', NaN));
    e.trigger  = double(get_field(info,'TriggerTime', NaN)); % ms usually
    e.acqTime  = get_field(info,'AcquisitionTime','');

    % slice ordering helpers
    e.slicePos = double(get_field(info,'SliceLocation', NaN));
    if isfield(info,'ImagePositionPatient')
        e.ipp = double(info.ImagePositionPatient(:))'; % [x y z]
    else
        e.ipp = [NaN NaN NaN];
    end

    e.modality = get_field(info,'Modality','');
    e.imageType = imageTypeStr;
    e.cineLikely = cineLikely;

    entries(end+1) = e; %#ok<AGROW>
end

if isempty(entries)
    error('No readable DICOM image files found in: %s', rootFolder);
end

fprintf('Kept %d DICOM image files.\n', numel(entries));

% ---- Group by SeriesInstanceUID ----
uids = {entries.seriesUID};
[uniqUIDs, ~, idxUID] = unique(uids);

fprintf('Found %d series.\n', numel(uniqUIDs));

for s = 1:numel(uniqUIDs)
    seriesEntries = entries(idxUID == s);

    % Optional: skip obviously non-cine series
    % Require at least some cine-likelihood or repeated trigger/instance patterns

    cineScore = mean([seriesEntries.cineLikely]);
    if cineScore < 0.2
        % still might be cine, but likely not
        % continue; % uncomment to skip
    end

    % Read first header for naming
    try
        info0 = dicominfo(seriesEntries(num_f).file);
    catch
        continue;
    end

    seriesDesc = safe_str(get_field(info0,'SeriesDescription','Series'));

    % Skip if series name doesn't contain the filter string
    if ~isempty(seriesNameFilter) && ~contains(lower(seriesDesc), lower(seriesNameFilter))
        fprintf('  Skipping series: %s (does not contain "%s")\n', seriesDesc, seriesNameFilter);
        continue;
    end


    studyDesc  = safe_str(get_field(info0,'StudyDescription','Study'));
    patName    = '';
    if isfield(info0,'PatientName')
        try
            patName = safe_str(info0.PatientName.FamilyName);
        catch
            patName = '';
        end
    end

    uidShort = uniqUIDs{s};
    if numel(uidShort) > 8, uidShort = uidShort(end-7:end); end

    baseName = sprintf('%s_Ser%g_%s_UID%s%s', ...
        patName, get_field(info0,'SeriesNumber',NaN), seriesDesc, uidShort, specialNotes);
    baseName = regexprep(baseName, '[^\w\-]+', '_');
    if isempty(baseName), baseName = sprintf('Series_%d_UID%s%s', s, uidShort, specialNotes); end

    fprintf('\n[%d/%d] Processing: %s  (n=%d)\n', s, numel(uniqUIDs), baseName, numel(seriesEntries));

    % ---- Build 4D: [Ny Nx Nz Nt] from single-frame files (common cine) ----
    % Determine slice key
    sliceKey = nan(numel(seriesEntries),1);
    for i = 1:numel(seriesEntries)
        if ~isnan(seriesEntries(i).ipp(3))
            sliceKey(i) = seriesEntries(i).ipp(3);
        elseif ~isnan(seriesEntries(i).slicePos)
            sliceKey(i) = seriesEntries(i).slicePos;
        else
            sliceKey(i) = NaN;
        end
    end

    % Determine time key (cardiac phase)
    timeKey = nan(numel(seriesEntries),1);
    for i = 1:numel(seriesEntries)
        if ~isnan(seriesEntries(i).trigger)
            timeKey(i) = seriesEntries(i).trigger;
        else
            timeKey(i) = seriesEntries(i).instance; % fallback
        end
    end

    % If this is a multi-frame DICOM (NumberOfFrames > 1), handle separately
    nFrames0 = double(get_field(info0,'NumberOfFrames',1));
    if nFrames0 > 1
        fprintf('  Detected multi-frame DICOM (NumberOfFrames=%d). Reading frames...\n', nFrames0);
        try
            V = dicomread(info0); % could be [Rows Cols Frames] or [Rows Cols 1 Frames]
        catch ME
            warning('  Failed reading multi-frame: %s', ME.message);
            continue;
        end

        V = squeeze(V);
        if ndims(V) == 3
            % [Ny Nx Nt]
            recon = normalize01(V);
            % make it [Ny Nx 1 Nt] so montage code can work
            recon = reshape(recon, size(recon,1), size(recon,2), 1, size(recon,3));
            export_gifs_from_recon(recon, outFolder, baseName, do_montage, do_perSlice, delay_sec, loopcount);
        else
            warning('  Unexpected multi-frame dimensions: %s. Skipping.', mat2str(size(V)));
        end
        continue;
    end

    % ---- Single-frame files: infer unique slices and phases ----
    % Use tolerance because IPP z can be float with tiny jitter
    tolZ = 1e-3;
    [zVals, zIdx] = unique_tol(sliceKey, tolZ);
    [tVals, ~, tIdx] = unique(timeKey);

    Nz = numel(zVals);
    Nt = numel(tVals);

    if Nz < 2 || Nt < 2
        warning('  Not enough slices/phases to look like cine (Nz=%d, Nt=%d). Trying anyway...', Nz, Nt);
    end

    Ny = seriesEntries(1).rows;
    Nx = seriesEntries(1).cols;

    % Preallocate volume
    vol = zeros(Ny, Nx, Nz, Nt, 'single');
    hit = false(Nz, Nt);

    % Fill
    for i = 1:numel(seriesEntries)
        z = zIdx(i);
        t = tIdx(i);

        % If duplicates happen, keep the first
        if hit(z,t), continue; end

        try
            img = dicomread(seriesEntries(i).file);
        catch
            continue;
        end

        img = single(squeeze(img));
        if ~isequal(size(img), [Ny Nx])
            % Sometimes images stored transposed / different size; skip
            continue;
        end

        vol(:,:,z,t) = img;
        hit(z,t) = true;
    end

    % Check completeness
    completeness = nnz(hit) / numel(hit);
    fprintf('  Filled %d/%d frames (%.1f%%).\n', nnz(hit), numel(hit), 100*completeness);

    % Sort slices (z) ascending, phases (t) ascending
    [~, zOrder] = sort(zVals, 'ascend');
    [~, tOrder] = sort(tVals, 'ascend');
    vol = vol(:,:,zOrder,tOrder);

    % Normalize like your recon (abs then global max), but here magnitude already
    recon = normalize01(vol);

    % Optional shift like your code (fftshift in y)
    % recon = fftshift(recon, 1);

    % Export gifs
    export_gifs_from_recon(recon, outFolder, baseName, do_montage, do_perSlice, delay_sec, loopcount);
end

fprintf('\nDone. GIFs saved to: %s\n', outFolder);

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
% sanitize to short string
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

function export_gifs_from_recon(recon, outFolder, baseName, do_montage, do_perSlice, delay_sec, loopcount)
% recon: [Ny Nx Nz Nt] in [0..1]
[Ny, Nx, Nz, Nt] = size(recon);

% Grid for montage tiling
ncol = ceil(sqrt(Nz));
nrow = ceil(Nz / ncol);

if do_montage
    gif_name = fullfile(outFolder, sprintf('%s_montage.gif', baseName));
    for t = 1:Nt
        vol = recon(:,:,:,t); % [Ny Nx Nz]
        frame2d = imtile(vol, 'GridSize', [nrow ncol], 'BorderSize', [2 2]);
        frame2d = im2uint8(mat2gray(frame2d));
        [A,map] = gray2ind(frame2d, 256);

        if t == 1
            imwrite(A, map, gif_name, 'gif', 'LoopCount', loopcount, 'DelayTime', delay_sec);
        else
            imwrite(A, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
        end
    end
    fprintf('  Saved montage GIF: %s (Nt=%d, Nz=%d tiled %dx%d)\n', gif_name, Nt, Nz, nrow, ncol);
end

if do_perSlice
    outdir = fullfile(outFolder, sprintf('%s_per_slice', baseName));
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    for z = 1:Nz
        gifz = fullfile(outdir, sprintf('%s_slice_%02d.gif', baseName, z));
        for t = 1:Nt
            img = recon(:,:,z,t);
            img = im2uint8(mat2gray(img));
            [A,map] = gray2ind(img, 256);

            if t == 1
                imwrite(A, map, gifz, 'gif', 'LoopCount', inf, 'DelayTime', delay_sec);
            else
                imwrite(A, map, gifz, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
            end
        end
    end
    fprintf('  Saved per-slice GIFs in: %s\n', outdir);
end
end

function [valsU, idx] = unique_tol(vals, tol)
% unique with tolerance for floating values; NaNs treated as one bin
vals = vals(:);
idx = zeros(size(vals));
valsU = [];

nanMask = isnan(vals);
if any(nanMask)
    % put all NaNs into one "slice"
    valsU = [valsU; NaN];
    idx(nanMask) = numel(valsU);
end

finiteMask = ~nanMask;
vf = vals(finiteMask);
if isempty(vf)
    return;
end

[vs, ord] = sort(vf);
group = 1;
rep = vs(1);
valsU = [valsU; rep]; %#ok<AGROW>
idx_f = zeros(size(vs));
idx_f(1) = numel(valsU);

for k = 2:numel(vs)
    if abs(vs(k) - rep) <= tol
        idx_f(k) = numel(valsU);
    else
        rep = vs(k);
        valsU = [valsU; rep]; %#ok<AGROW>
        idx_f(k) = numel(valsU);
        group = group + 1; %#ok<NASGU>
    end
end

% map back to original order
idx_tmp = zeros(size(vf));
idx_tmp(ord) = idx_f;

idx(finiteMask) = idx_tmp;
end