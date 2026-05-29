function outputs = run_lrt_ddce_pipeline(config)
% run_lrt_ddce_pipeline  Fit LRT T1 maps, dDCE maps, and synthetic late enhancement.
%
%   outputs = run_lrt_ddce_pipeline(config)
%
% See default_lrt_ddce_pipeline_config for the expected config fields.

if nargin < 1 || isempty(config)
    config = default_lrt_ddce_pipeline_config();
else
    config = mergeStructs(default_lrt_ddce_pipeline_config(), config);
end

config = validatePipelineConfig(config);
ensureDir(config.outputDir);
ensureDir(fullfile(config.outputDir, 'masks'));
ensureDir(fullfile(config.outputDir, 'qc'));

if exist(config.dceCodeDir, 'dir')
    addpath(config.dceCodeDir);
else
    error('run_lrt_ddce_pipeline:MissingDceCodeDir', ...
        'DCE code directory not found: %s', config.dceCodeDir);
end

save(fullfile(config.outputDir, 'pipeline_config_used.mat'), 'config');

[t1Stage, dceInput, dceStage] = getPipelineStages(config);
synthStage = runAndSaveSynthesis(dceStage, config);
synthDicomReport = exportSynthDicomIfRequested(synthStage, config);

if config.saveQc
    fprintf('LRT dDCE pipeline: writing QC images...\n');
    writePipelineQc(t1Stage, dceStage, synthStage, config);
end

outputs = struct();
outputs.config = config;
outputs.t1Stage = t1Stage;
outputs.dceInput = dceInput;
outputs.dceStage = dceStage;
outputs.synthStage = synthStage;
outputs.synthDicomReport = synthDicomReport;

fprintf('LRT dDCE pipeline complete: %s\n', config.outputDir);

end

function [t1Stage, dceInput, dceStage] = getPipelineStages(config)

switch lower(config.startStage)
    case 't1'
        t1Stage = runAndSaveT1Stage(config);
        [dceInput, dceStage] = runAndSaveDceStage(t1Stage, config);
    case 'dce'
        fprintf('LRT dDCE pipeline: loading T1 maps from previous run...\n');
        t1Stage = loadT1Stage(config);
        [dceInput, dceStage] = runAndSaveDceStage(t1Stage, config);
    case 'synthesis'
        t1Stage = struct();
        fprintf('LRT dDCE pipeline: loading dDCE outputs from previous run...\n');
        [dceInput, dceStage] = loadDceStage(config);
    otherwise
        error('run_lrt_ddce_pipeline:InvalidStartStage', ...
            'config.startStage must be ''t1'', ''dce'', or ''synthesis''.');
end

end

function t1Stage = runAndSaveT1Stage(config)

fprintf('LRT dDCE pipeline: fitting T1 maps...\n');
t1Stage = fitLrtT1Maps(config);
t1MapMs = t1Stage.t1MapMs;
t1Mask = t1Stage.t1Mask;
dictionary = t1Stage.dictionary;
sequence = t1Stage.sequence;
save(fullfile(config.outputDir, 't1_maps.mat'), ...
    't1MapMs', 't1Mask', 'dictionary', 'sequence', 't1Stage', '-v7.3');

end

function [dceInput, dceStage] = runAndSaveDceStage(t1Stage, config)

fprintf('LRT dDCE pipeline: converting T1 maps to Gd concentration...\n');
dceInput = prepareDceInput(t1Stage, config);
dceStage = runDceStage(dceInput, config);
Gdcon = dceStage.Gdcon;
dR1 = dceStage.dR1;
precontrastT1Ms = dceStage.precontrastT1Ms;
dceSliceIndices = dceStage.dceSliceIndices;
roiMasksBySlice = dceStage.roiMasksBySlice;
sliceResults = dceStage.sliceResults;
save(fullfile(config.outputDir, 'gd_concentration.mat'), ...
    'dceInput', 'Gdcon', 'dR1', 'precontrastT1Ms', ...
    'dceSliceIndices', 'roiMasksBySlice', 'sliceResults', '-v7.3');
aifBySlice = dceStage.aifBySlice;
curvesBySlice = dceStage.curvesBySlice;
save(fullfile(config.outputDir, 'aif_fit.mat'), 'aifBySlice', 'curvesBySlice', '-v7.3');

if isfield(dceStage, 'maps2CXM')
    maps2CXM = dceStage.maps2CXM;
    save(fullfile(config.outputDir, 'ddce_2cxm_maps.mat'), 'maps2CXM', '-v7.3');
end

if isfield(dceStage, 'mapsETK')
    mapsETK = dceStage.mapsETK;
    save(fullfile(config.outputDir, 'ddce_etk_maps.mat'), 'mapsETK', '-v7.3');
end

end

function synthStage = runAndSaveSynthesis(dceStage, config)

if isfield(dceStage, 'maps2CXM')
    fprintf('LRT dDCE pipeline: synthesizing late enhancement...\n');
    synthStage = synthesizeLateEnhancement(dceStage, config);
    synthTimeMinutes = synthStage.timeMinutes;
    synthBySlice = synthStage.synthBySlice;
    save(fullfile(config.outputDir, 'synth_lge_2cxm.mat'), ...
        'synthTimeMinutes', 'synthBySlice', 'synthStage', '-v7.3');
else
    synthStage = struct();
end

end

function t1Stage = fitLrtT1Maps(config)

data = load(config.reconMatPath, 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'params');
requiredFields = {'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'params'};
assertFields(data, requiredFields, 'reconstruction MAT file');

sequence = getSequencePreset(config.sequencePreset, data.params, config.sequence);
if isempty(config.sliceIndices)
    sliceIndices = 1:data.Nz;
else
    sliceIndices = config.sliceIndices(:).';
end
sliceIndices = validateIndices(sliceIndices, data.Nz, 'sliceIndices');

if isempty(config.timeIndices)
    timeIndices = 1:size(data.Phi, 5);
else
    timeIndices = config.timeIndices(:).';
end
timeIndices = validateIndices(timeIndices, size(data.Phi, 5), 'timeIndices');
cardiacPhases = validateIndices(config.cardiacPhases(:).', size(data.Phi, 3), 'cardiacPhases');
respPhase = validateIndices(config.respPhase, size(data.Phi, 4), 'respPhase');

t1Mask = loadOrCreateT1Mask(data, sliceIndices, cardiacPhases(1), respPhase, sequence, config);
dict = buildT1Dictionary(sequence, data.params);

t1MapMs = zeros(data.Ny, data.Nx, data.Nz, numel(timeIndices), numel(cardiacPhases));
Uimg = reshape(data.U, data.Ny, data.Nx, data.Nz, []);

for ntIdx = 1:numel(timeIndices)
    nt = timeIndices(ntIdx);
    for sliceIdx = 1:numel(sliceIndices)
        slc = sliceIndices(sliceIdx);
        dispim = @(x) fftshift(x(:,:,slc,:), 1);
        sliceBasis = reshape(dispim(Uimg), [], data.L);
        for cardIdx = 1:numel(cardiacPhases)
            card = cardiacPhases(cardIdx);
            coef = data.Gr \ reshape(data.Phi(:,:,card,respPhase,nt), data.L, []);
            recon = reshape(sliceBasis * coef, data.Ny, data.Nx, [], sequence.reconNEco);
            recon = squeeze(recon(:,:,:,1));
            fitTrain = cropRecoveryTrain(recon, sequence);
            mask1d = t1Mask(:,:,slc);
            t1MapMs(:,:,slc,ntIdx,cardIdx) = fitT1Dictionary(abs(fitTrain), mask1d(:), dict, data.Ny, data.Nx);
        end
    end
end

t1MapMs = fftshift(t1MapMs, 3);
maskShifted = fftshift(t1Mask, 3);

t1Stage = struct();
t1Stage.t1MapMs = t1MapMs;
t1Stage.t1Mask = maskShifted;
t1Stage.dictionary = dict;
t1Stage.sequence = sequence;
t1Stage.sliceIndices = sliceIndices;
t1Stage.timeIndices = timeIndices;
t1Stage.cardiacPhases = cardiacPhases;
t1Stage.respPhase = respPhase;

end

function dceInput = prepareDceInput(t1Stage, config)

if isstruct(t1Stage)
    t1MapMs = t1Stage.t1MapMs;
    if isfield(t1Stage, 't1Mask')
        t1Mask = t1Stage.t1Mask;
    else
        t1Mask = [];
    end
else
    t1MapMs = t1Stage;
    t1Mask = [];
end
t1ForDce = applyOrientationTransform(t1MapMs, config.orientationTransform);
if ndims(t1ForDce) == 5
    t1ForDce = mean(t1ForDce, 5);
end
if ndims(t1ForDce) == 3
    t1ForDce = reshape(t1ForDce, size(t1ForDce, 1), size(t1ForDce, 2), size(t1ForDce, 3), 1);
end
if ~isempty(t1Mask)
    t1MaskForDce = logical(applyOrientationTransform(t1Mask, config.orientationTransform));
else
    t1MaskForDce = [];
end

preT1 = loadPrecontrastT1(t1ForDce, config);
postT1 = t1ForDce;
if numel(config.timeMinutes) < size(postT1, 4)
    postT1 = postT1(:,:,:,1:numel(config.timeMinutes));
end

prePostT1 = cat(4, preT1, postT1);
dR1 = (1 ./ (prePostT1(:,:,:,2:end) + eps) - 1 ./ (prePostT1(:,:,:,1) + eps)) * 1e3;
Gdcon = dR1 / config.gdRelaxivity;
Gdcon(Gdcon < 0) = 0;
Gdcon(~isfinite(Gdcon)) = 0;

dceInput = struct();
dceInput.t1MapMs = t1ForDce;
dceInput.precontrastT1Ms = preT1;
dceInput.dR1 = dR1;
dceInput.Gdcon = Gdcon;
dceInput.t1Mask = t1MaskForDce;

end

function dceStage = runDceStage(dceInput, config)

dceTimer = tic;
Gdcon = dceInput.Gdcon;
if isempty(config.dceSliceIndices)
    dceSliceIndices = chooseDceSlicesFromMask(dceInput, size(Gdcon, 3));
else
    dceSliceIndices = validateIndices(config.dceSliceIndices(:).', size(Gdcon, 3), 'dceSliceIndices');
end

run2CXM = any(strcmpi(config.models, '2CXM'));
runETK = any(strcmpi(config.models, 'ETK'));

if run2CXM
    maps2CXM = initSliceMapStruct(size(Gdcon), {'F', 'PS', 'vp', 've', 'Rsq', 'RMSE'});
end
if runETK
    mapsETK = initSliceMapStruct(size(Gdcon), {'Ktrans', 'Kep', 'Ve', 'Vp', 'Rsq', 'RMSE'});
end

sliceResults = struct([]);
aifBySlice = struct([]);
curvesBySlice = struct([]);
roiMasksBySlice = struct([]);
roiTimingSec = zeros(1, numel(dceSliceIndices));

fprintf('LRT dDCE pipeline: collecting ROIs for %d DCE slice(s)...\n', numel(dceSliceIndices));
for idx = 1:numel(dceSliceIndices)
    slc = dceSliceIndices(idx);
    roiTimer = tic;
    GdSlice = squeeze(Gdcon(:,:,slc,:));
    fprintf('  ROI slice %d (%d/%d)...\n', slc, idx, numel(dceSliceIndices));
    roiMasksBySlice(idx).sliceIndex = slc;
    roiMasksBySlice(idx).roiMasks = loadOrCreateDceRois(GdSlice, slc, config);
    roiTimingSec(idx) = toc(roiTimer);
    fprintf('  ROI slice %d complete in %.1f sec.\n', slc, roiTimingSec(idx));
end

fprintf('LRT dDCE pipeline: fitting requested DCE model(s): %s\n', strjoin(config.models, ', '));
for idx = 1:numel(dceSliceIndices)
    slc = dceSliceIndices(idx);
    sliceTimer = tic;
    timingSec = struct('roi', roiTimingSec(idx), 'prep', 0, 'fit2CXM', 0, ...
        'fitETK', 0, 'total', 0);
    GdSlice = squeeze(Gdcon(:,:,slc,:));
    roiMasks = roiMasksBySlice(idx).roiMasks;
    qualityMask = dceFitQualityMask(roiMasks);

    fprintf('  DCE slice %d (%d/%d): preparing AIF/curves...\n', slc, idx, numel(dceSliceIndices));
    prepTimer = tic;
    [Ctoi, Cp, time, fitPointInd] = buildDceCurves(GdSlice, roiMasks, config);
    [aif, aifFit] = fitAif(Cp, time, config);
    curves = summarizeDceCurves(GdSlice, roiMasks, config);
    timingSec.prep = toc(prepTimer);

    sliceResults(idx).sliceIndex = slc;
    sliceResults(idx).fitPointInd = fitPointInd;
    sliceResults(idx).timeMinutes = time;
    sliceResults(idx).Cp = Cp;
    sliceResults(idx).aif = aif;

    if run2CXM
        fprintf('  DCE slice %d: fitting 2CXM...\n', slc);
        fitTimer = tic;
        fit2 = fit2CxmSlice(Ctoi, aif, time, roiMasks.HeartMask, qualityMask, config);
        timingSec.fit2CXM = toc(fitTimer);
        maps2CXM.F(:,:,slc) = fit2.F;
        maps2CXM.PS(:,:,slc) = fit2.PS;
        maps2CXM.vp(:,:,slc) = fit2.vp;
        maps2CXM.ve(:,:,slc) = fit2.ve;
        maps2CXM.Rsq(:,:,slc) = fit2.Rsq;
        maps2CXM.RMSE(:,:,slc) = fit2.RMSE;
        sliceResults(idx).fit2CXM = fit2;
        fprintf('  DCE slice %d: 2CXM fit complete in %.1f sec.\n', slc, timingSec.fit2CXM);
    end

    if runETK
        fprintf('  DCE slice %d: fitting ETK...\n', slc);
        fitTimer = tic;
        fitE = fitEtkSlice(Ctoi, aif, time, roiMasks.HeartMask, qualityMask, config);
        timingSec.fitETK = toc(fitTimer);
        mapsETK.Ktrans(:,:,slc) = fitE.Ktrans;
        mapsETK.Kep(:,:,slc) = fitE.Kep;
        mapsETK.Vp(:,:,slc) = fitE.Vp;
        mapsETK.Ve(:,:,slc) = fitE.Ve;
        mapsETK.Rsq(:,:,slc) = fitE.Rsq;
        mapsETK.RMSE(:,:,slc) = fitE.RMSE;
        sliceResults(idx).fitETK = fitE;
        fprintf('  DCE slice %d: ETK fit complete in %.1f sec.\n', slc, timingSec.fitETK);
    end

    aifBySlice(idx).sliceIndex = slc;
    aifBySlice(idx).fit = aifFit;
    curvesBySlice(idx).sliceIndex = slc;
    curvesBySlice(idx).curves = curves;
    timingSec.total = roiTimingSec(idx) + toc(sliceTimer);
    sliceResults(idx).timingSec = timingSec;
    fprintf('  DCE slice %d complete in %.1f sec (ROI %.1f, prep %.1f, 2CXM %.1f, ETK %.1f).\n', ...
        slc, timingSec.total, timingSec.roi, timingSec.prep, timingSec.fit2CXM, timingSec.fitETK);
end

dceStage = struct();
dceStage.Gdcon = Gdcon;
dceStage.dR1 = dceInput.dR1;
dceStage.precontrastT1Ms = dceInput.precontrastT1Ms;
dceStage.dceSliceIndices = dceSliceIndices;
dceStage.sliceResults = sliceResults;
dceStage.aifBySlice = aifBySlice;
dceStage.curvesBySlice = curvesBySlice;
dceStage.roiMasksBySlice = roiMasksBySlice;
if run2CXM
    dceStage.maps2CXM = maps2CXM;
end
if runETK
    dceStage.mapsETK = mapsETK;
end
dceStage.timingSec = struct('total', toc(dceTimer));
fprintf('LRT dDCE pipeline: DCE stage complete in %.1f sec.\n', dceStage.timingSec.total);

end

function synthStage = synthesizeLateEnhancement(dceStage, config)

timeMinutes = 0:config.synthesisStepMinutes:config.synthesisEndMinutes;
synthBySlice = struct([]);
scaleFactors = config.synthesisScaleFactors(:).';

for idx = 1:numel(dceStage.dceSliceIndices)
    slc = dceStage.dceSliceIndices(idx);
    roiMasks = dceStage.roiMasksBySlice(idx).roiMasks;
    aifFit = dceStage.aifBySlice(idx).fit;
    timeInterp = 0:(1/60):max(timeMinutes);
    synthCell = cell(numel(scaleFactors), 1);

    for sfIdx = 1:numel(scaleFactors)
        aifParams = aifFit.fittedParamsStruct_gauss;
        aifParams.Sf = scaleFactors(sfIdx);
        aif = gaussianExpP1Only(aifParams, {aifFit.fittedParamsStruct.a1, ...
            aifFit.fittedParamsStruct.sigma1, aifFit.fittedParamsStruct.mu1, ...
            aifFit.fittedParamsStruct.a2, aifFit.fittedParamsStruct.sigma2, ...
            aifFit.fittedParamsStruct.mu2, timeInterp});
        aif(1) = 0;

        synthImg = zeros(size(dceStage.Gdcon, 1), size(dceStage.Gdcon, 2), numel(timeMinutes));
        heartIdx = find(roiMasks.HeartMask);
        Fmap = dceStage.maps2CXM.F(:,:,slc);
        PSmap = dceStage.maps2CXM.PS(:,:,slc);
        vpMap = dceStage.maps2CXM.vp(:,:,slc);
        veMap = dceStage.maps2CXM.ve(:,:,slc);
        for c = 1:numel(heartIdx)
            linIdx = heartIdx(c);
            F = Fmap(linIdx);
            PS = PSmap(linIdx);
            vp = vpMap(linIdx);
            ve = veMap(linIdx);
            curve = cxmBoundCurve([F * scaleFactors(sfIdx), PS, vp, ve], {timeMinutes(:), aif(:), 1});
            [r, col] = ind2sub(size(roiMasks.HeartMask), linIdx);
            synthImg(r,col,:) = reshape(curve, 1, 1, []);
        end
        synthCell{sfIdx} = synthImg;
    end

    synthBySlice(idx).sliceIndex = slc;
    synthBySlice(idx).scaleFactors = scaleFactors;
    synthBySlice(idx).synthLGE = synthCell;
end

synthStage = struct();
synthStage.timeMinutes = timeMinutes;
synthStage.synthBySlice = synthBySlice;

end

function report = exportSynthDicomIfRequested(synthStage, config)

report = struct('enabled', config.exportSynthDicom, 'filesWritten', 0, ...
    'filesSkipped', 0, 'outputDir', config.synthDicomOutputDir);
if ~config.exportSynthDicom
    return;
end
if ~isfield(synthStage, 'synthBySlice') || isempty(synthStage.synthBySlice)
    fprintf('LRT dDCE pipeline: synthetic DICOM export requested, but no synthetic LGE images are available.\n');
    return;
end

fprintf('LRT dDCE pipeline: exporting synthetic LGE 2CXM DICOM...\n');
headers = synthDicomLgeHeaders(config);
ensureDir(config.synthDicomOutputDir);
[filesWritten, filesSkipped] = writeSynthLgeDicom(synthStage, headers, config);
report.filesWritten = filesWritten;
report.filesSkipped = filesSkipped;
fprintf('LRT dDCE pipeline: synthetic DICOM export complete: wrote %d, skipped %d existing. Output: %s\n', ...
    filesWritten, filesSkipped, config.synthDicomOutputDir);

end

function headers = synthDicomLgeHeaders(config)

fid = extractFidForSynthDicom(config);
subjectId = synthDicomSubjectId(fid, config);
subjectDir = fullfile(config.synthDicomRoot, subjectId);
if ~exist(subjectDir, 'dir')
    error('run_lrt_ddce_pipeline:MissingSynthDicomSubjectDir', ...
        'No DICOM subject folder for FID%s subject %s: %s', fid, subjectId, subjectDir);
end

dicomFiles = dir(fullfile(subjectDir, '**', '*'));
dicomFiles = dicomFiles(~[dicomFiles.isdir]);
dicomFiles = dicomFiles(~cellfun(@(name) strncmp(name, '._', 2), {dicomFiles.name}));
lgeHeaders = {};
readCount = 0;
for k = 1:numel(dicomFiles)
    filePath = fullfile(dicomFiles(k).folder, dicomFiles(k).name);
    try
        info = dicominfo(filePath);
    catch
        continue;
    end
    readCount = readCount + 1;
    if isfield(info, 'SeriesDescription') && ...
            contains(strtrim(info.SeriesDescription), config.synthDicomLgeSeriesDescriptionContains, 'IgnoreCase', true)
        lgeHeaders{end+1} = info; %#ok<AGROW>
    end
end

if isempty(lgeHeaders)
    error('run_lrt_ddce_pipeline:MissingSynthDicomLgeHeaders', ...
        'Read %d DICOM header(s) for subject %s, but found no LGE series containing "%s".', ...
        readCount, subjectId, config.synthDicomLgeSeriesDescriptionContains);
end

headers = dedupeHeadersBySlice(sortHeadersBySlice(lgeHeaders));
fprintf('LRT dDCE pipeline: selected %d LGE header slice(s) for synthetic DICOM export from subject %s.\n', ...
    numel(headers), subjectId);

end

function fid = extractFidForSynthDicom(config)

candidates = {config.reconMatPath, config.outputDir};
fid = '';
for k = 1:numel(candidates)
    tok = regexp(candidates{k}, 'FID(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        fid = tok{1};
        return;
    end
end
error('run_lrt_ddce_pipeline:MissingSynthDicomFid', ...
    'Could not extract FID from config.reconMatPath or config.outputDir.');

end

function subjectId = synthDicomSubjectId(fid, config)

lookupMap = buildSynthDicomLookupMap(config.synthDicomLookupFile, 3, 1);
lookupKey = canonicalDigits(fid);
if ~isKey(lookupMap, lookupKey)
    error('run_lrt_ddce_pipeline:MissingSynthDicomLookupRow', ...
        'No lookup row in %s for FID%s.', config.synthDicomLookupFile, fid);
end
subjectId = lookupMap(lookupKey);

end

function lookupMap = buildSynthDicomLookupMap(lookupFile, matFIDColumn, dicomSubjectColumn)

T = readtable(lookupFile, 'ReadVariableNames', true, 'PreserveVariableNames', true);
matVals = getTableColumn(T, matFIDColumn);
subjVals = getTableColumn(T, dicomSubjectColumn);

lookupMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:height(T)
    matKey = canonicalDigits(matVals(i));
    subj = canonicalSubjectId(subjVals(i));
    if ~isempty(matKey) && ~isempty(subj)
        lookupMap(matKey) = subj;
    end
end

end

function col = getTableColumn(T, columnId)

if isnumeric(columnId)
    col = T{:, columnId};
else
    col = T.(columnId);
end

end

function [filesWritten, filesSkipped] = writeSynthLgeDicom(synthStage, headers, config)

filesWritten = 0;
filesSkipped = 0;
timeMinutes = synthStage.timeMinutes(:).';
sliceIndices = [synthStage.synthBySlice.sliceIndex];
scaleFactors = synthStage.synthBySlice(1).scaleFactors(:).';
fid = extractFidForSynthDicom(config);

for sfIdx = 1:numel(scaleFactors)
    scaleLabel = dicomLabel(scaleFactors(sfIdx));
    scaleDir = fullfile(config.synthDicomOutputDir, sprintf('LGE_synth_Scale_%s', scaleLabel));
    ensureDir(scaleDir);
    seriesUID = dicomuid();

    for sliceOrderIdx = 1:numel(synthStage.synthBySlice)
        slc = synthStage.synthBySlice(sliceOrderIdx).sliceIndex;
        sliceDir = fullfile(scaleDir, sprintf('Slice_%02d', slc));
        ensureDir(sliceDir);
        synthCell = synthStage.synthBySlice(sliceOrderIdx).synthLGE;
        if sfIdx > numel(synthCell)
            continue;
        end
        synthImg = synthCell{sfIdx};
        sliceScale = synthDicomSliceScale(synthImg);
        headerIdx = synthDicomHeaderIndex(slc, sliceOrderIdx, sliceIndices, numel(headers));

        for timeIdx = 1:numel(timeMinutes)
            timeLabel = dicomLabel(timeMinutes(timeIdx));
            outFile = fullfile(sliceDir, sprintf('FID%s_SynthLGE2CXM_Scale%s_Slice%02d_Time%smin.dcm', ...
                fid, scaleLabel, slc, timeLabel));
            if exist(outFile, 'file') == 2
                filesSkipped = filesSkipped + 1;
                continue;
            end

            meta = headers{headerIdx};
            meta.SeriesDescription = sprintf('LRT_SynthLGE_2CXM_Scale%s', scaleLabel);
            if isfield(meta, 'SeriesNumber') && isnumeric(meta.SeriesNumber)
                meta.SeriesNumber = double(meta.SeriesNumber) + 5000 + sfIdx;
            end
            meta.SeriesInstanceUID = seriesUID;
            meta.SOPInstanceUID = dicomuid();
            meta.MediaStorageSOPInstanceUID = meta.SOPInstanceUID;
            meta.InversionTime = timeMinutes(timeIdx) * 60 * 1000;
            meta.InstanceNumber = (sliceOrderIdx - 1) * numel(timeMinutes) + timeIdx;
            meta.AcquisitionNumber = meta.InstanceNumber;

            frameImg = prepareSynthDicomImage(synthImg(:,:,timeIdx), meta, config, sliceScale);
            meta = updateSynthDicomPixelMetadata(meta, frameImg, config);
            dicomwrite(frameImg, outFile, meta, 'CreateMode', 'copy');
            filesWritten = filesWritten + 1;
        end
    end
end

end

function idx = synthDicomHeaderIndex(sliceIndex, sliceOrderIdx, sliceIndices, nHeaders)

if nHeaders >= max(sliceIndices) && sliceIndex >= 1
    idx = sliceIndex;
elseif nHeaders >= numel(sliceIndices)
    idx = sliceOrderIdx;
else
    idx = min(sliceOrderIdx, nHeaders);
end
idx = min(max(idx, 1), nHeaders);

end

function sliceScale = synthDicomSliceScale(synthImg)

sliceScale = max(abs(synthImg(:)));
if ~isfinite(sliceScale) || sliceScale <= 0
    sliceScale = 1;
end

end

function frameImg = prepareSynthDicomImage(frameImg, meta, config, sliceScale)

frameImg = abs(frameImg);
if isfield(meta, 'Rows') && size(frameImg, 1) > double(meta.Rows)
    targetRows = double(meta.Rows);
    rowStart = floor((size(frameImg, 1) - targetRows) / 2) + 1;
    rowEnd = rowStart + targetRows - 1;
    frameImg = frameImg(rowStart:rowEnd, :);
end

frameImg = imrotate(frameImg, 90);
frameImg = flip(frameImg, 2);
frameImg = frameImg ./ max(sliceScale, eps);
frameImg = min(max(frameImg, 0), 1);
frameImg = uint16(frameImg * config.synthDicomScaleMax);

end

function meta = updateSynthDicomPixelMetadata(meta, frameImg, config)

meta.SmallestImagePixelValue = min(frameImg(:));
meta.LargestImagePixelValue = max(frameImg(:));
meta.WindowCenter = config.synthDicomWindowCenter;
meta.WindowWidth = config.synthDicomWindowWidth;
meta.Rows = size(frameImg, 1);
meta.Columns = size(frameImg, 2);

end

function label = dicomLabel(value)

label = sprintf('%g', value);
label = strrep(label, '-', 'm');
label = strrep(label, '.', 'p');

end

function headers = sortHeadersBySlice(headers)

locs = nan(1, numel(headers));
for k = 1:numel(headers)
    locs(k) = dicomSliceLocation(headers{k});
end
[~, order] = sort(locs, 'ascend');
headers = headers(order);

end

function headersOut = dedupeHeadersBySlice(headers)

headersOut = {};
locKeys = {};
for k = 1:numel(headers)
    loc = dicomSliceLocation(headers{k});
    if isnan(loc)
        key = sprintf('no_loc_%04d', k);
    else
        key = sprintf('%.4f', loc);
    end
    if ~any(strcmp(locKeys, key))
        locKeys{end+1} = key; %#ok<AGROW>
        headersOut{end+1} = headers{k}; %#ok<AGROW>
    end
end

end

function loc = dicomSliceLocation(info)

loc = NaN;
if isfield(info, 'SliceLocation') && ~isempty(info.SliceLocation)
    loc = double(info.SliceLocation);
elseif isfield(info, 'ImagePositionPatient') && numel(info.ImagePositionPatient) >= 3
    loc = double(info.ImagePositionPatient(3));
end

end

function key = canonicalDigits(value)

s = scalarToText(value);
pieces = regexp(s, '\d+', 'match');
if isempty(pieces)
    key = '';
else
    key = pieces{1};
    key = regexprep(key, '^0+', '');
    if isempty(key)
        key = '0';
    end
end

end

function subjectId = canonicalSubjectId(value)

digitsOnly = regexp(scalarToText(value), '\d', 'match');
if isempty(digitsOnly)
    subjectId = '';
    return;
end

s = [digitsOnly{:}];
if numel(s) > 10
    s = s(end-9:end);
elseif numel(s) < 10
    s = [repmat('0', 1, 10 - numel(s)), s];
end
subjectId = s;

end

function s = scalarToText(value)

if iscell(value)
    if isempty(value) || isempty(value{1})
        s = '';
        return;
    end
    value = value{1};
end

if isnumeric(value)
    if isempty(value) || any(isnan(value(:)))
        s = '';
    else
        s = sprintf('%.0f', value(1));
    end
elseif ischar(value)
    s = value;
elseif supportsStringArrays() && isstring(value)
    if ismissing(value)
        s = '';
    else
        s = char(value);
    end
else
    try
        s = char(string(value));
    catch
        s = '';
    end
end

end

function t1Stage = loadT1Stage(config)

pathIn = requireRestartFile(config, 't1_maps.mat', 'dce');
loaded = load(pathIn);
if isfield(loaded, 't1Stage')
    t1Stage = loaded.t1Stage;
else
    assertFields(loaded, {'t1MapMs', 't1Mask', 'dictionary', 'sequence'}, 't1_maps.mat');
    t1Stage = struct();
    t1Stage.t1MapMs = loaded.t1MapMs;
    t1Stage.t1Mask = loaded.t1Mask;
    t1Stage.dictionary = loaded.dictionary;
    t1Stage.sequence = loaded.sequence;
end
assertFields(t1Stage, {'t1MapMs'}, 'loaded t1Stage');
t1Stage = refreshLoadedT1Mask(t1Stage, config);

end

function t1Stage = refreshLoadedT1Mask(t1Stage, config)

maskPath = currentT1MaskPath(config);
if isempty(maskPath)
    return;
end

mask = loadT1MaskFile(maskPath);
mask = fftshift(logical(mask), 3);
targetSize = size(t1Stage.t1MapMs);
targetSize = targetSize(1:3);
if ~isequal(size(mask), targetSize)
    warning('run_lrt_ddce_pipeline:T1MaskSizeMismatch', ...
        'Current T1 mask size %s does not match loaded T1 maps size %s; using mask saved inside t1_maps.mat.', ...
        mat2str(size(mask)), mat2str(targetSize));
    return;
end

t1Stage.t1Mask = mask;
fprintf('LRT dDCE pipeline: refreshed T1 mask from %s\n', maskPath);

end

function maskPath = currentT1MaskPath(config)

if ~isempty(config.t1MaskPath) && isfile(config.t1MaskPath)
    maskPath = config.t1MaskPath;
    return;
end

candidate = fullfile(config.outputDir, 'masks', 't1_mask.mat');
if isfile(candidate)
    maskPath = candidate;
else
    maskPath = '';
end

end

function mask = loadT1MaskFile(maskPath)

loaded = load(maskPath);
if isfield(loaded, 'mask')
    mask = loaded.mask;
elseif isfield(loaded, 't1Mask')
    mask = loaded.t1Mask;
else
    error('run_lrt_ddce_pipeline:InvalidT1Mask', ...
        'T1 mask file must contain mask or t1Mask: %s', maskPath);
end

end

function [dceInput, dceStage] = loadDceStage(config)

gdPath = requireRestartFile(config, 'gd_concentration.mat', 'synthesis');
aifPath = requireRestartFile(config, 'aif_fit.mat', 'synthesis');
maps2Path = requireRestartFile(config, 'ddce_2cxm_maps.mat', 'synthesis');

gd = load(gdPath);
aif = load(aifPath);
maps2 = load(maps2Path);
assertFields(gd, {'dceInput', 'Gdcon', 'dR1', 'precontrastT1Ms', ...
    'roiMasksBySlice', 'sliceResults'}, 'gd_concentration.mat');
assertFields(aif, {'aifBySlice', 'curvesBySlice'}, 'aif_fit.mat');
assertFields(maps2, {'maps2CXM'}, 'ddce_2cxm_maps.mat');

dceInput = gd.dceInput;
dceStage = struct();
dceStage.Gdcon = gd.Gdcon;
dceStage.dR1 = gd.dR1;
dceStage.precontrastT1Ms = gd.precontrastT1Ms;
dceStage.dceSliceIndices = loadedDceSliceIndices(gd);
dceStage.sliceResults = gd.sliceResults;
dceStage.aifBySlice = aif.aifBySlice;
dceStage.curvesBySlice = aif.curvesBySlice;
dceStage.roiMasksBySlice = gd.roiMasksBySlice;
dceStage.maps2CXM = maps2.maps2CXM;

mapsEPath = fullfile(config.outputDir, 'ddce_etk_maps.mat');
if isfile(mapsEPath)
    mapsE = load(mapsEPath);
    if isfield(mapsE, 'mapsETK')
        dceStage.mapsETK = mapsE.mapsETK;
    end
end
dceStage = filterLoadedDceStageSlices(dceStage, config);

end

function dceSliceIndices = loadedDceSliceIndices(gd)

if isfield(gd, 'dceSliceIndices') && ~isempty(gd.dceSliceIndices)
    dceSliceIndices = gd.dceSliceIndices(:).';
elseif ~isempty(gd.roiMasksBySlice) && isfield(gd.roiMasksBySlice, 'sliceIndex')
    dceSliceIndices = [gd.roiMasksBySlice.sliceIndex];
else
    error('run_lrt_ddce_pipeline:MissingDceSliceIndices', ...
        'Could not infer dceSliceIndices from gd_concentration.mat.');
end

end

function dceStage = filterLoadedDceStageSlices(dceStage, config)

if isempty(config.dceSliceIndices)
    return;
end

requested = validateIndices(config.dceSliceIndices(:).', size(dceStage.Gdcon, 3), 'dceSliceIndices');
[isAvailable, positions] = ismember(requested, dceStage.dceSliceIndices);
if any(~isAvailable)
    error('run_lrt_ddce_pipeline:MissingRestartSlice', ...
        'Starting at ''synthesis'' requested DCE slices %s, but saved DCE artifacts only contain slices %s.', ...
        mat2str(requested(~isAvailable)), mat2str(dceStage.dceSliceIndices));
end

dceStage.dceSliceIndices = requested;
dceStage.sliceResults = dceStage.sliceResults(positions);
dceStage.aifBySlice = dceStage.aifBySlice(positions);
dceStage.curvesBySlice = dceStage.curvesBySlice(positions);
dceStage.roiMasksBySlice = dceStage.roiMasksBySlice(positions);

end

function pathIn = requireRestartFile(config, fileName, startStage)

pathIn = fullfile(config.outputDir, fileName);
if ~isfile(pathIn)
    error('run_lrt_ddce_pipeline:MissingRestartFile', ...
        'Starting at ''%s'' requires %s in config.outputDir: %s', ...
        startStage, fileName, pathIn);
end

end

function config = validatePipelineConfig(config)

config.startStage = lower(char(config.startStage));
if ~any(strcmp(config.startStage, {'t1', 'dce', 'synthesis'}))
    error('run_lrt_ddce_pipeline:InvalidStartStage', ...
        'config.startStage must be ''t1'', ''dce'', or ''synthesis''.');
end
config.models = validateDceModels(config.models);
if ~(islogical(config.dceOptionalRois) || isnumeric(config.dceOptionalRois)) || ...
        ~isscalar(config.dceOptionalRois)
    error('run_lrt_ddce_pipeline:InvalidDceOptionalRois', ...
        'config.dceOptionalRois must be true or false.');
end
config.dceOptionalRois = logical(config.dceOptionalRois);
if isempty(config.reconMatPath) || ~isfile(config.reconMatPath)
    error('run_lrt_ddce_pipeline:MissingReconMat', ...
        'config.reconMatPath must point to an existing MAT file.');
end
if isempty(config.outputDir)
    [reconDir, reconName] = fileparts(config.reconMatPath);
    config.outputDir = fullfile(reconDir, [reconName '_lrt_ddce_pipeline']);
end
config.exportSynthDicom = validateScalarLogical(config.exportSynthDicom, 'exportSynthDicom');
if config.exportSynthDicom && isempty(config.synthDicomOutputDir)
    synthDicomFid = extractFidForSynthDicom(config);
    config.synthDicomOutputDir = fullfile('/Volumes/Extreme SSD/Anzhen/AMI/DICOM_converted', ...
        sprintf('FID%s', synthDicomFid));
end
config.synthDicomScaleMax = validatePositiveScalar(config.synthDicomScaleMax, 'synthDicomScaleMax');
config.synthDicomWindowCenter = validatePositiveScalar(config.synthDicomWindowCenter, 'synthDicomWindowCenter');
config.synthDicomWindowWidth = validatePositiveScalar(config.synthDicomWindowWidth, 'synthDicomWindowWidth');
if config.exportSynthDicom
    if ~exist(config.synthDicomRoot, 'dir')
        error('run_lrt_ddce_pipeline:MissingSynthDicomRoot', ...
            'config.synthDicomRoot does not exist: %s', config.synthDicomRoot);
    end
    if ~isfile(config.synthDicomLookupFile)
        error('run_lrt_ddce_pipeline:MissingSynthDicomLookupFile', ...
            'config.synthDicomLookupFile does not exist: %s', config.synthDicomLookupFile);
    end
end
if isempty(config.timeMinutes)
    error('run_lrt_ddce_pipeline:MissingTimeMinutes', ...
        'config.timeMinutes must contain postcontrast acquisition times in minutes.');
end
if isempty(config.fitPoints)
    config.fitPoints = 1:numel(config.timeMinutes);
end
config.fitPoints = validateIndices(config.fitPoints(:).', numel(config.timeMinutes), 'fitPoints');
if ~isnumeric(config.precontrastT1Ms) || any(config.precontrastT1Ms(:) <= 0)
    error('run_lrt_ddce_pipeline:InvalidPreT1', ...
        'config.precontrastT1Ms must be positive when a precontrast map is not used.');
end
if ~isempty(config.precontrastT1MapPath) && ~isfile(config.precontrastT1MapPath)
    error('run_lrt_ddce_pipeline:MissingPreT1Map', ...
        'config.precontrastT1MapPath was provided but does not exist.');
end

end

function models = validateDceModels(models)

if ischar(models)
    models = {models};
elseif supportsStringArrays() && isstring(models)
    models = cellstr(models);
end
if ~iscell(models) || isempty(models)
    error('run_lrt_ddce_pipeline:InvalidModels', ...
        'config.models must be ''2CXM'', ''ETK'', or a cell array containing one or both.');
end

validated = {};
for n = 1:numel(models)
    if ~(ischar(models{n}) || (supportsStringArrays() && isstring(models{n})))
        error('run_lrt_ddce_pipeline:InvalidModels', ...
            'config.models entries must be ''2CXM'' or ''ETK''.');
    end
    name = lower(char(models{n}));
    switch name
        case '2cxm'
            canonical = '2CXM';
        case 'etk'
            canonical = 'ETK';
        otherwise
            error('run_lrt_ddce_pipeline:InvalidModels', ...
                'Unknown DCE model: %s. Use ''2CXM'', ''ETK'', or both.', char(models{n}));
    end
    if ~any(strcmp(validated, canonical))
        validated{end+1} = canonical; %#ok<AGROW>
    end
end

models = validated;

end

function tf = supportsStringArrays()

tf = exist('isstring', 'builtin') || exist('isstring', 'file');

end

function value = validateScalarLogical(value, fieldName)

if ~(islogical(value) || isnumeric(value)) || ~isscalar(value)
    error('run_lrt_ddce_pipeline:InvalidConfigValue', ...
        'config.%s must be true or false.', fieldName);
end
value = logical(value);

end

function value = validatePositiveScalar(value, fieldName)

if ~isnumeric(value) || ~isscalar(value) || value <= 0
    error('run_lrt_ddce_pipeline:InvalidConfigValue', ...
        'config.%s must be a positive scalar.', fieldName);
end

end

function sequence = getSequencePreset(name, reconParams, overrides)

if nargin < 3 || isempty(overrides)
    overrides = struct();
end

switch lower(name)
    case {'standard192', 'lrt192', 'seg192', 'standard', 'lrt'}
        sequence = struct('name', name, 'isVTR', false, 'Nseg', 192, ...
            'TR', 0.0114 + 0.0035, 'TI', 0.0105, 'flipAngleDeg', 5, ...
            'cutoff', 21, 'reps', 20, 'minT1', 0.1, 'maxT1', 2, ...
            'numT1', 401, 'reconNEco', getfieldwithdefault(reconParams, 'NEco', 1), ...
            'dictOptions', struct(), 'params', struct(), 'fitParams', struct());
    case {'seg500_sg4', 'seg500', 'sg4'}
        sequence = struct('name', name, 'isVTR', false, 'Nseg', 500, ...
            'TR', 0.00357, 'TI', 0.0105, 'flipAngleDeg', 5, ...
            'cutoff', 21, 'reps', 20, 'minT1', 0.1, 'maxT1', 2, ...
            'numT1', 401, 'reconNEco', getfieldwithdefault(reconParams, 'NEco', 1), ...
            'dictOptions', struct(), 'params', struct(), 'fitParams', struct());
    case {'vtr192', 'vtr'}
        sequence = struct('name', name, 'isVTR', true, 'Nseg', 192, ...
            'cutoff', 10, 'reconNEco', 1, 'dictOptions', struct(), ...
            'params', struct(), 'fitParams', struct());
    otherwise
        vtrMatch = regexp(lower(name), '^vtr(\d+)$', 'tokens', 'once');
        segMatch = regexp(lower(name), '^(?:seg|lrt|standard)(\d+)$', 'tokens', 'once');
        if ~isempty(vtrMatch)
            sequence = struct('name', name, 'isVTR', true, 'Nseg', str2double(vtrMatch{1}), ...
                'cutoff', 10, 'reconNEco', 1, 'dictOptions', struct(), ...
                'params', struct(), 'fitParams', struct());
        elseif ~isempty(segMatch)
            sequence = struct('name', name, 'isVTR', false, 'Nseg', str2double(segMatch{1}), ...
                'TR', 0.0114 + 0.0035, 'TI', 0.0105, 'flipAngleDeg', 5, ...
                'cutoff', 21, 'reps', 20, 'minT1', 0.1, 'maxT1', 2, ...
                'numT1', 401, 'reconNEco', getfieldwithdefault(reconParams, 'NEco', 1), ...
                'dictOptions', struct(), 'params', struct(), 'fitParams', struct());
        else
            error('run_lrt_ddce_pipeline:UnknownPreset', ...
                'Unknown sequencePreset: %s', name);
        end
end

sequence = applySequenceOverrides(sequence, overrides);
sequence = normalizeSequence(sequence);

end

function dict = buildT1Dictionary(sequence, reconParams)

if sequence.isVTR
    [paramsDefault, fitParams, dictOptions] = defaultT1DictionaryVTRInputs();
    paramsDefault = mergeStructs(paramsDefault, reconParams);
    paramsDefault = mergeStructs(paramsDefault, sequence.params);
    fitParams = mergeStructs(fitParams, sequence.fitParams);
    dictOptions = mergeStructs(dictOptions, sequence.dictOptions);
    paramsDefault.linesPerShot = sequence.Nseg;
    paramsDefault.Nseg = sequence.Nseg;
    paramsDefault.Necho = sequence.reconNEco;
    paramsDefault.NEco = sequence.reconNEco;
    dictOptions.cutoff = sequence.cutoff;
    dictRaw = genT1Dictionary_VTR(paramsDefault, fitParams, dictOptions);
    dict = struct();
    dict.T1sSec = dictRaw.T1s(:).';
    dict.T1sMs = dict.T1sSec * 1000;
    dict.curves = dictRaw.Mz_dict_norm_abs_truc.';
    dict.raw = dictRaw;
else
    alpha = sequence.flipAngleDeg * pi / 180;
    E1 = @(t, R1) exp(-t * R1);
    T1sSec = logspace(log10(sequence.minT1), log10(sequence.maxT1), sequence.numT1);
    R1s = 1 ./ T1sSec;
    curves = zeros(numel(T1sSec), sequence.Nseg);
    for k = 1:numel(R1s)
        R1 = R1s(k);
        M0 = 1;
        M00 = -M0;
        Mz = zeros(sequence.Nseg * sequence.reps, 1);
        for r = 1:sequence.reps
            M01 = M00 * E1(sequence.TI, R1) + M0 * (1 - E1(sequence.TI, R1));
            for i = 1:sequence.Nseg
                Mz(i + (r-1) * sequence.Nseg) = M01 * cos(alpha) * E1(sequence.TR, R1) + M0 * (1 - E1(sequence.TR, R1));
                M01 = Mz(i + (r-1) * sequence.Nseg);
            end
            M00 = -M01;
        end
        train = Mz(end-sequence.Nseg+1:end).';
        curves(k,:) = abs(train ./ max(train));
    end
    dict = struct();
    dict.T1sSec = T1sSec;
    dict.T1sMs = T1sSec * 1000;
    dict.curves = curves(:, sequence.cutoff+1:end);
end

end

function sequence = applySequenceOverrides(sequence, overrides)

if ~isstruct(overrides)
    error('run_lrt_ddce_pipeline:InvalidSequenceOverride', ...
        'config.sequence must be a struct.');
end

if hasOverride(overrides, 'linesPerShot')
    sequence.Nseg = overrides.linesPerShot;
end
if hasOverride(overrides, 'Nseg')
    sequence.Nseg = overrides.Nseg;
end
if hasOverride(overrides, 'cutoff')
    sequence.cutoff = overrides.cutoff;
end
if hasOverride(overrides, 'TR')
    sequence.TR = overrides.TR;
end
if hasOverride(overrides, 'TR1') || hasOverride(overrides, 'TR2')
    tr1 = getfieldwithdefault(overrides, 'TR1', 0);
    tr2 = getfieldwithdefault(overrides, 'TR2', 0);
    if tr1 > 0 && tr2 > 0
        sequence.TR = tr1 + tr2;
        sequence.TR1 = tr1;
        sequence.TR2 = tr2;
    end
end
if hasOverride(overrides, 'TI')
    sequence.TI = overrides.TI;
end
if hasOverride(overrides, 'flipAngleDeg')
    sequence.flipAngleDeg = overrides.flipAngleDeg;
end
if hasOverride(overrides, 'reps')
    sequence.reps = overrides.reps;
end
if hasOverride(overrides, 'minT1')
    sequence.minT1 = overrides.minT1;
end
if hasOverride(overrides, 'maxT1')
    sequence.maxT1 = overrides.maxT1;
end
if hasOverride(overrides, 'numT1')
    sequence.numT1 = overrides.numT1;
end
if hasOverride(overrides, 'reconNEco')
    sequence.reconNEco = overrides.reconNEco;
end
if isfield(overrides, 'dictOptions') && isstruct(overrides.dictOptions)
    sequence.dictOptions = mergeStructs(sequence.dictOptions, overrides.dictOptions);
end
if isfield(overrides, 'params') && isstruct(overrides.params)
    sequence.params = mergeStructs(sequence.params, overrides.params);
end
if isfield(overrides, 'fitParams') && isstruct(overrides.fitParams)
    sequence.fitParams = mergeStructs(sequence.fitParams, overrides.fitParams);
end

end

function tf = hasOverride(s, fieldName)

tf = isfield(s, fieldName) && ~isempty(s.(fieldName));

end

function sequence = normalizeSequence(sequence)

numericFields = {'Nseg', 'cutoff', 'reconNEco'};
for k = 1:numel(numericFields)
    fieldName = numericFields{k};
    if ~isfield(sequence, fieldName) || isempty(sequence.(fieldName)) || ...
            ~isnumeric(sequence.(fieldName)) || sequence.(fieldName) <= 0 || ...
            mod(sequence.(fieldName), 1) ~= 0
        error('run_lrt_ddce_pipeline:InvalidSequenceField', ...
            'Sequence field %s must be a positive integer.', fieldName);
    end
end
if sequence.cutoff >= sequence.Nseg
    error('run_lrt_ddce_pipeline:InvalidCutoff', ...
        'Sequence cutoff (%d) must be smaller than linesPerShot/Nseg (%d).', ...
        sequence.cutoff, sequence.Nseg);
end

if ~sequence.isVTR
    requiredPositive = {'TR', 'TI', 'flipAngleDeg', 'reps', 'minT1', 'maxT1', 'numT1'};
    for k = 1:numel(requiredPositive)
        fieldName = requiredPositive{k};
        if ~isfield(sequence, fieldName) || isempty(sequence.(fieldName)) || ...
                ~isnumeric(sequence.(fieldName)) || sequence.(fieldName) <= 0
            error('run_lrt_ddce_pipeline:InvalidSequenceField', ...
                'Sequence field %s must be positive.', fieldName);
        end
    end
    if sequence.minT1 >= sequence.maxT1
        error('run_lrt_ddce_pipeline:InvalidT1Grid', ...
            'sequence.minT1 must be smaller than sequence.maxT1.');
    end
    if mod(sequence.reps, 1) ~= 0 || mod(sequence.numT1, 1) ~= 0
        error('run_lrt_ddce_pipeline:InvalidSequenceField', ...
            'sequence.reps and sequence.numT1 must be positive integers.');
    end
end

end

function fitTrain = cropRecoveryTrain(recon, sequence)

if size(recon, 3) < sequence.Nseg
    error('run_lrt_ddce_pipeline:ShortRecoveryTrain', ...
        'Reconstructed recovery train has %d frames, but preset %s expects %d.', ...
        size(recon, 3), sequence.name, sequence.Nseg);
end
fitTrain = recon(:,:,sequence.cutoff+1:sequence.Nseg);

end

function t1Map2d = fitT1Dictionary(ipt3d, mask1d, dict, Ny, Nx)

ipt2d = reshape(ipt3d, [], size(ipt3d, 3));
t1Map1d = zeros(1, size(ipt2d, 1));
for ii = 1:size(ipt2d, 1)
    if mask1d(ii)
        ipt = ipt2d(ii,:);
        denom = max(ipt);
        if denom > 0
            iptNorm = ipt ./ denom;
            diffs = dict.curves - repmat(iptNorm, size(dict.curves, 1), 1);
            [~, bestIdx] = min(sum(diffs.^2, 2));
            t1Map1d(ii) = dict.T1sMs(bestIdx);
        end
    end
end
t1Map2d = reshape(t1Map1d, Ny, Nx);

end

function mask = loadOrCreateT1Mask(data, sliceIndices, card, resp, sequence, config)

maskPath = config.t1MaskPath;
if isempty(maskPath)
    reconDir = fileparts(config.reconMatPath);
    candidates = {fullfile(config.outputDir, 'masks', 't1_mask.mat'), ...
        fullfile(reconDir, 'mask_rect.mat')};
    for c = 1:numel(candidates)
        if isfile(candidates{c})
            maskPath = candidates{c};
            break;
        end
    end
end

if ~isempty(maskPath) && isfile(maskPath)
    loaded = load(maskPath);
    if isfield(loaded, 'mask')
        mask = loaded.mask;
    elseif isfield(loaded, 't1Mask')
        mask = loaded.t1Mask;
    else
        error('run_lrt_ddce_pipeline:InvalidT1Mask', ...
            'T1 mask file must contain mask or t1Mask: %s', maskPath);
    end
    return;
end

if ~strcmpi(config.maskMode, 'reuse_or_draw')
    error('run_lrt_ddce_pipeline:MissingT1Mask', ...
        'No T1 mask found and maskMode is %s.', config.maskMode);
end

mask = false(data.Ny, data.Nx, data.Nz);
Uimg = reshape(data.U, data.Ny, data.Nx, data.Nz, []);
for slc = sliceIndices
    dispim = @(x) fftshift(x(:,:,slc,:), 1);
    sliceBasis = reshape(dispim(Uimg), [], data.L);
    coef = data.Gr \ reshape(data.Phi(:,:,card,resp,1), data.L, []);
    recon = reshape(sliceBasis * coef, data.Ny, data.Nx, [], sequence.reconNEco);
    ref = abs(recon(:,:,min(size(recon, 3), sequence.Nseg),1));
    figure('Name', sprintf('Draw T1 fitting mask, slice %d', slc));
    imagesc(ref); axis image off; colormap gray;
    title(sprintf('Draw T1 fitting mask, slice %d', slc));
    roi = drawpolygon;
    mask(:,:,slc) = createMask(roi);
    close(gcf);
end

t1Mask = mask;
save(fullfile(config.outputDir, 'masks', 't1_mask.mat'), 'mask', 't1Mask');

end

function roiMasks = loadOrCreateDceRois(GdSlice, slc, config)

roiPath = fullfile(config.outputDir, 'masks', sprintf('MASK_ROI_Slice%d.mat', slc));
if isfile(roiPath)
    loaded = load(roiPath);
    needed = {'Blood', 'HeartMask'};
    assertFields(loaded, needed, roiPath);
    roiMasks = normalizeDceRoiMasks(loaded);
    return;
end

if ~strcmpi(config.maskMode, 'reuse_or_draw')
    error('run_lrt_ddce_pipeline:MissingDceRoi', ...
        'ROI file not found for slice %d and maskMode is %s.', slc, config.maskMode);
end

refFrame = max(1, size(GdSlice, 3) - 2);
figure('Name', sprintf('Draw dDCE ROIs, slice %d', slc));
imagesc(GdSlice(:,:,refFrame)); axis image off; colormap gray; climAuto();
title(sprintf('Slice %d: draw blood pool', slc));
Blood = roipoly;
title(sprintf('Slice %d: draw heart mask', slc));
HeartMask = roipoly;
if config.dceOptionalRois
    close(gcf);
    save(roiPath, 'Blood', 'HeartMask');
    roiMasks = struct('Blood', logical(Blood), 'HeartMask', logical(HeartMask));
    return;
end

title(sprintf('Slice %d: draw LV endocardium (optional)', slc));
LVendo = roipoly;
title(sprintf('Slice %d: draw LV epicardium (optional)', slc));
LVepi = roipoly;
Myomask = ~LVendo & LVepi;
title(sprintf('Slice %d: draw MI mask (optional)', slc));
MIMask = roipoly;
title(sprintf('Slice %d: draw remote mask (optional)', slc));
RemoteMask = roipoly;
close(gcf);

save(roiPath, 'Blood', 'HeartMask', 'LVendo', 'LVepi', 'Myomask', 'MIMask', 'RemoteMask');
roiMasks = struct('Blood', logical(Blood), 'HeartMask', logical(HeartMask), ...
    'LVendo', logical(LVendo), 'LVepi', logical(LVepi), 'Myomask', logical(Myomask), ...
    'MIMask', logical(MIMask), 'RemoteMask', logical(RemoteMask));

end

function roiMasks = normalizeDceRoiMasks(loaded)

roiMasks = struct();
roiMasks.Blood = logical(loaded.Blood);
roiMasks.HeartMask = logical(loaded.HeartMask);
optionalNames = {'LVendo', 'LVepi', 'Myomask', 'MIMask', 'RemoteMask'};
for n = 1:numel(optionalNames)
    name = optionalNames{n};
    if isfield(loaded, name)
        roiMasks.(name) = logical(loaded.(name));
    end
end
if ~isfield(roiMasks, 'Myomask') && isfield(roiMasks, 'LVendo') && isfield(roiMasks, 'LVepi')
    roiMasks.Myomask = ~roiMasks.LVendo & roiMasks.LVepi;
end

end

function qualityMask = dceFitQualityMask(roiMasks)

if isfield(roiMasks, 'Myomask') && any(roiMasks.Myomask(:))
    qualityMask = roiMasks.Myomask;
else
    qualityMask = roiMasks.HeartMask;
end

end

function [Ctoi, Cp, time, fitPointInd] = buildDceCurves(GdSlice, roiMasks, config)

fitPointInd = config.fitPoints(:).';
fitPointInd = fitPointInd(fitPointInd <= size(GdSlice, 3));
timeIn = config.timeMinutes(:).';
time = timeIn(fitPointInd);
Cp = zeros(1, numel(fitPointInd));
Ctoi = zeros(sum(roiMasks.HeartMask(:)), numel(fitPointInd));
for n = 1:numel(fitPointInd)
    temp = GdSlice(:,:,fitPointInd(n));
    Cp(n) = meanFinite(temp(roiMasks.Blood));
    Ctoi(:,n) = temp(roiMasks.HeartMask);
end
Cp(~isfinite(Cp)) = 0;
Ctoi(~isfinite(Ctoi)) = 0;

end

function [aif, aifFit] = fitAif(Cp, time, config)

timeRound = round(time * 60);
timeInMin = timeRound / 60;
timeInterp = (0:1:max(timeRound)) / 60;
base = config.aifInitial;

initialParams = [base.B, base.m1, base.m2, base.tc, base.Sf];
lb = [0, 0, 0, 0, 1];
ub = [Inf, Inf, Inf, Inf, 1];
fitFunction = @(params, paramsCell) gaussianExpP1Only(struct( ...
    'B', params(1), 'm1', params(2), 'm2', params(3), ...
    'tc', params(4), 'Sf', params(5)), paramsCell);
options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);

if numel(Cp) >= 3
    CpGauss = [0 Cp(3:end)];
    timeGauss = [0 timeInMin(3:end)];
else
    CpGauss = [0 Cp(:).'];
    timeGauss = [0 timeInMin(:).'];
end
paramsCell = {base.a1, base.sigma1, base.mu1, base.a2, base.sigma2, base.mu2, timeGauss};
fittedParamsGauss = lsqcurvefit(fitFunction, initialParams, paramsCell, CpGauss, lb, ub, options);
fittedParamsStructGauss = struct('B', fittedParamsGauss(1), 'm1', fittedParamsGauss(2), ...
    'm2', fittedParamsGauss(3), 'tc', fittedParamsGauss(4), 'Sf', fittedParamsGauss(5));

initialParams = [base.a1, base.a2, base.sigma1, base.sigma2];
lb = [0, 0, 0, 0];
ub = [Inf, Inf, Inf, Inf];
fitFunction = @(params, paramsCell) gaussianP23OnlyFixedMu(struct( ...
    'a1', params(1), 'a2', params(2), 'sigma1', params(3), ...
    'sigma2', params(4)), paramsCell);
CpExp = Cp(1:min(3, numel(Cp)));
timeExp = timeInMin(1:min(3, numel(timeInMin)));
paramsCell = {fittedParamsGauss(1), fittedParamsGauss(2), fittedParamsGauss(3), ...
    fittedParamsGauss(4), fittedParamsGauss(5), base.mu1, base.mu2, timeExp};
fittedParamsExp = lsqcurvefit(fitFunction, initialParams, paramsCell, CpExp, lb, ub, options);

fittedParamsStruct = struct('B', base.B, 'a1', base.a1, 'a2', base.a2, ...
    'm1', base.m1, 'm2', base.m2, 'sigma1', base.sigma1, 'sigma2', base.sigma2, ...
    'mu1', base.mu1, 'mu2', base.mu2, 'tc', base.tc);
fittedParamsStructExp = struct('a1', fittedParamsExp(1), 'a2', fittedParamsExp(2), ...
    'sigma1', fittedParamsExp(3), 'sigma2', fittedParamsExp(4));
paramsCellInterp = {fittedParamsGauss(1), fittedParamsGauss(2), fittedParamsGauss(3), ...
    fittedParamsGauss(4), fittedParamsGauss(5), base.mu1, base.mu2, timeInterp};
aif = gaussianP23OnlyFixedMu(fittedParamsStructExp, paramsCellInterp);
aif(1) = 0;

aifFit = struct();
aifFit.timeInputMinutes = timeInMin;
aifFit.Cp = Cp;
aifFit.timeInterpMinutes = timeInterp;
aifFit.aif = aif;
aifFit.fittedParams_gauss = fittedParamsGauss;
aifFit.fittedParams_exp = fittedParamsExp;
aifFit.fittedParamsStruct = fittedParamsStruct;
aifFit.fittedParamsStruct_gauss = fittedParamsStructGauss;
aifFit.fittedParamsStruct_exp = fittedParamsStructExp;

end

function fit = fit2CxmSlice(Ctoi, aif, time, HeartMask, Myomask, config)

settings = config.fit2CXM;
numVox = size(Ctoi, 1);
fitVector = zeros(numVox, 6);
Rsq = -ones(numVox, 1);
x0c = repmat(settings.x0, numVox, 1);
timeInMin = round(time * 60) / 60;

for iter = 1:settings.maxIterations
    for c = 1:numVox
        if Rsq(c) < settings.targetR2
            x0 = x0c(c,:);
            if iter > 1
                x0 = x0 .* (1 - (rand(1, numel(x0)) - 0.5) * 1.5);
                x0 = min(max(x0, settings.lb), settings.ub);
            end
            tempFit = fitdcemri_etk_fft(Ctoi(c,:)', aif', timeInMin', x0, settings.lb, settings.ub, 'cxm_bound');
            if tempFit(end) > Rsq(c)
                fitVector(c,1:numel(tempFit)) = tempFit(:).';
                x0c(c,:) = tempFit(1:4);
                Rsq(c) = tempFit(end);
            end
        end
    end
    if qualityReached(HeartMask, Myomask, fitVector(:,end), settings.targetR2)
        break;
    end
end

fit = vectorTo2CxmMaps(fitVector, HeartMask);

end

function fit = fitEtkSlice(Ctoi, aif, time, HeartMask, Myomask, config)

settings = config.fitETK;
numVox = size(Ctoi, 1);
fitVector = zeros(numVox, 5);
Rsq = -ones(numVox, 1);
x0c = repmat(settings.x0, numVox, 1);
x0temp = x0c;
timeInMin = round(time * 60) / 60;
omitIdx = settings.omitTimeIndex;
if omitIdx > 0 && omitIdx <= numel(timeInMin)
    timeFit = timeInMin;
    timeFit(omitIdx) = [];
else
    timeFit = timeInMin;
    omitIdx = [];
end

for iter = 1:settings.maxIterations
    for c = 1:numVox
        if Rsq(c) < settings.targetR2
            sigin = Ctoi(c,:);
            if ~isempty(omitIdx) && omitIdx <= numel(sigin)
                sigin(omitIdx) = [];
            end
            tempFit = fitdcemri_etk_resample(sigin', aif', timeFit', x0temp(c,:), settings.lb, settings.ub, 'etk');
            if tempFit(end) > Rsq(c)
                fitVector(c,1:numel(tempFit)) = tempFit(:).';
                Rsq(c) = tempFit(end);
                x0c(c,:) = x0temp(c,:);
            end
            x0temp(c,:) = max(settings.lb, x0c(c,:) .* (1 - (rand(1,3) - 0.5) * 1.5));
            x0temp(c,:) = min(settings.ub, x0temp(c,:));
        end
    end
    if qualityReached(HeartMask, Myomask, fitVector(:,end), settings.targetR2)
        break;
    end
end

fit = vectorToEtkMaps(fitVector, HeartMask);

end

function ok = qualityReached(HeartMask, Myomask, fitQVector, targetR2)

fitQ = zeros(size(HeartMask));
fitQ(HeartMask) = fitQVector;
myoQ = fitQ(Myomask & fitQ > 0);
if isempty(myoQ)
    ok = false;
    return;
end
ok = prctile(myoQ, 5) >= targetR2;

end

function fit = vectorTo2CxmMaps(fitVector, HeartMask)

fit = blank2dMaps(size(HeartMask), {'F', 'PS', 'vp', 've', 'RMSE', 'Rsq'});
names = {'F', 'PS', 'vp', 've', 'RMSE', 'Rsq'};
for n = 1:numel(names)
    temp = zeros(size(HeartMask));
    temp(HeartMask) = fitVector(:,n);
    fit.(names{n}) = temp;
end

end

function fit = vectorToEtkMaps(fitVector, HeartMask)

fit = blank2dMaps(size(HeartMask), {'Ktrans', 'Kep', 'Vp', 'RMSE', 'Rsq'});
names = {'Ktrans', 'Kep', 'Vp', 'RMSE', 'Rsq'};
for n = 1:numel(names)
    temp = zeros(size(HeartMask));
    temp(HeartMask) = fitVector(:,n);
    fit.(names{n}) = temp;
end
fit.Ve = zeros(size(HeartMask));
valid = fit.Kep ~= 0;
fit.Ve(valid) = fit.Ktrans(valid) ./ fit.Kep(valid);
fit.Ve = fit.Ve .* HeartMask;

end

function curves = summarizeDceCurves(GdSlice, roiMasks, config)

names = {'Blood', 'Myomask', 'MIMask', 'RemoteMask'};
outNames = {'blood', 'myocardium', 'mi', 'remote'};
curves = struct();
curves.timeMinutes = config.timeMinutes(1:size(GdSlice, 3));
for n = 1:numel(names)
    if ~isfield(roiMasks, names{n}) || ~any(roiMasks.(names{n})(:))
        continue;
    end
    vals = zeros(1, size(GdSlice, 3));
    mask = roiMasks.(names{n});
    for t = 1:size(GdSlice, 3)
        img = GdSlice(:,:,t);
        vals(t) = meanFinite(img(mask));
    end
    curves.(outNames{n}) = vals;
end

end

function preT1 = loadPrecontrastT1(t1ForDce, config)

if ~isempty(config.precontrastT1MapPath)
    loaded = load(config.precontrastT1MapPath);
    fieldName = findFirstNumericField(loaded);
    preT1 = loaded.(fieldName);
    preT1 = applyOrientationTransform(preT1, config.orientationTransform);
    if ndims(preT1) == 4
        preT1 = preT1(:,:,:,1);
    elseif ndims(preT1) == 5
        preT1 = squeeze(mean(preT1(:,:,:,1,:), 5));
    end
    if ~isequal(size(preT1), size(t1ForDce(:,:,:,1)))
        error('run_lrt_ddce_pipeline:PreT1SizeMismatch', ...
            'Precontrast T1 map size %s does not match postcontrast map size %s.', ...
            mat2str(size(preT1)), mat2str(size(t1ForDce(:,:,:,1))));
    end
else
    preT1 = ones(size(t1ForDce, 1), size(t1ForDce, 2), size(t1ForDce, 3)) * config.precontrastT1Ms;
end

end

function arr = applyOrientationTransform(arr, transforms)

for k = 1:numel(transforms)
    switch lower(transforms{k})
        case 'rot90ccw'
            arr = rot90FirstTwoDims(arr, 1);
        case 'rot90cw'
            arr = rot90FirstTwoDims(arr, -1);
        case 'flip_lr'
            arr = flip(arr, 2);
        case 'flip_ud'
            arr = flip(arr, 1);
        case {'none', ''}
        otherwise
            error('run_lrt_ddce_pipeline:UnknownOrientationTransform', ...
                'Unknown orientation transform: %s', transforms{k});
    end
end

end

function arr = rot90FirstTwoDims(arr, direction)

if ndims(arr) < 3
    arr = rot90(arr, direction);
    return;
end

order = 1:ndims(arr);
order(1:2) = [2 1];
arr = permute(arr, order);
if direction > 0
    arr = flip(arr, 1);
else
    arr = flip(arr, 2);
end

end

function maps = initSliceMapStruct(volSize, names)

maps = struct();
for n = 1:numel(names)
    maps.(names{n}) = zeros(volSize(1), volSize(2), volSize(3));
end

end

function maps = blank2dMaps(mapSize, names)

maps = struct();
for n = 1:numel(names)
    maps.(names{n}) = zeros(mapSize);
end

end

function writePipelineQc(t1Stage, dceStage, synthStage, config)

qcDir = fullfile(config.outputDir, 'qc');
if isfield(t1Stage, 't1MapMs')
    writeMontagePng(squeeze(t1Stage.t1MapMs(:,:,:,end,1)), fullfile(qcDir, 't1_montage.png'), [0 1000]);
end
if ~isempty(dceStage.dceSliceIndices)
    for idx = 1:numel(dceStage.dceSliceIndices)
        slc = dceStage.dceSliceIndices(idx);
        writeGif(squeeze(dceStage.Gdcon(:,:,slc,:)), fullfile(qcDir, sprintf('gdcon_slice%d.gif', slc)), [0 2]);
        writeAifPlot(dceStage.aifBySlice(idx), fullfile(qcDir, sprintf('aif_slice%d.png', slc)));
    end
end
if isfield(dceStage, 'maps2CXM')
    writeParamMontage(dceStage.maps2CXM, fullfile(qcDir, 'maps_2cxm.png'), {'F', 'PS', 'vp', 've', 'Rsq'});
end
if isfield(dceStage, 'mapsETK')
    writeParamMontage(dceStage.mapsETK, fullfile(qcDir, 'maps_etk.png'), {'Ktrans', 'Kep', 'Ve', 'Vp', 'Rsq'});
end
if isfield(synthStage, 'synthBySlice') && ~isempty(synthStage.synthBySlice)
    for idx = 1:numel(synthStage.synthBySlice)
        slc = synthStage.synthBySlice(idx).sliceIndex;
        synthCell = synthStage.synthBySlice(idx).synthLGE;
        scaleFactors = synthStage.synthBySlice(idx).scaleFactors;
        for sfIdx = 1:numel(synthCell)
            synth = synthCell{sfIdx};
            suffix = synthesisQcSuffix(scaleFactors, sfIdx);
            writeMontagePng(synth, fullfile(qcDir, sprintf('synth_lge_slice%d%s.png', slc, suffix)), [0 2]);
            writeGif(synth, fullfile(qcDir, sprintf('synth_lge_slice%d%s.gif', slc, suffix)), [0 2]);
        end
    end
end

end

function suffix = synthesisQcSuffix(scaleFactors, sfIdx)

if numel(scaleFactors) <= 1
    suffix = '';
    return;
end
label = sprintf('%g', scaleFactors(sfIdx));
label = strrep(label, '-', 'm');
label = strrep(label, '.', 'p');
suffix = sprintf('_scale%s', label);

end

function writeMontagePng(vol, pathOut, displayRange)

fig = figure('Visible', 'off');
try
    montage(vol, 'DisplayRange', displayRange);
catch
    num = size(vol, 3);
    rows = ceil(sqrt(num));
    cols = ceil(num / rows);
    for i = 1:num
        subplot(rows, cols, i);
        imagesc(vol(:,:,i), displayRange); axis image off; colormap gray;
    end
end
saveas(fig, pathOut);
close(fig);

end

function writeParamMontage(maps, pathOut, names)

fig = figure('Visible', 'off', 'Position', [100 100 1200 800]);
for n = 1:numel(names)
    subplot(2, ceil(numel(names)/2), n);
    vol = maps.(names{n});
    imagesc(max(vol, [], 3)); axis image off; colorbar; title(names{n}, 'Interpreter', 'none');
end
saveas(fig, pathOut);
close(fig);

end

function writeGif(vol, pathOut, displayRange)

for idx = 1:size(vol, 3)
    frame = mat2gray(vol(:,:,idx), displayRange);
    [imind, cm] = gray2ind(frame, 256);
    if idx == 1
        imwrite(imind, cm, pathOut, 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, pathOut, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

end

function writeAifPlot(aifBySlice, pathOut)

fig = figure('Visible', 'off');
plot(aifBySlice.fit.timeInputMinutes, aifBySlice.fit.Cp, 'or', 'MarkerFaceColor', 'r');
hold on;
plot(aifBySlice.fit.timeInterpMinutes, aifBySlice.fit.aif, '-b', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Gd concentration (mM)');
legend({'Measured blood', 'Fitted AIF'}, 'Location', 'best');
grid on;
saveas(fig, pathOut);
close(fig);

end

function Cp = gaussianExpP1Only(obj, paramsCell)

a1 = paramsCell{1};
sigma1 = paramsCell{2};
mu1 = paramsCell{3};
a2 = paramsCell{4};
sigma2 = paramsCell{5};
mu2 = paramsCell{6};
time = paramsCell{7};
p1 = (obj.B .* exp(-obj.m1 .* time)) ./ (1 + exp(-obj.m2 .* (time - obj.tc)));
p2 = (a1 / sigma1 / sqrt(2 * pi)) .* exp(-((time - mu1) ./ sigma1 .* sqrt(2)).^2);
p3 = (a2 / sigma2 / sqrt(2 * pi)) .* exp(-((time - mu2) ./ sigma2 .* sqrt(2)).^2);
Cp = obj.Sf * (p1 + p2 + p3);

end

function Cp = gaussianP23OnlyFixedMu(obj, paramsCell)

B = paramsCell{1};
m1 = paramsCell{2};
m2 = paramsCell{3};
tc = paramsCell{4};
Sf = paramsCell{5};
mu1 = paramsCell{6};
mu2 = paramsCell{7};
time = paramsCell{8};
p1 = (B .* exp(-m1 .* time)) ./ (1 + exp(-m2 .* (time - tc)));
p2 = (obj.a1 / obj.sigma1 / sqrt(2 * pi)) .* exp(-((time - mu1) ./ obj.sigma1 .* sqrt(2)).^2);
p3 = (obj.a2 / obj.sigma2 / sqrt(2 * pi)) .* exp(-((time - mu2) ./ obj.sigma2 .* sqrt(2)).^2);
Cp = Sf * (p1 + p2 + p3);

end

function Ctoi = cxmBoundCurve(beta, X)

time = X{1};
Cp = X{2};
startPoint = X{3};
F = beta(1);
PS = beta(2);
vp = beta(3);
ve = beta(4);
a = (F + PS) / max(vp, eps);
b = PS / max(ve, eps);
c = F / max(vp, eps);
disc = max((a + b).^2 - 4 * b * c, 0);
M1 = 0.5 * (a + b + sqrt(disc));
M2 = 0.5 * (a + b - sqrt(disc));
B = (M2 - c) / (M2 - M1);
tFine = (time(1):1/60:time(end)).';
H = B * exp(-M1 * tFine) + (1 - B) * exp(-M2 * tFine);
Ctoi = F * customConvolutionResampled(H, Cp, time);
Ctoi = Ctoi(startPoint:end);

end

function y = customConvolutionResampled(H, Cp, t)

t = t(:);
H = H(:);
Cp = Cp(:);
dtFine = 1 / 60;
tFine = (t(1):dtFine:t(end)).';
if numel(Cp) ~= numel(tFine)
    Cp = interp1(linspace(tFine(1), tFine(end), numel(Cp)), Cp, tFine, 'linear', 'extrap');
end
if numel(H) ~= numel(tFine)
    H = interp1(linspace(tFine(1), tFine(end), numel(H)), H, tFine, 'linear', 'extrap');
end
N = numel(tFine);
L = 2^nextpow2(2 * N - 1);
convResult = real(ifft(fft(H, L) .* fft(Cp, L)));
convResult = convResult(1:N) * dtFine;
y = interp1(tFine, convResult, t, 'linear');

end

function dceSliceIndices = chooseDceSlicesFromMask(dceInput, numSlices)

if isfield(dceInput, 't1Mask') && ~isempty(dceInput.t1Mask)
    t1Mask = logical(dceInput.t1Mask);
    if size(t1Mask, 3) ~= numSlices
        error('run_lrt_ddce_pipeline:T1MaskSizeMismatch', ...
            'T1 mask has %d slices after orientation transform, but Gd concentration has %d slices.', ...
            size(t1Mask, 3), numSlices);
    end
    sliceHasMask = squeeze(any(any(t1Mask, 1), 2));
    dceSliceIndices = find(sliceHasMask).';
    if ~isempty(dceSliceIndices)
        fprintf('LRT dDCE pipeline: using nonzero T1 mask slices for DCE: %s\n', ...
            mat2str(dceSliceIndices));
        return;
    end
    warning('run_lrt_ddce_pipeline:EmptyT1Mask', ...
        'T1 mask has no nonzero slices after orientation transform; falling back to default DCE slice selection.');
else
    fprintf('LRT dDCE pipeline: no T1 mask available for DCE slice selection; using default DCE slice.\n');
end

dceSliceIndices = chooseDefaultDceSlices(dceInput.Gdcon);
fprintf('LRT dDCE pipeline: using default DCE slice: %s\n', mat2str(dceSliceIndices));

end

function dceSliceIndices = chooseDefaultDceSlices(Gdcon)
% max energy doesn't work
% try min energy
energy = squeeze(sum(sum(min(Gdcon, [], 4), 1), 2));
[~, best] = min(energy);
if isempty(best) || best < 1
    best = 1;
end
dceSliceIndices = best;

end

function idx = validateIndices(idx, maxIdx, fieldName)

idx = idx(:).';
if any(idx < 1) || any(idx > maxIdx) || any(mod(idx, 1) ~= 0)
    error('run_lrt_ddce_pipeline:InvalidIndex', ...
        '%s must contain integer indices from 1 to %d.', fieldName, maxIdx);
end

end

function value = meanFinite(values)

values = values(isfinite(values));
if isempty(values)
    value = 0;
else
    value = mean(values);
end

end

function out = mergeStructs(defaults, overrides)

out = defaults;
if isempty(overrides)
    return;
end
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    if isstruct(overrides.(name)) && isfield(out, name) && isstruct(out.(name))
        out.(name) = mergeStructs(out.(name), overrides.(name));
    else
        out.(name) = overrides.(name);
    end
end

end

function value = getfieldwithdefault(s, fieldName, defaultValue)

if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end

end

function assertFields(s, fields, sourceName)

for i = 1:numel(fields)
    if ~isfield(s, fields{i})
        error('run_lrt_ddce_pipeline:MissingField', ...
            '%s is missing required field %s.', sourceName, fields{i});
    end
end

end

function fieldName = findFirstNumericField(s)

names = fieldnames(s);
for i = 1:numel(names)
    if isnumeric(s.(names{i}))
        fieldName = names{i};
        return;
    end
end
error('run_lrt_ddce_pipeline:NoNumericPreT1', ...
    'Precontrast T1 map MAT file does not contain a numeric array.');

end

function ensureDir(pathIn)

if ~exist(pathIn, 'dir')
    mkdir(pathIn);
end

end

function climAuto()

try
    clim([0 2]);
catch
    caxis([0 2]);
end

end
