function test_lrt_ddce_pipeline_static()
% test_lrt_ddce_pipeline_static  Lightweight static checks for the pipeline API.
%
% This test does not require a real reconstruction MAT file.

config = default_lrt_ddce_pipeline_config();
required = {'reconMatPath', 'outputDir', 'dceCodeDir', 'sequencePreset', ...
    'sequence', 'cardiacPhases', 'respPhase', 'sliceIndices', 'maskMode', ...
    'precontrastT1Ms', 'precontrastT1MapPath', 'timeMinutes', 'fitPoints', ...
    'gdRelaxivity', 'models', 'synthesisEndMinutes', 'synthesisStepMinutes', ...
    'synthesisScaleFactors', 'startStage', 'dceOptionalRois', ...
    'exportSynthDicom', 'synthDicomRoot', 'synthDicomLookupFile', ...
    'synthDicomOutputDir', 'synthDicomLgeSeriesDescriptionContains', ...
    'synthDicomScaleMax', 'synthDicomWindowCenter', 'synthDicomWindowWidth', ...
    'saveQc'};

for k = 1:numel(required)
    assert(isfield(config, required{k}), 'Missing default config field: %s', required{k});
end
assert(strcmp(config.startStage, 't1'), 'Default startStage should be t1.');
assert(config.dceOptionalRois == false, 'Default dceOptionalRois should be false.');
assert(config.exportSynthDicom == false, 'Default exportSynthDicom should be false.');

config.startStage = 'bad_stage';
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:InvalidStartStage');
end
assert(didError, 'Invalid startStage should produce a clear InvalidStartStage error.');

config.reconMatPath = fullfile(tempdir, 'definitely_missing_lrt_recon.mat');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingReconMat');
end
assert(didError, 'Missing reconstruction MAT should produce a clear MissingReconMat error.');

config = default_lrt_ddce_pipeline_config();
config.models = {'2CXM'};
config.reconMatPath = fullfile(tempdir, 'definitely_missing_lrt_recon.mat');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingReconMat');
end
assert(didError, '2CXM-only model selection should validate cleanly before file loading.');

config = default_lrt_ddce_pipeline_config();
config.models = {'ETK'};
config.reconMatPath = fullfile(tempdir, 'definitely_missing_lrt_recon.mat');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingReconMat');
end
assert(didError, 'ETK-only model selection should validate cleanly before file loading.');

config = default_lrt_ddce_pipeline_config();
config.models = {'bad_model'};
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:InvalidModels');
end
assert(didError, 'Invalid DCE model should produce a clear InvalidModels error.');

config = default_lrt_ddce_pipeline_config();
config.exportSynthDicom = true;
config.synthDicomRoot = fullfile(tempdir, 'definitely_missing_synth_dicom_root');
config.reconMatPath = fullfile(tempdir, 'FID00001_dummy_recon_for_synth_dicom.mat');
dummy = 1;
save(config.reconMatPath, 'dummy');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingSynthDicomRoot');
end
assert(didError, 'Missing synthetic DICOM root should produce a clear MissingSynthDicomRoot error.');

config = default_lrt_ddce_pipeline_config();
config.exportSynthDicom = true;
config.synthDicomRoot = tempdir;
config.synthDicomLookupFile = fullfile(tempdir, 'definitely_missing_synth_dicom_lookup.xlsx');
config.reconMatPath = fullfile(tempdir, 'FID00001_dummy_recon_for_synth_dicom.mat');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingSynthDicomLookupFile');
end
assert(didError, 'Missing synthetic DICOM lookup should produce a clear MissingSynthDicomLookupFile error.');

config = default_lrt_ddce_pipeline_config();
config.startStage = 'synthesis';
config.dceCodeDir = tempdir;
config.outputDir = tempname;
mkdir(config.outputDir);
config.reconMatPath = fullfile(config.outputDir, 'dummy_recon.mat');
dummy = 1;
save(config.reconMatPath, 'dummy');
didError = false;
try
    run_lrt_ddce_pipeline(config);
catch ME
    didError = strcmp(ME.identifier, 'run_lrt_ddce_pipeline:MissingRestartFile');
end
assert(didError, 'Synthesis restart should report missing restart artifacts clearly.');

config = default_lrt_ddce_pipeline_config();
config.sequencePreset = 'vtr240';
config.sequence.linesPerShot = 240;
config.sequence.cutoff = 12;
assert(strcmp(config.sequencePreset, 'vtr240'));
assert(config.sequence.linesPerShot == 240);
assert(config.sequence.cutoff == 12);

fprintf('test_lrt_ddce_pipeline_static passed.\n');

end
