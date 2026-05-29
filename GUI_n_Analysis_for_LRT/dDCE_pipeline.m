cd('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/GUI_n_Analysis_for_LRT')

config = default_lrt_ddce_pipeline_config();

config.reconMatPath = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID111140/FID11114_MID01040_FID111140_10_8mm_USR28%_L8_results_Echo_1_2026_05_14_23_13.mat';
config.outputDir = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID111140';

% config.reconMatPath = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID151790/FID15179_MID02102_FID151790_10_8mm_USR28%_L8_results_Echo_1_2026_04_21_15_47.mat';
% config.outputDir = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID151790';

config.reconMatPath = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID174088/FID17408_MID01067_FID174088_10_8mm_USR28%_L16_results_Echo_1_2026_04_20_14_32.mat';
config.outputDir = '/Volumes/Extreme SSD/Anzhen/AMI/ProcessedData/FID174088';


config.sequencePreset = 'vtr';
config.sequence.linesPerShot = 240;   % change to your sequence
config.sequence.cutoff = 20;          % early frames to skip

config.cardiacPhases = 12;            % choose diastolic phase
config.respPhase = 1;
config.sliceIndices = [];             % [] = all slices for T1 fitting
config.dceSliceIndices = [];          % [] = auto-pick one DCE slice

config.precontrastT1Ms = 1350;        % fallback scalar precontrast T1
% Optional:
% config.precontrastT1MapPath = '/path/to/precontrast_t1_map.mat';


%config.startStage = 'synthesis';
% config.startStage = 'dce';
config.synthesisScaleFactors = [1];
config.dceOptionalRois = true;
config.models = {'2CXM'}; 
config.exportSynthDicom = true;
config.synthDicomOutputDir = '';

outputs = run_lrt_ddce_pipeline(config);