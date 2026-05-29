function [params, FitParams, dictOptions] = defaultT1DictionaryVTRInputs()
% defaultT1DictionaryVTRInputs  Minimal inputs for genT1Dictionary_VTR.
%
% These defaults are intended as a runnable template for VTR T1 dictionary
% generation after multitasking reconstruction. Replace the sequence timing
% values with the values from your reconstructed FitParams.params.
%
% Example:
%   [params, FitParams, dictOptions] = defaultT1DictionaryVTRInputs();
%   dict = genT1Dictionary_VTR(params, FitParams, dictOptions);

%% Required scan timing and readout fields
params.linesPerShot = 192;          % Nseg: number of readout lines per recovery train
params.lEchoSpacing = 11.4e-3;      % First VTR spacing, seconds
params.alTE_seconds = [1.6 5.1]*1e-3; % Echo times, seconds; used to compute TRnav
params.Necho = numel(params.alTE_seconds);
params.SGBlock = 2;                 % VTR model currently supports SGBlock == 2
params.isVTR = 1;

%% Flip angle schedule
params.adFlipAngleDegree = 5;       % scalar or vector of flip angles in degrees
FitParams.numFA = numel(params.adFlipAngleDegree);

%% Preparation schedule for T1 recovery
params.moduleLength = 1;
params.numSRns = 0;
params.numIRns = 1;
params.numT2prep = 0;
params.numT1rhoPrep = 0;
params.numDiffPrep = 0;
params.numDiffPrepDirs = 0;
params.isPrep2Seg = 0;

%% Optional preparation flags, kept for compatibility with the Bloch model
params.isT2IR = 1;
params.isT1rhoIR = 1;
params.isDiffPrepIR = 1;

%% Minimal FitParams fields used by genT1Dictionary_VTR
FitParams.params = params;
FitParams.Nz = 1;                   % Use >1 for 3D; affects SMS/B1 scaling branch
FitParams.MBfactor = 1;
FitParams.ScanType = 'IR';          % 'IR' or 'SR'; controls default BIR grid
FitParams.flagUseBT2 = true;

%% Dictionary controls
dictOptions = struct();
dictOptions.cutoff = 10;            % Drop early recovery points before fitting
dictOptions.minT1 = 100e-3;         % seconds
dictOptions.maxT1 = 2.0;            % seconds
dictOptions.T1s = logspace(log10(dictOptions.minT1), log10(dictOptions.maxT1), 101);

dictOptions.Betas = linspace(0.2, 1.0, 5);
dictOptions.BIRs = linspace(0.6, 1.0, 3);
dictOptions.BT2s = linspace(0.6, 1.0, 3);

% Select which B1/inversion-efficiency dictionary slice is used for direct
% T1 map fitting. The full grid is still retained in dict.curves.
dictOptions.dictBeta = 1.0;
dictOptions.dictBIR = 1.0;
dictOptions.dictBT2 = 1.0;

% Optional override if params.alTE_seconds is unavailable:
% dictOptions.TRnav = 3.5e-3;

end
