% MR Multitasking reconstruction code v0.1 - VB/VE/XA multitasking sequence support
% This is the unified branch, composed by Hsu-Lei Lee
% It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories
% T1/T2/T1rho/VFA contrast
% This scipt is mainly for Human Study

%% Set up paths
mainpath='/home/leeh5/recon/multitaskingReconTool';                       % Put your path here.
addpath(genpath(fullfile(mainpath,'supporting','recon')));      % Supporting reconstruction scripts
addpath(genpath(fullfile(mainpath,'supporting','utils')));      % Supporting utilities
addpath(genpath(fullfile(mainpath,'supporting','colormaps')));  % Colormaps

% -- Third-party tool paths are below:

% -- Siemens rawdata import code. Included. 
addpath(fullfile(mainpath,'supporting','mapVBVD'));


% ======================================================================= 
%  finufft (fast NUFFT tool for radial data reconstruction)
% ======================================================================= 
% --  Linux --
% If you are NOT using recon3:
% Before first use, remove the mainpath/supporting/nufft/finufft folder
% Open a terminal in mainpath/supporting/nufft/ and run
% > git clone https://github.com/flatironinstitute/finufft.git
% > cd finufft && mkdir build && cd build
% > cmake -D FINUFFT_USE_CUDA=ON -D FINUFFT_BUILD_MATLAB=ON ..
% > cmake --build . -j
% Make sure the 'matlab' folder is under mainpath/supporting/nufft/finufft
% and the mex files are inside
%
% If need to build for a specific compute capability, find GPU info with
% > nvidia-smi --query-gpu=compute_cap --format=csv,noheader
% and add CMAKE_CUDA_ARCHITECTURES flag to the above cmake command
% > cmake -D FINUFFT_USE_CUDA=ON -D FINUFFT_BUILD_MATLAB=ON -D CMAKE_CUDA_ARCHITECTURES=75 ..
%
% --  Windows 10 (CPU-based)  --
% Before first use, download the finufft package from 
% https://github.com/flatironinstitute/finufft
% and copy the 'matlab' folder to mainpath\supporting\nufft\finufft\
%
% Download the below mex file to mainpath\supporting\nufft\finufft\matlab\
% https://users.flatironinstitute.org/~ahb/codes/finufft-binaries/2.0.2/win/finufft.mexw64
%
% If using Windows 10, copy all dll files under mainpath\supporting\nufft\finufft\winlib
% to C:\Windows\system32\
%
% ======================================================================= 
%  Matlab cufinufft interface (Linux only)
% ======================================================================= 
% Written by Junzhou Chen 
% Junzhou.Chen@med.usc.edu
%
% If CUDA path is not already in your bash LD_LIBRARY_PATH,
% run the following lines:
% (change CUDA_LIB to the proper CUDA library path on yur machine)
% > CUDA_LIB='/usr/local/cuda/lib64'
% > echo "export LD_LIBRARY_PATH=$CUDA_LIB:\$LD_LIBRARY_PATH" >> ~/.bashrc
% > source ~/.bashrc
%
% If you are NOT using recon3:
% Before first use, open mainpath/supporting/nufft/MT_CUDA/compile_mex.m,
% change the cuda settings according to your system environment, 
% then run the following lines in Matlab command window:
% > compile_mex
% > copyfile *.mexa64 ../finufft/matlab
%
% If cufinufft keeps failing, try open Matlab by running
% mainpath/matlab_cufinufft.sh
% (remember to check the paths in the script)
% =======================================================================
if exist(fullfile(mainpath,'supporting','nufft','finufft','matlab'),'dir')
    addpath(fullfile(mainpath,'supporting','nufft','finufft','matlab')); 
    fprintf('-- finufft path added.\n');
end


% ======================================================================= 
%  MIRT (optional, only use this if finufft is not available)
% ======================================================================= 
% Before first use, after the following addpath command, run
% > ir_mex_build
% to build the mex files
% =======================================================================
if exist(fullfile(mainpath,'supporting','nufft','irt'), 'dir')
    currentdir = pwd;
    irtdir = fullfile(mainpath,'supporting','nufft','irt');
    cd(irtdir); setup; cd(currentdir);
    clear currentdir irtdir;
end


% ======================================================================= 
%  gpuNUFFT (GPU-based NUFFT tool for radial data reconstruction)
% ======================================================================= 
% --  Linux Only  --
% Before first use, open a terminal in mainpath/supporting/nufft/ and run
% > git clone https://github.com/andyschwarzl/gpuNUFFT.git
% > cd gpuNUFFT/CUDA && mkdir -p build && cd build
% > cmake ..
% > make
% =======================================================================
if gpuDeviceCount && exist(fullfile(mainpath,'supporting','nufft','gpuNUFFT','gpuNUFFT'),'dir')
    addpath(fullfile(mainpath,'supporting','nufft','gpuNUFFT','gpuNUFFT'));  
    fprintf('-- gpuNUFFT path added.\n');
end


% ======================================================================= 
%  BART (optional)
% ======================================================================= 
% --  Linux Only  --
% Before first use, open a terminal in mainpath/supporting/nufft/
% and run
% > sudo apt-get install make gcc libfftw3-dev liblapacke-dev libpng-dev libopenblas-dev
% > git clone https://github.com/mrirecon/bart.git
% > cd bart
% > make
% =======================================================================
if exist(fullfile(mainpath,'supporting','nufft','bart','matlab'),'dir')
    addpath(fullfile(mainpath,'supporting','nufft','bart','matlab'));   
    setenv('TOOLBOX_PATH', fullfile(mainpath,'supporting','nufft','bart'));
    fprintf('-- BART path added.\n');
end


%% Initialize variable struct

Params = [];
ReconOptions = [];
TwixObj = [];
DataArray = [];
TemporalBasis = [];
SpatialCoeff = [];
ReconstructedImages = [];
FitParams = [];


%% Load data and preprocess

% -- Add break point when error occurs
% -- remove the following line if not debugging
dbstop if error

% -- store mainpath in ReconOptions
ReconOptions.mainpath = mainpath;

% -- Turn on image display when run in command line mode
ReconOptions.flagCommandLine = true;

% -- Set NUFFT package to use for radial data recon
% if
% flagUseBart = 1                      -- BART
% else
% flagUsefinufft = 1 && flagUseGPU = 1 -- cufinufft
% flagUsefinufft = 1 && flagUseGPU = 0 -- finufft
% flagUsefinufft = 0 && flagUseGPU = 1 -- gpuNUFFT
% flagUsefinufft = 0 && flagUseGPU = 0 -- irt

% -- Use BART?
ReconOptions.flagUseBart = 0; 

% -- Use finufft/cufinufft?
ReconOptions.flagUsefinufft = 1;

% -- Set GPU mode
ReconOptions.flagUseGPU = true;  % (exist('gpuNUFFT','file')>1); % true if using GPU to do nufft. false if using finufft/MIRT on CPU.
if ~ReconOptions.flagUsefinufft
    if ReconOptions.flagUseGPU && gpuDeviceCount && ((exist('gpuNUFFT','file')>1) || (exist('bart','file')>1))
        chooseGPU; 		% Give MATLAB main software access to the GPU before the gpuNUFFT software locks it up
    else
        ReconOptions.flagUseGPU = false;
    end
end

% -- Adjust NUFFT memory usage
% 0: over-sampling factor = 2; 1: 1.25 
ReconOptions.nufftlowmem = 0;     

% -- Select data file
Params = selectData;
[Params, ReconOptions, TwixObj, DataArray] = loadHeader(Params, ReconOptions);


%% Check multitasking parameters and reassign the values if necessary
% -- Should not need this if data was acquired with the unified sequence

% Params.ScanType       = 'T2IR';     % options: {'IR','IR_VFA','T2IR','T2IR_VFA','SR','Cine'}
% Params.numIRns        = 0;
% Params.numSRns        = 0;
% Params.numIRsel       = 0;
% Params.numT2prep      = 4;
% Params.T2prepDuration = [30 55 80 0];
% Params.numT1rhoPrep   = 0;
% Params.T1rhoDuration  = [];
% Params.isT2IR         = true;
% Params.numFA          = 1;
% for n = 1:Params.numFA
% 	Params.flipAngleArray(n) = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{n}];
% end
% Params.moduleLength   = lcm((Params.numIRns+Params.numSRns+Params.numT2prep+Params.numT1rhoPrep),Params.numFA);
% Params.SGBlock        = 4;
% Params.linesPerShot   = 144;
% Params.MBfactor       = 1;
% Params.navShift       = 1;
% Params.isTinyGoldenAngle = true;
% Params.isVTR          = false;

%% Pre-processing
% -- Load data into variable struct
ReconOptions.totalTime = inf;   % length of data to use, in seconds. "inf" to use all.
[Params, ReconOptions, DataArray] = loadData(Params, ReconOptions, TwixObj, DataArray);

% -- Setup trajectories
% Optional correction for radial data
ReconOptions.flagGradientDelayCorr = 1;     % 0: no gradient delay correction; 1: old method;   2: new method (slower)
ReconOptions.flagEddyCurrentCorr   = 1;     % 0: no Eddy current correction;   1: do correction
[Params, ReconOptions, DataArray] = setupTrajectories(Params, ReconOptions, DataArray);

% -- Pre-whiten
% Optional FOV crop for Cartesian data
ReconOptions.cropFOV = 1;      % FOVnew = FOV/cropFOV in readout direction
[Params, ReconOptions, DataArray] = preWhiten(Params, ReconOptions, DataArray, ReconOptions.cropFOV);

% -- Back-projection
[Params, ReconOptions, DataArray] = calcfbp(Params, ReconOptions, DataArray);

%%
% -- Coil compression params
ReconOptions.minNewCoils = 4 ;          % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.95;      % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = 16; %min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


%% Initial reconstruction

% -- Which slice to display?
ReconOptions.dispSlice = floor(Params.Nz/2) + 1;

% -- Estimate coil sensitivities
ReconOptions.SEmethod = 'CBD';  % options: {'CBD', 'Walsh', 'ESPIRiT'}
[Params, ReconOptions, DataArray] = estimateSensitivities(Params, ReconOptions, DataArray);

% -- Real-time recon
ReconOptions.L_rt = 4;         % Rank for real-time recon
[Params,ReconOptions,DataArray,TemporalBasis] = realtimeSubspace(Params, ReconOptions, DataArray);

ReconOptions.flagUseInitialGuess = false; 
ReconOptions.flagUseToeplitz = false;
ReconOptions.L = ReconOptions.L_rt;
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis);

% -- Show real-time images
ReconstructedImages = displayRealtime(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,ReconOptions.dispSlice,[],20);

%% Generate relaxation subspace curves

ReconOptions.flagDataDriven = true;

FitParams.minT1 = 100e-3;
FitParams.maxT1 = 3000e-3;
FitParams.minT2 = 10e-3;
FitParams.maxT2 = 300e-3;
[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis,5);


%% Binning: Respiratory

disp('Start binning...')

% -- Set # of respiratory bins
ReconOptions.rbins = 1;    

% -- Breathing rate filter range (cycle/min)
ReconOptions.BRlow  = 3;
ReconOptions.BRhigh = 50;

% -- Do respiratory binning
DataArray = binningResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);

% -- Swap R-bins to make end-of-respiration bin #1 if necessary
% DataArray = flipRidx(DataArray,true);

% -- Motion correction (reference: Ridx=1)
% [Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Binning: Cardiac

% -- Set # of cardiac bins
ReconOptions.cbins = 1; 

% -- Heart rate filter range (cycle/min)
ReconOptions.HRlow  = 40;    
ReconOptions.HRhigh = 130;

% -- Do cardiac binning
if ReconOptions.totalTime > 30
    % -- Change the last input to true to manually draw heart ROI
    [DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);
else
    % -- Another cardiac binning method assuming minimum heart rate variation
    [DataArray,TemporalBasis] = binningCardSelfGate(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
end

% -- Shift end-of-diastole phase to the first bin
% DataArray.diastoleIdx = 1;     %str2double(inputdlg({'Select end-of-diastole phase:'},'User Input',[1 35],{'1'}));
% DataArray = shiftHidx(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

disp('done.')

%% Tensor subspace estimation

ReconOptions.lowmem = 1;                    % 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.L_tensor = 12;                  % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;
ReconOptions.tempL = 500;
ReconOptions.flagUseFirstEcho = false;

ReconOptions.flagUseInitialGuess = false;   
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseGPUToeplitz = true; 

ReconOptions.flagAutoLambdaTensor = true;  % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

% -- Overwrite lr and sms if necessary
ReconOptions.tensor.lr  = 1e-3;
% ReconOptions.tensor.sms = [3e-6 1e-6];

% -- Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);

% -- Least-squares tensor recon
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% -- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);


%% TV recon

ReconOptions.flagAutoLambdaTV = true;     % automatically determines lambda?

% -- Manually set lambda, only used if flagAutoLambdaWavelet is false
% -- Increase the value to suppress noise, decrease if image is blurred
ReconOptions.tv.lambda = 3e-10;
ReconOptions.tv.alpha  = 1;

% --Do wavelet recon
ReconOptions.tv.flagForce2D  = 0;
ReconOptions.tv.flagContinue = 0;          % 1: use SpatialCoeff.U as initial value; 0: use U_tensor
ReconOptions.tv.maxIter = 25;              % maximum number of iterations    
   
[ReconOptions,SpatialCoeff] = tvRecon3D_aniso(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


% -- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = 144;
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);

% -- Show black-blood images
% ReconstructedImages = displayTensorBlackBlood(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages);

% -- Plot image curves
% curve = plotSelectedVoxel(ReconstructedImages.reconRecovery,[],[],[],[],1);   


%% Collect parameters needed for parametric fitting

FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);


%% Save workspace

saveMultitaskingWorkspace;
% save(Params.fileString,'Params','ReconOptions','DataArray','TemporalBasis','SpatialCoeff','TwixObj','FitParams','-v7.3'); 
% disp('Saved.')

if strcmp(get(0,'Diary'),'on')  
    diary off
end


%% Parametric fitting
% 
% -- load necessary variables for fitting from saved .mat file
% load 'fitParams.mat' 'FitParams';

% -- Set the time points in the recovery curve to be used in fitting
FitParams.imstart = 3;
FitParams.imskip  = 2;

% -- 0: k-means; 1~2: more pixels; 3: No mask; 4: BET (brain extraction from FSL)
FitParams.flagLargeROI = 4;
FitParams.BETfrac = 0.3;       % BET: larger value -> smaller mask

% -- which respiratory/cardiac phases/slice to fit?
FitParams.rphase  = 1;
FitParams.cphases = 1; 
FitParams.fitSlice = 1:FitParams.Nzdisp;

% -- Flip angle adjustment
FitParams.initGlobalBeta = 0.9;

% -- initial values
FitParams.initT1     = 900e-3;
FitParams.initT2     = 90e-3;
FitParams.initBIR    = 0.8;
FitParams.initBT2    = 0.8;

% -- lower & upper bounds
FitParams.minT1 = 100e-3;
FitParams.maxT1 = 3000e-3;
FitParams.minT2 = 20e-3;
FitParams.maxT2 = 300e-3;
FitParams.minBIR = 0.5;
FitParams.maxBIR = 1;
FitParams.minBT2 = 0.5; 
FitParams.maxBT2 = 1;
FitParams.minBeta = 0.1;
FitParams.maxBeta = 1;

% BIR  - RF efficiency for inversion pulses
% BT2  - RF efficiency for T2 prep
% Beta - RF efficiency for excitation pulses

% -- Fitting option
% --- for T1T2 1FA measurements:
% --- 0 - fit T1/T2/BIR/BT2/Beta;  1 - constant BT2;  2: BT2 = BIR; 
% --- 3 - constant BT2 & Beta;     4 - BT2 = BIR, constant Beta
FitParams.ChooseConst = 2;

% -- Do fitting
FitResult = paramFit(FitParams);


% % -- Check fitted curves
% curve = plotFitResult(FitResult,FitParams.Nzdisp/2,1);   


%% QSM
% if you want to do QSM with your own MEDI package, the complex-valued multi-echo images
% and required parameters can be generated with
[reconME,reconME_iso,QSMparams] = genMEimagesfromFitParams(FitParams);

% Otherwise use the following lines to generate QSM, SWI, and tSWI maps
FitParams.methodUnwrap = 'Laplacian';
FitParams.methodBackgroundFieldRemoval = 'VSHARP';     % Choose from {'VSHARP';'PDF'}
FitParams.methodQSM = 'MEDI';
FitParams.flagDoSWI = true;
FitParams.filterSize = 65;      % phase low-pass filter size for SWI calculation
FitParams.power_weighting = 4;  % phase weighting power for SWI calculation
FitParams.threshold = 1;        % phase mask threshold for SWI calculation
FitParams.chi_1 = 0;        % lower threshold for tSWI calculation
FitParams.chi_2 = 0.3;      % upper threshold for tSWI calculation

FitResult.QSMResult = runQSM(FitParams);

%% Save fit results

fitResultFileString = [FitParams.filePath '/multitasking_FitResult.mat'];
save(fitResultFileString,'FitResult','-v7.3'); 
disp('FitResult Saved.')

dbclear all
