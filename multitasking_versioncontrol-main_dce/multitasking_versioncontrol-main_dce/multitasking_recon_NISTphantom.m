% MR Multitasking reconstruction code v0.1 - VB and VE multitasking sequence support
% This is the unified branch, composed by Hsu-Lei Lee
% It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories
% T1/T2/VFA contrast
% This scipt is majorly for Human Study

%% Set up paths
mainpath = '/home/leeh5/recon/multitaskingReconTool';                       % Put your path here.
addpath(genpath(fullfile(mainpath,'supporting','recon')));      % Supporting reconstruction scripts
addpath(genpath(fullfile(mainpath,'supporting','utils')));      % Supporting utilities
addpath(genpath(fullfile(mainpath,'supporting','colormaps')));  % Colormaps

% -- Third-party paths are below:

% -- Siemens rawdata import code. Included. 
addpath(fullfile(mainpath,'supporting','mapVBVD'));


% ======================================================================= 
%  finufft (fast CPU-based NUFFT tool for radial data reconstruction)
% ======================================================================= 
% --  Linux  --
% Before first use, open a terminal and run
% > git clone https://github.com/flatironinstitute/finufft.git
% > cd finufft
% > make clean && make test -j && make matlab
% Then copy the 'matlab' folder to mainpath/supporting/nufft/finufft
%
% --  Windows 10  --
% Before first use, download the finufft package from 
% https://github.com/flatironinstitute/finufft
% and copy the 'matlab' folder to mainpath\supporting\nufft\finufft\
%
% Download the below mex file to mainpath\supporting\nufft\finufft\matlab\
% https://users.flatironinstitute.org/~ahb/codes/finufft-binaries/2.0.2/win/finufft.mexw64
%
% If using Windows 10, copy all dll files under mainpath\supporting\nufft\finufft\winlib
% to C:\Windows\system32\
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
% dbstop if error

% -- store mainpath in ReconOptions
ReconOptions.mainpath = mainpath;

% -- Turn on image display when run in command line mode
ReconOptions.flagCommandLine = true;

% -- Set GPU mode
ReconOptions.flagUseGPU  = false;%(exist('gpuNUFFT','file')>1); % true if using gpuNUFFT. false if using finufft/MIRT.
if ReconOptions.flagUseGPU && gpuDeviceCount
    chooseGPU; 		% Give MATLAB main software access to the GPU before the gpuNUFFT software locks it up
else
    ReconOptions.flagUseGPU = false;
end

% -- Set NUFFT package to use for radial data recon in CPU mode
% 1: finufft(default, much faster); 0: irt
ReconOptions.flagUsefinufft = 1;

% -- Adjust NUFFT memory usage
% 0: over-sampling factor = 2; 1: 1.25 
ReconOptions.nufftlowmem = 0;                    

% -- Use BART?
ReconOptions.flagUseBart = 1;

% -- Select data file
Params = selectData;
[Params, ReconOptions, TwixObj, DataArray] = loadHeader(Params, ReconOptions);


%% Check multitasking parameters and reassign the values if necessary

Params.ScanType       = 'T2IR';
Params.numIRns        = 0;
Params.numT2prep      = 4;
Params.T2prepDuration = [30 40 60 80];
Params.numT1rhoPrep   = 0;
Params.T1rhoDuration  = [];
Params.numFA          = 1;
Params.flipAngleArray = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{1}];
Params.moduleLength   = 4;
Params.SGBlock        = 16;
Params.linesPerShot   = Params.lSegments;
Params.MBfactor       = 1;
Params.navShift       = 1;
Params.isTinyGoldenAngle = false;
Params.isVTR          = false;

%% Pre-processing
% -- Load data into variable struct
ReconOptions.totalTime = inf;   % amount of data to use, in seconds. Use "inf" to use all.
[Params, ReconOptions, DataArray] = loadData(Params, ReconOptions, TwixObj, DataArray);

% -- Setup trajectories
% Optional correction for radial data
ReconOptions.flagGradientDelayCorr = 1;     % 0: no gradient delay correction; 1: old method;   2: new method (slower)
ReconOptions.flagEddyCurrentCorr   = 1;     % 0: no eddy current correction;   1: do correction
[Params, ReconOptions, DataArray] = setupTrajectories(Params, ReconOptions, DataArray);

% -- Pre-whiten
ReconOptions.cropFOV = 1;               % readFOVnew = readFOV/cropFOV
[Params, ReconOptions, DataArray] = preWhiten(Params, ReconOptions, DataArray, ReconOptions.cropFOV);

% -- Back-projection
[Params, ReconOptions, DataArray] = calcfbp(Params, ReconOptions, DataArray);

% -- Coil compression params
ReconOptions.minNewCoils = 12;          % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.99;      % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = 12; %min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


%% Initial reconstruction

ReconOptions.dispSlice = 49;

% -- Estimate coil sensitivities
ReconOptions.SEmethod = 'ESPIRiT';  % options: {'CBD', 'Walsh', 'ESPIRiT'}
[Params, ReconOptions, DataArray] = estimateSensitivities(Params, ReconOptions, DataArray);

% -- Real-time recon
ReconOptions.L_rt = 5;         % Rank for real-time recon
[Params,ReconOptions,DataArray,TemporalBasis] = realtimeSubspace(Params, ReconOptions, DataArray);

ReconOptions.flagUseInitialGuess = false; 
ReconOptions.flagUseToeplitz = false;
ReconOptions.L = ReconOptions.L_rt;
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis);

% -- Show real-time images
ReconstructedImages = displayRealtime(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,ReconOptions.dispSlice,[],20);
%saveGif(DataArray.reconRealtime(:,:,1,:),Params.filePath,'_reconRealtime.gif',1/DataArray.reconRealtimeFrameTime);

% -- Generate relaxation subspace curves
FitParams.minT1 = 500e-3;
FitParams.maxT1 = 2000e-3;
FitParams.minT2 = 10e-3;
FitParams.maxT2 = 300e-3;
FitParams.minT1rho = 30e-3;
FitParams.maxT1rho = 100e-3;
[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis);



%% binning
disp('Starting binning...')
ReconOptions.rbins = 1;  %Set # of respiratory bins here

ReconOptions.cbins = 1; %Set #do_binning; of cardiac bins here

% ReconOptions.flagUseEMD = false;    % For testing only, do not use

% Respiratory binning
DataArray = binningResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);

% Swap R-bins if necessary
% DataArray = flipRidx(DataArray);

% Motion correction (reference: Ridx=1)
%[Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Cardiac binning
if ReconOptions.totalTime > 30
    [DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);

    % Test different parameters and choose the best result. Very slow (and does not work well...)
    %[DataArray,TemporalBasis] = binningCardSweep(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
else
    % Another cardiac binning method assuming minimum heart rate variation
    [DataArray,TemporalBasis] = binningCardSelfGate(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
end

% Shift end-of-diastole phase to the first bin
% DataArray.diastoleIdx = 1;     %str2double(inputdlg({'Select end-of-diastole phase:'},'User Input',[1 35],{'1'}));
% DataArray = shiftHidx(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


%% Tensor subspace estimation

ReconOptions.lowmem = 0;                    % 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.L_tensor = 8;                 % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;

ReconOptions.flagAutoLambdaTensor = false;  % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

% Overwrite lr and sms if necessary
% ReconOptions.tensor.lr  = 7e-10;
% ReconOptions.tensor.sms = [3e-6 1e-6];
ReconOptions.tensor.sms = [0 0];

% Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);

% Least-squares tensor recon 
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseInitialGuess = false;  
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);


%% Wavelet recon

% -- Manually set lambda, only used if flagAutoLambdaWavelet is false
% -- Increase the value to suppress noise, decrease if image is blurred
%ReconOptions.wavelet.lambda  = 2.5808e-15;
ReconOptions.wavelet.lambda = 2e-12;
ReconOptions.wavelet.alpha  = 1;

ReconOptions.flagAutoLambdaWavelet = false;     % automatically determines lambda?
ReconOptions.wavelet.flagContinue = 1;      % 1: use current U as initial value; 0: use U from tensor recon
ReconOptions.wavelet.maxIter = 10;          % maximum number of iterations         
[ReconOptions,SpatialCoeff]  = waveletRecon(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);


%% Save cine images as dicom
%genDicomfromSiemens(abs(DataArray.reconCard), TwixObj, 1, 'Cine');

%% Collect parameters needed for parametric fitting

FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);

%% Save workspace

saveMultitaskingWorkspace;
disp('Saved.')
diary off

%% Parametric fitting

%load 'multitasking.mat' 'FitParams';

% Set the time points in the recovery curve to be used in fitting
FitParams.imstart = 3;
FitParams.imskip  = 10;

% 0: k-means; 1: force ROI in central FOV/4; 2: force ROI in central FOV/2 
FitParams.flagLargeROI = 0;

% which respiratory/cardiac phases to fit?
FitParams.rphase  = 1;
FitParams.cphases = 1:ReconOptions.cbins;
FitParams.fitSlice = 49; %1:Params.Nz;

% -- Flip angle adjustment
FitParams.initGlobalBeta = 1;

% initial values
FitParams.initT1     = 1;
FitParams.initT2     = 60e-3;
FitParams.initBIR    = 1;
FitParams.initBT2    = 1;

% lower & upper bounds
FitParams.minT1 = 100e-3;
FitParams.maxT1 = 5;
FitParams.minT2 = 20e-3;
FitParams.maxT2 = 300e-3;
FitParams.minBIR = 0.5;
FitParams.maxBIR = 1;
FitParams.minBeta = 1;   %FitParams.initGlobalBeta;
FitParams.maxBeta = 1;   %FitParams.initGlobalBeta;
FitParams.minBT2 = 0.5;
FitParams.maxBT2 = 1;

% Fitting option
% for T1T2 1FA measurements:
% 0 - fit T1/T2/BIR/BT2/Beta;  1 - constant BT2;  2: BT2 = BIR; 
% 3 - constant BT2 & Beta;     4 - BT2 = BIR, constant Beta
FitParams.ChooseConst = 4;

% Do fitting
FitResult = paramFit(FitParams);

% Display maps as movie
%T1map_color = displayMap(FitResult.T1map,inferno(256),[0 3],20,'T1map_rphase1.gif',FitParams.filePath);
%T2map_color = displayMap(FitResult.T2map,inferno(256),[0 300],20,'T1map_rphase1.gif',FitParams.filePath);

% Display voxel value change
curve = plotSelectedVoxel(ReconstructedImages.reconRecovery,[],[],[],FitResult,1);    % T1 display range = [0 3000]; slice = 1;

%% Save fit results

fitResultFileString = [FitParams.filePath '/multitasking_FitResult.mat'];
save(fitResultFileString,'FitResult','-v7.3'); 
disp('FitResult Saved.')

%% Export T1map to dicom
% load 'multitasking.mat' 'TwixObj';
% load 'multitasking_FitResult';

genDicomfromSiemens(abs(round(wmedfilt2(1000*FitResult.T1map))), TwixObj, [1 1], 0, 'T1map');

