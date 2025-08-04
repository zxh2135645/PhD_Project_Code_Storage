% MR Multitasking reconstruction code v0.1 - VB and VE multitasking sequence support
% This is the unified branch, composed by Hsu-Lei Lee
% It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories
% T1/T2/VFA contrast
% This scipt is majorly for Human Study

%% Set up paths
mainpath = '/home/leeh5/recon/multitaskingReconTool';           % Put your path here.
addpath(genpath(fullfile(mainpath,'supporting','recon')));      % Supporting reconstruction scripts
addpath(genpath(fullfile(mainpath,'supporting','utils')));      % Supporting utilities
addpath(genpath(fullfile(mainpath,'supporting','colormaps')));  % Colormaps

% -- Third-party tool paths are below:

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
    currentdir = pwd; cd(mainpath); mainpath = pwd;
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
    addpath(genpath(fullfile(mainpath,'supporting','nufft','bart','matlab')));   
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


%% Recon package license ID for Showa University Hospital

ReconOptions.licenseID = 'UTSW_MT';

%% Select data and load header

% -- Add break point when error occurs
% -- comment out the next line if not debugging
% dbstop if error

% -- store mainpath in ReconOptions
ReconOptions.mainpath = mainpath;

% -- Turn on image display when run in command line mode
ReconOptions.flagCommandLine = true;

% -- Set GPU mode
ReconOptions.flagUseGPU = false;  % (exist('gpuNUFFT','file')>1); % true if using GPU to do nufft. false if using finufft/MIRT on CPU.
if ReconOptions.flagUseGPU && gpuDeviceCount && ((exist('gpuNUFFT','file')>1) || (exist('bart','file')>1))
    chooseGPU; 		% Give MATLAB main software access to the GPU before the gpuNUFFT software locks it up
else
    ReconOptions.flagUseGPU = false;
end

% -- Set NUFFT package to use for radial data recon in CPU mode
% 1: finufft(default, much faster); 0: irt
ReconOptions.flagUsefinufft = 1;

% -- Adjust NUFFT memory usage
% 0: over-sampling factor = 2.0; otherwise: 1.25 
ReconOptions.nufftlowmem = 0;                    

% -- Use BART?
ReconOptions.flagUseBart = 0;                   

% -- Select data file
Params = selectData;
[Params, ReconOptions, TwixObj, DataArray] = loadHeader(Params, ReconOptions);


%% Check multitasking parameters and reassign the values if necessary

% For spiralVIBE scans
% Params.ScanType       = 'Cine';
% Params.numSRns        = 0;
% Params.numIRns        = 0;
% Params.numT2prep      = 0;
% Params.T2prepDuration = [];
% Params.numT1rhoPrep   = 0;
% Params.T1rhoDuration  = [];
% Params.numFA          = 1;
% for n = 1:Params.numFA
%     Params.flipAngleArray(n) = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{n}];
% end
% Params.moduleLength   = 1;
% Params.SGBlock        = 1;
% Params.linesPerShot   = Params.lSegments;
% Params.MBfactor       = 1;
% Params.navShift       = 1;
% Params.Ntrajs         = 0;
% Params.trajInc        = 0;
% Params.isTinyGoldenAngle = false;
% Params.isVTR          = false;

%% Pre-processing

% Load data into variable struct
ReconOptions.totalTime = inf; % amount of data to use, in seconds. Use "inf" to use all.
[Params, ReconOptions, DataArray] = loadData(Params, ReconOptions, TwixObj, DataArray);

% -- Setup trajectories
% -- Crop recon FOV in radial casecd
% -- in abdomen cases use 2 for A-P dimension to reduce recon loading
ReconOptions.cropFOVy = 1;              % FOVynew = FOVy/cropFOVy for radial trajectory
ReconOptions.cropFOVx = 1;              % FOVxnew = FOVx/cropFOVy for radial trajectory
% -- Optional correction for radial data
ReconOptions.flagGradientDelayCorr = 0;     % 0: no gradient delay correction; 1: global grad delay; 2: angle-dependent (slower)
ReconOptions.flagEddyCurrentCorr   = 0;     % 0: no eddy current correction;   1: do correction
[Params, ReconOptions, DataArray] = setupTrajectories(Params, ReconOptions, DataArray);

% -- Pre-whiten
ReconOptions.cropRO = 1;               % FOVnew = FOV/cropFOV
[Params, ReconOptions, DataArray] = preWhiten(Params, ReconOptions, DataArray, ReconOptions.cropRO);

% -- Back-projection
[Params, ReconOptions, DataArray] = calcfbp(Params, ReconOptions, DataArray);

%%
% -- Coil compression params
ReconOptions.minNewCoils = 4;          % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.99;      % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = 4; %min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


%% Initial reconstruction

ReconOptions.dispSlice = 80; %floor(Params.Nz/2) + 1;    % Use central slice for display

% Estimate coil sensitivities
ReconOptions.SEmethod = 'CBD';  % {'CBD', 'Walsh', 'ESPIRiT'}
[Params, ReconOptions, DataArray] = estimateSensitivities(Params, ReconOptions, DataArray);

% Real-time recon
ReconOptions.L_rt = 16;         % Rank for real-time recon
[Params,ReconOptions,DataArray,TemporalBasis] = realtimeSubspace(Params, ReconOptions, DataArray);

ReconOptions.flagUseInitialGuess = false; 
ReconOptions.flagUseToeplitz = false;
ReconOptions.L = ReconOptions.L_rt;
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis);

% Show real-time images
ReconstructedImages = displayRealtime(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,ReconOptions.dispSlice,[],20);

%%

% -- Generate relaxation subspace curves
FitParams.minT1 = 300e-3;
FitParams.maxT1 = 3000e-3;
[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis);


%% binning
disp('Starting binning...')

% Respiratory binning
% ReconOptions.rbins = 1;
% [DataArray,ReconOptions] = binningResp_DCE(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% -- Set # of respiratory bins
ReconOptions.rbins = 6;    

% -- Breathing rate filter range (cycle/min)
ReconOptions.BRlow  = 3;
ReconOptions.BRhigh = 30;

% -- Do respiratory binning
DataArray = binningResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);

% Swap R-bins if necessary
% DataArray = flipRidx(Params,DataArray);

%saveGif(1.2*abs(DataArray.binsRespMean),'.','binRespMean.gif',2);

% Motion correction (reference: Ridx=1)
%[Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


% Cardiac binning
ReconOptions.cbins = 1;
[DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);


%% Tensor subspace estimation

ReconOptions.lowmem = 1;                    % 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.L_tensor = 8;                  % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;

ReconOptions.NDCEAve = 2;                      % DCE temporal resolution = NDCEAve * TR

ReconOptions.flagAutoLambdaTensor = true;   % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

% Overwrite lr and sms if necessary
% ReconOptions.tensor.lr  = 1e-4;
% ReconOptions.tensor.sms = [5e-6 1e-6];
% ReconOptions.tensor.sms = ReconOptions.tensor.sms * 1000;

% Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);

% Least-squares tensor recon
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseInitialGuess = false;  
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% Show tensor images
DispParams.scale = 1;
DispParams.tIdx = 150;
DispParams.cIdx = 1;
DispParams.rIdx = 1;
DispParams.DCEIdx = 1;
DispParams.eIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice,DispParams.DCEIdx,DispParams.eIdx);

% Save workspace
saveMultitaskingWorkspace;


%% Wavelet recon

ReconOptions.flagAutoLambdaWavelet = false;     % automatically determines lambda?
ReconOptions.flagUseBart = false;

% -- Manually set lambda, only used if flagAutoLambdaWavelet is false
% -- Increase the value to suppress noise, decrease if image is blurred
% ReconOptions.wavelet.lambda = 2e-10;
% ReconOptions.wavelet.alpha  = 1;

% --Do wavelet recon
ReconOptions.wavelet.flagContinue = 0;     % 1: use current U as initial value; 0: use U from tensor recon
ReconOptions.wavelet.maxIter = 5;          % maximum number of iterations         
[ReconOptions,SpatialCoeff]  = waveletRecon(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Show tensor images
DispParams.scale = 1.2;
DispParams.tIdx = 120;
DispParams.cIdx = 1;
DispParams.rIdx = 1;
DispParams.DCEIdx = 100;
DispParams.eIdx = 1;
ReconOptions.dispSlice =16;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice,DispParams.DCEIdx,DispParams.eIdx);

%% Collect parameters needed for parametric fitting
FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);

%% Save workspace

saveMultitaskingWorkspace;

if strcmp(get(0,'Diary'),'on')  
    diary off
end

%% Parametric fitting

%load 'multitasking.mat' 'FitParams';

% Collect parameters needed for parametric fitting
%FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);


% Set the time points in the recovery curve to be used in fitting
FitParams.imstart = 1;
FitParams.imskip  = 1;

% 0: k-means; 1: force ROI in central FOV/4; 2: force ROI in central FOV/2 
FitParams.flagLargeROI = 0;

% which respiratory/cardiac phases to fit?
FitParams.rphase = 1;
FitParams.cphase = 1;
FitParams.fitSlice = 1:Params.Nz;

% Fitting option
% for T1T2 1FA measurements:
% 0 - fit T1/BIR/beta;  1 - constant beta
FitParams.ChooseConst = 1;

% initial values
initT1     = 1;         % sec
initBIR    = 1;         % SR

% lower & upper bounds
FitParams.minT1 = 50e-3;
FitParams.maxT1 = 5;
FitParams.minBIR = 0.5;
FitParams.maxBIR = 1.5;
FitParams.minBeta = 0.6;
FitParams.maxBeta = 1;

% Flip angle adjustment
FitParams.initGlobalBeta = 0.6;

% Assign the ROI mask if you have one
%FitParams.mask = [];

% DCE parameters
FitParams.preCutoff   = 50;
FitParams.contrastInj = 35;

% Do fitting
% DCE: fit T1
%FitResult1 = DCE_fit_T1(FitParams);

% ME: fit T1 and estimate fat fraction
FitResult1 = paramFit(FitParams);
FitResult2 = fit_ME_FF(FitParams);

% Plot all voxels
%figure;plot(FitResult.T1');

% Display maps as movie
% nphases = numel(FitResult1.T1map)/numel(FitResult1.mask);
% T1map = zeros(numel(FitResult1.mask), nphases);
% T1map(logical(FitResult1.mask(:)),:) = reshape(FitResult1.T1map,[],nphases);
% T1map = reshape(T1map,[size(FitResult1.mask) nphases]);

T1map = FitResult1.T1map;
dispSlice = floor(size(T1map,3)/2) + 1;
T1map_color = displayMap(reshape(FitResult2.T2star,size(FitResult2.T2star,1),size(FitResult2.T2star,2),1,[]),inferno(256),[0 0.15],5,'temp.gif',FitParams.filePath);
 
% % Display voxel value change
% curve = plotSelectedVoxel(FitResult.T1map(:,:,:)*1000,[0 3000],'inferno',1);    % T1 display range = [0 3000]; slice = 1;

%% Save fit results

fitResultFileString = [FitParams.filePath '/multitasking_FitResult_MID00084_FID06974_biascorrection.mat'];
save(fitResultFileString,'FitResult1','FitResult2','-v7.3'); 
disp('FitResult Saved.');

%% Export T1map to dicom
load 'multitasking.mat' 'TwixObj';
load 'multitasking_FitResult';

genDicomfromSiemens(abs(round(1000*FitResult.T1map)), TwixObj, 0, 'T1map');

