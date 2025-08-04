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

ReconOptions.licenseID = 'UTSW_MT_DCE';

%% Select data and load header

% -- Add break point when error occurs
% -- comment out the next line if not debugging
% dbstop if error

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
% 0: over-sampling factor = 2.0; otherwise: 1.25 
ReconOptions.nufftlowmem = 1;                    


% -- Select data file
Params = selectData;
[Params, ReconOptions, TwixObj, DataArray] = loadHeader(Params, ReconOptions);


%% Check multitasking parameters and reassign the values if necessary

% check whether parameters correct
% Params.ScanType       = 'SR';
% Params.numSRns        = 1;
% Params.numIRns        = 0;
% Params.numT2prep      = 0;
% Params.T2prepDuration = [];
% Params.numT1rhoPrep   = 0;
% Params.T1rhoDuration  = [];
% Params.numFA          = 1;
% Params.flipAngleArray = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{1}];
% Params.moduleLength   = 1;
% Params.SGBlock        = 8;
% Params.linesPerShot   = Params.lSegments;
% Params.MBfactor       = 1;
% Params.navShift       = 1;
% Params.Ntrajs         = 0;
% Params.trajInc        = 0;
% Params.isTinyGoldenAngle = true;
% Params.isVTR          = false;

%% Pre-processing
% -- Load data into variable struct
ReconOptions.totalTime = inf;               % amount of data to use, in seconds. Use "inf" to use all.
[Params, ReconOptions, DataArray] = loadData(Params, ReconOptions, TwixObj, DataArray);

% -- Setup trajectories
% -- Optional correction for radial data
ReconOptions.flagGradientDelayCorr = 1;     % 0: no gradient delay correction; 1: old method;   2: new method (slower)
ReconOptions.flagEddyCurrentCorr   = 1;     % 0: no Eddy current correction;   1: do correction
[Params, ReconOptions, DataArray] = setupTrajectories(Params, ReconOptions, DataArray);

% -- Pre-whiten k-space data
ReconOptions.cropFOV = 1;                   % readFOVnew = readFOV/cropFOV
[Params, ReconOptions, DataArray] = preWhiten(Params, ReconOptions,DataArray, ReconOptions.cropFOV);

% -- Back-projection
[Params, ReconOptions, DataArray] = calcfbp(Params, ReconOptions, DataArray);

% -- Coil compression params
ReconOptions.minNewCoils = 12;              % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.99;          % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = 16; %min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


%% Initial reconstruction

% -- Display option
ReconOptions.dispSlice = floor(Params.Nz/2) + 1;

% -- Estimate coil sensitivities
ReconOptions.SEmethod = 'CBD';      % Options: {'CBD', 'Walsh', 'ESPIRiT'}
[Params, ReconOptions, DataArray] = estimateSensitivities(Params, ReconOptions, DataArray);

% -- Real-time recon
ReconOptions.L_rt = 16;             % Rank for real-time recon
[Params,ReconOptions,DataArray,TemporalBasis] = realtimeSubspace(Params, ReconOptions, DataArray, TemporalBasis);

ReconOptions.flagUseInitialGuess = false; 
ReconOptions.flagUseToeplitz = false;
ReconOptions.L = ReconOptions.L_rt;
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis);

% -- Show real-time images
ReconstructedImages = displayRealtime(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,ReconOptions.dispSlice);
%saveGif(DataArray.reconRealtime(:,:,1,:),Params.filePath,'reconRealtime.gif',1/DataArray.reconRealtimeFrameTime);

% -- Generate relaxation subspace curves
[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis);


%% binning
disp('Starting binning...')

% Respiratory binning
ReconOptions.rbins = 6;
[DataArray,ReconOptions] = binningResp_DCE(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Swap R-bins if necessary
% DataArray = flipRidx(DataArray);

% Motion correction (reference: Ridx=1)
[Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


% Cardiac binning
ReconOptions.cbins = 1;
[DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);



%% Tensor subspace estimation

ReconOptions.lowmem = 1;                    % 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.L_tensor = 16;                 % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;

ReconOptions.flagAutoLambdaTensor = true;  % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

% Overwrite lr and sms if necessary
% ReconOptions.tensor.lr  = 1e-4;
% ReconOptions.tensor.sms = [5e-6 1e-6];
%ReconOptions.tensor.sms = [2e-7 1.7e-7];

% Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);

% Least-squares tensor recon
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseInitialGuess = true;  
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% Show tensor images
DispParams.scale = 1;
DispParams.tIdx = 168;
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor_DCE(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);

% Save workspace
saveMultitaskingWorkspace;


%% Wavelet recon


ReconOptions.flagAutoLambdaWavelet = false;     % automatically determines lambda?
ReconOptions.flagUseBart = false;

% -- Manually set lambda, only used if flagAutoLambdaWavelet is false
% -- Increase the value to suppress noise, decrease if image is blurred
%ReconOptions.wavelet.lambda  = 2.5808e-15;
ReconOptions.wavelet.lambda = 2e-11;
ReconOptions.wavelet.alpha  = 1;

% --Do wavelet recon
ReconOptions.wavelet.flagContinue = 0;      % 1: use current U as initial value; 0: use U from tensor recon
ReconOptions.wavelet.maxIter = 10;          % maximum number of iterations         
[ReconOptions,SpatialCoeff]  = waveletRecon(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% Show tensor images
DispParams.scale = 1;
DispParams.tIdx = 160;
DispParams.cIdx = 1;
DispParams.rIdx = 1;
DispParams.DCEIdx = 11;
ReconOptions.dispSlice = 14;

ReconstructedImages = displayTensor_DCE(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice,DispParams.DCEIdx,3);


%% Save cine images as dicom
%genDicomfromSiemens(abs(DataArray.reconVolume), TwixObj, [1 1], 1, 'Phase0');


%% Collect parameters needed for parametric fitting

%FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);

%% Save workspace

saveMultitaskingWorkspace;

if strcmp(get(0,'Diary'),'on')  
    diary off
end
