% MR Multitasking reconstruction code v0.1 - VB/VE/XA multitasking sequence support
% This is the unified branch, composed by Hsu-Lei Lee
% It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories
% T1/T2/T1rho/VFA contrast
% This scipt is mainly for Human Study

%% Set up paths
mainpath='/home/leeh5/recon/multitaskingReconTool';             % Put your path here.
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

% -- Choose NUFFT package for radial data recon
% if
% flagUseBart = 1                      -- BART
% else
% flagUsefinufft = 1 && flagUseGPU = 1 -- cufinufft (fastest)
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

% Params.ScanType       = 'CEST';     % will be automatically determined in loadData()   
%                                     % options: {'IR','IR_VFA','T2IR','T2IR_VFA','SR','Cine','CEST'}
% Params.SGBlock        = 8;
% Params.linesPerShot   = 8;
% Params.numIRns        = 0;
% Params.numSRns        = 0;
% Params.numT2prep      = 0;
% Params.T2prepDuration = [];
% Params.numT1rhoPrep   = 0;
% Params.T1rhoDuration  = [];
% Params.numFA          = 1;
% for n = 1:Params.numFA
% 	Params.flipAngleArray(n) = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{n}];
% end
% if strcmp(Params.ScanType,'CEST')
%     Params.CESTMetabolite = 'GAG';
%     Params.CESTSatFA      = 500;
% %     Params.CESTSatFreqOffsetppmList = [ 300, 300, 300, -40.0, -30.0, -20.0, -15.0, -10.0, -9.0, -8.0, -7.0, ...
% %                                       -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.7, -1.4, -1.2, ...
% %                                       -1.0, -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, ...
% %                                        1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, ...
% %                                        7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 300, 300, 300];	% in ppm;
%     Params.CESTSatFreqOffsetppmList = [ 300.0, 300.0, 300.0, -100.0, -40.0, -30.0, -20.0, -15.0, -10.0, -9.0, ...
%                                         -8.0, -7.0, -6.5, -6.0, -5.5, -5, -4.5, -4.0, -3.5, -3, -2.5, -2.0, -1.5, ...
%                                         -1.0, -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, ...
%                                         1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, ...
%                                         7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 100.0, 300.0, 300.0, 300.0 ];
% 
%     Params.CESTnumFreqOffsets = numel(Params.CESTSatFreqOffsetppmList);
%     Params.CESTNumRep = numel(TwixObj.image.Lin)/Params.CESTnumFreqOffsets/Params.linesPerShot;
% else
%     try
%         Params.moduleLength = lcm((Params.numIRns+Params.numSRns+Params.numT2prep+Params.numT1rhoPrep),Params.numFA);
%     catch
%         Params.moduleLength = Params.numFA;
%     end
% end
% 
% Params.MBfactor       = 1;
% Params.navShift       = 1;
% Params.isTinyGoldenAngle = false;
% Params.isVTR             = false;


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

%% Coil Compression

% -- Coil compression params
ReconOptions.minNewCoils = 8;          % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.97;      % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


%% Initial reconstruction

% -- Which slice to display?
ReconOptions.dispSlice = floor(Params.Nz/2) + 1;

% -- Estimate coil sensitivities
ReconOptions.SEmethod = 'CBD';  % options: {'CBD', 'Walsh', 'ESPIRiT'}
[Params, ReconOptions, DataArray] = estimateSensitivities(Params, ReconOptions, DataArray);

% -- Real-time recon
ReconOptions.L_rt = 16;         % Rank for real-time recon
[Params,ReconOptions,DataArray,TemporalBasis] = realtimeSubspace(Params, ReconOptions, DataArray);

ReconOptions.flagUseInitialGuess = false; 
ReconOptions.flagUseToeplitz = false;
ReconOptions.L = ReconOptions.L_rt;
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis);

% -- Show real-time images
ReconstructedImages = displayRealtime(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,ReconOptions.dispSlice,[],20);
%saveGif(DataArray.reconRealtime(:,:,1,:),Params.filePath,'_reconRealtime.gif',1/DataArray.reconRealtimeFrameTime);

%% Signal contrast model

% -- Generate relaxation subspace curves
ReconOptions.flagDataDriven = false;

% Currently supporting {'APT','Creatine','Glucose','GAG'}
ReconOptions.CESTMetabolite = Params.CESTMetabolite;  
% N-pool model: currently supporting 2/3/4
% water - CEST - MT - NOE
ReconOptions.CESTnPools = 4;  

[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis, 5, DataArray);


%% Binning: Respiratory

disp('Start binning...')

% -- Set # of respiratory bins
ReconOptions.rbins = 6;    

% -- Breathing rate filter range (cycle/min)
ReconOptions.BRlow  = 3;
ReconOptions.BRhigh = 50;

% -- Use PMU signal for binning
ReconOptions.usePMUdata = true;

% -- Do respiratory binning
DataArray = binningResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% -- Swap R-bins to make end-of-respiration bin #1 if necessary
% DataArray = flipRidx(Params,DataArray,true);

% -- Motion correction (reference: Ridx=1)
% [Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


%% Binning: Cardiac

% -- Set # of cardiac bins
ReconOptions.cbins = 20; 

% -- Heart rate filter range (cycle/min)
ReconOptions.HRlow  = 45;    
ReconOptions.HRhigh = 100;

% -- Do cardiac binning
% if ReconOptions.totalTime > 30
%     % -- Change the last input to true to manually draw heart ROI
%     [DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,false);
% else
    % -- Another cardiac binning method assuming minimum heart rate variation
    [DataArray,TemporalBasis] = binningCardSelfGate(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
% end

% -- Shift end-of-diastole phase to the first bin
% DataArray.diastoleIdx = 11;     %str2double(inputdlg({'Select end-of-diastole phase:'},'User Input',[1 35],{'1'}));
% DataArray = shiftHidx(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


%% Tensor subspace estimation
% -- Low memory options
% 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.lowmem = 2;          % cardiac rbins = 6, cbins = 20
%ReconOptions.lowmem = 1;         % rbins = cbins = 1

ReconOptions.L_tensor = 32;                  % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;

ReconOptions.flagUseInitialGuess = false;   
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseGPUToeplitz = false; 

ReconOptions.flagAutoLambdaTensor = false;    % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

% -- Overwrite lr and sms if necessary
% ReconOptions.tensor.lr  = 1e-4;
ReconOptions.tensor.sms = [1.7e-5 3.6e-5];  % cardiac; lowmem = 2
% ReconOptions.tensor.sms = [7e-5];     % rbins = cbins = 1; flagDataDriven

% -- Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);

% -- Least-squares tensor recon
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% -- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[]);


%% Wavelet recon

ReconOptions.flagAutoLambdaWavelet = false;     % automatically determines lambda?
                                                % (doesn't work well)
% -- Overwrite lambda & alpha if necessary
ReconOptions.wavelet.lambda = 3e-15;
ReconOptions.wavelet.alpha  = 1;

% -- Do wavelet recon
ReconOptions.wavelet.flagContinue = 0;          % 1: use current U as initial value; 0: use U_tensor
ReconOptions.wavelet.maxIter      = 30;         % maximum number of iterations
[ReconOptions,SpatialCoeff] = waveletRecon(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% -- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[]);


% -- Plot voxel z spectrum
curve = plotSelectedVoxel(ReconstructedImages.reconCEST(:,:,1,:),[0 1]);   


%% Export CEST images as dicom

% -- interpolate the image if necessary
% [ReconstructedImages.reconCard,ReconOptions.interpFactor] = interp2dimage(abs(ReconstructedImages.reconCEST),[0.9375 0.9375],'resolution',TwixObj);

% -- Write dicom file
% genDicomfromSiemens(ReconstructedImages.reconCEST, TwixObj, ReconOptions.interpFactor, 1, 'CEST');


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
% load 'multitasking.mat' 'FitParams';
% load 'multitasking.mat' 'TwixObj';

% -- 0: k-means; 1~2: more pixels; 3: No mask 
FitParams.flagLargeROI = 0;

% -- which respiratory/cardiac phases/slice to fit?
FitParams.rphase  = 1;
FitParams.cphases = 1; %1:ReconOptions.cbins;
FitParams.fitSlice = 1:FitParams.Nz;


% -- Fitting option
% --- for T1T2 1FA measurements:
% --- 0 - fit T1/T2/BIR/BT2/Beta;  1 - constant BT2;  2: BT2 = BIR; 
% --- 3 - constant BT2 & Beta;     4 - BT2 = BIR, constant Beta
FitParams.ChooseConst = 3;

% -- Do fitting
FitResult = paramFit(FitParams);
 

% -- Check fitted curves
curve = plotFitResult(FitResult,1,1);   

%%
% -- Display color maps and save as gif
T1map_color = displayMap(FitResult.T1map,warmmetal(256),[0 FitResult.maxT1],20,'T1map_const4.gif', FitResult.fitParams.filePath);
T2map_color = displayMap(FitResult.T2map,warmmetal(256),[0 FitResult.maxT2],20,'T2map_const4.gif',FitResult.fitParams.filePath);  


% -- Save as Dicom files
genDicomfromSiemens(abs(round(1000*FitResult.T1map)), TwixObj, ReconOptions.interpFactor, 0, 'T1map_const4');
genDicomfromSiemens(abs(round(1000*FitResult.T2map)), TwixObj, ReconOptions.interpFactor, 0, 'T2map_const4');


%% Save fit results

fitResultFileString = [FitParams.filePath '/multitasking_FitResult_sliceprofile.mat'];
save(fitResultFileString,'FitResult','-v7.3'); 
disp('FitResult Saved.')

dbclear all
