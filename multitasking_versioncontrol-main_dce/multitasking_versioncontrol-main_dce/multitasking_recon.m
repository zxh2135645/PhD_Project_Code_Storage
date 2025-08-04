% MR Multitasking reconstruction code v0.1 - VB/VE/XA multitasking sequence support
% This is the unified branch, composed by Hsu-Lei Lee
% It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories
% T1/T2/T1rho/VFA contrast
% This scipt is mainly for Human Study

%% Set up paths
mainpath=pwd;                       % Put your path here.
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
% If cufinufft keeps failing, try open matlab by running
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


%% Configure recon and load data 

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
ReconOptions.flagUseGPU = false;  % (exist('gpuNUFFT','file')>1); % true if using GPU to do nufft. false if using finufft/MIRT on CPU.
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

 Params.ScanType       = 'IR';     % will be automatically determined in loadData()   
                                      % options: {'IR','IR_VFA','T2IR','T2IR_VFA','SR','Cine'}
% Params.numIRns        = 1;
% Params.numSRns        = 0;
% Params.numT2prep      = 4;
% Params.T2prepDuration = [0 30 55 80];
% Params.numT1rhoPrep   = 0;
% Params.T1rhoDuration  = [];
% Params.numFA          = 1;
% for n = 1:Params.numFA
% 	Params.flipAngleArray(n) = [TwixObj.hdr.MeasYaps.adFlipAngleDegree{n}];
% end
% Params.moduleLength   = lcm((Params.numIRns+Params.numSRns+Params.numT2prep+Params.numT1rhoPrep),Params.numFA);

% Params.SGBlock        = 2;
% Params.linesPerShot   = 710;
% Params.MBfactor       = 1;
% Params.navShift       = 1;
% Params.isTinyGoldenAngle = true;
% Params.isVTR          = false;


%% Pre-processing
% -- Load data into variable struct
ReconOptions.totalTime = inf;%inf;   % length of data to use, in seconds. "inf" to use all.
[Params, ReconOptions, DataArray] = loadData(Params, ReconOptions, TwixObj, DataArray);

% -- Setup trajectories
% -- Crop recon FOV in radial case
ReconOptions.cropFOVy = 1;              % FOVynew = FOVy/cropFOVy for radial trajectory
ReconOptions.cropFOVx = 1;              % FOVxnew = FOVx/cropFOVy for radial trajectory
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


%% Coil compression

% -- Coil compression params
ReconOptions.minNewCoils = 12;          % # of coils to keep after compression
ReconOptions.mincoilEnergy = 0.99;      % percentage of total energy to keep after compression
% -- New coil number 
ReconOptions.newCoils = 12;%min(Params.Ncoils,max(find(DataArray.coilEnergycumsum>ReconOptions.mincoilEnergy,1),ReconOptions.minNewCoils)); 
selectROI = 1;
% -- Do coil compression
[Params, ReconOptions, DataArray] = coilCompression(Params, ReconOptions, DataArray);


temp = sqrt(sum(abs(DataArray.fbp(floor(Params.Ny/2)-floor(Params.Nydisp/2) + (1:Params.Nydisp), floor(Params.Nx/2)-floor(Params.Nxdisp/2) + (1:Params.Nxdisp),:,:)).^2,4));
cw = prctile(temp,95,'all');
%implay(temp./cw);
implay(temp./cw.*DataArray.roi_weighting(floor(Params.Ny/2)-floor(Params.Nydisp/2) + (1:Params.Nydisp), floor(Params.Nx/2)-floor(Params.Nxdisp/2) + (1:Params.Nxdisp),:));


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

%%
% -- Generate relaxation subspace curves
% FitParams.minT1 = 500e-3;
% FitParams.maxT1 = 2000e-3;
% FitParams.minT2 = 30e-3;
% FitParams.maxT2 = 100e-3;
% FitParams.minT1rho = 30e-3;
% FitParams.maxT1rho = 100e-3;
FitParams.minT1 = 100e-3;
FitParams.maxT1 = 2000e-3;
[TemporalBasis,FitParams] = genBlochSubspace(Params, FitParams, ReconOptions, TemporalBasis);


%% Binning: Respiratory

disp('Start binning...')

% -- Set # of respiratory bins
ReconOptions.rbins = 4;    

% -- Breathing rate filter range (cycle/min)
ReconOptions.BRlow  = 3;
ReconOptions.BRhigh = 30;

% -- Do respiratory binning
% = binningResp_PMU(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff, false, true);

% -- Swap R-bins to make end-of-respiration bin #1 if necessary
% DataArray = flipRidx(DataArray,true);

% -- Motion correction (reference: Ridx=1)
%[Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% -- Do respiratory binning
DataArray = binningResp_edge2(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);
%[DataArray,TemporalBasis] = binningResp_Yang(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,TwixObj,true);
%DataArray = binningResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);
% -- Swap R-bins to make end-of-respiration bin #1 if necessary
%DataArray = flipRidx(DataArray,true);

% -- Motion correction (reference: Ridx=1)
%[Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = mocoResp(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
implayZoom(DataArray.binsCard{1});
%% Binning: Cardiac

% -- Set # of cardiac bins
ReconOptions.cbins = 10; 

% -- Heart rate filter range (cycle/min)
ReconOptions.HRlow  = 50;    
ReconOptions.HRhigh = 140;
[DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);

% % -- Do cardiac binning
% if ReconOptions.totalTime > 30
%     % -- Change the last input to true to manually draw heart ROI
%     %[DataArray,TemporalBasis] = binningCard_old(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);
%     [DataArray,TemporalBasis] = binningCard(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);
% else
%     % -- Another cardiac binning method assuming minimum heart rate variation
%     [DataArray,TemporalBasis] = binningCardSelfGate(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);
% end
% 
% % -- Shift end-of-diastole phase to the first bin
% DataArray.diastoleIdx = 11;     %str2double(inputdlg({'Select end-of-diastole phase:'},'User Input',[1 35],{'1'}));
% DataArray = shiftHidx(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);


%% Tensor subspace estimation

ReconOptions.lowmem = 0;                    % 0: no reduction; 1: reduce to nav lines; 2: reduce to cL
ReconOptions.L_tensor = 32;              %try high rank with 1 ngd dur   % Rank for tensor subspace
ReconOptions.L = ReconOptions.L_tensor;

ReconOptions.flagUseInitialGuess = false;   
ReconOptions.flagUseToeplitz = false;
ReconOptions.flagUseGPUToeplitz = true; 

ReconOptions.Ngd_dur = 22;%32
Params.ScanType       = 'IR';
ReconOptions.flagAutoLambdaTensor = false;  % automatically determines lr and sms?
[ReconOptions,DataArray,TempStructTensor] = setTensorParams(Params,ReconOptions,DataArray,TemporalBasis);

%
% -- Overwrite lr and sms if necessary
% ReconOptions.tensor.lr  = 7e-10;
% ReconOptions.tensor.sms = [3e-6 1e-6];

% -- Tensor subspace
[DataArray,TemporalBasis] = tensorSubspace(TempStructTensor,Params,ReconOptions,DataArray,TemporalBasis);
%
ReconOptions.ls_lowmem = 1;
% -- Least-squares tensor recon
[ReconOptions,DataArray,TemporalBasis,SpatialCoeff] = reconLeastSquares(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,true);
tensorDisplay
%-- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
TemporalBasis.Phi_ten=[];
TemporalBasis.Phi_ten=TemporalBasis.Phi;
%% Mean signal comparison between realtime & tensor

[realtime_sig, temp_LV]= realtime_mean_sig(Params, SpatialCoeff, TemporalBasis, ReconstructedImages);
%%
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);
%
cw = prctile(abs(ReconstructedImages.reconRealtime(:)),99);
realtime_sig = realify(realtime_sig);
ReconstructedImages.Realtime_mean_lv = realtime_sig*cw;
ReconstructedImages.ROILV = temp_LV;
% figure; hold on
% plot(ReconstructedImages.Realtime_mean_lv./abs(min(ReconstructedImages.Realtime_mean_lv(:))),'k','LineWidth',2.5,'DisplayName','LV Realtime mean signal');
ROILV = find(ReconstructedImages.ROILV);
ReconstructedImages.tensorRecovery = ReconstructedImages.reconRecovery;
tensorrecon = reshape(ReconstructedImages.tensorRecovery,size(ReconstructedImages.tensorRecovery,1)*size(ReconstructedImages.tensorRecovery,2),[]);
tensorrecon = tensorrecon(ROILV(:),:);
tensorrecon = mean(tensorrecon,1);
ReconstructedImages.cw = prctile(abs(ReconstructedImages.tensorRecovery(:)),99);
tensorrecon = tensorrecon.*ReconstructedImages.cw;
ReconstructedImages.tensorrecon_lv = realify(tensorrecon);


figure; hold on
plot(ReconstructedImages.tensorrecon_lv./ReconstructedImages.tensorrecon_lv(end),'k','LineWidth',2.5,'DisplayName','LV Tensor Recon mean signal');
plot(ReconstructedImages.Realtime_mean_lv./max(ReconstructedImages.Realtime_mean_lv(:)),'r','LineWidth',2.5,'DisplayName','LV Realtime mean signal');
legend;

% figure; hold on
% plot(ReconstructedImages.tensorrecon_lv./abs(ReconstructedImages.tensorrecon_lv(1)),'k','LineWidth',2.5,'DisplayName','LV Recon signal');
% legend;
%% Wavelet recon
close all hidden
ReconOptions.flagAutoLambdaWavelet = false;     % automatically determines lambda?
                                                % (doesn't work well)
% -- Overwrite lambda & alpha if necessary
ReconOptions.wavelet.lambda = 2.5808e-11;
ReconOptions.wavelet.alpha  = 1;

% -- Do wavelet recon
ReconOptions.wavelet.flagContinue = 0;          % 1: use current U as initial value; 0: use U_tensor
ReconOptions.wavelet.maxIter      = 5;         % maximum number of iterations
[ReconOptions,SpatialCoeff] = waveletRecon(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff);

% -- Show tensor images
DispParams.scale = 1;
DispParams.tIdx = [];
DispParams.cIdx = 1;
DispParams.rIdx = 1;
ReconstructedImages = displayTensor(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages,DispParams.scale,DispParams.tIdx,DispParams.cIdx,DispParams.rIdx,[],ReconOptions.dispSlice);

% ---how black-blood images
% ReconstructedImages = displayTensorBlackBlood(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,ReconstructedImages);


%% Export cine images as dicom

% -- interpolate the image if necessary
% [ReconstructedImages.reconCard,ReconOptions.interpFactor] = interp2dimage(abs(ReconstructedImages.reconCard),[0.9375 0.9375],'resolution',TwixObj);

% -- Write dicom files
genDicomfromSiemens(ReconstructedImages.reconCard, TwixObj, ReconOptions.interpFactor, 1, 'Cine');


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
close all hidden
% -- Set the time points in the recovery curve to be used in fitting
FitParams.imstart = 1;
FitParams.imskip  = 10;

% -- 0: k-means; 1: force ROI in central FOV/4; 2: force ROI in central FOV/2 
FitParams.flagLargeROI = 0;

% -- which respiratory/cardiac phases/slice to fit?
FitParams.rphase  = 1;
FitParams.cphases = 1; %1:ReconOptions.cbins;
FitParams.fitSlice = 1:FitParams.Nz;

% -- Flip angle adjustment
FitParams.initGlobalBeta = 0.1;
%FitParams.B1_coeff = 1;

% -- initial values
FitParams.initT1     = 1;
FitParams.initT2     = 40e-3;
FitParams.initBIR    = 1;
FitParams.initBT2    = 1;

% -- lower & upper bounds
FitParams.minT1 = 100e-3;
FitParams.maxT1 = 3000e-3;
FitParams.minT2 = 20e-3;
FitParams.maxT2 = 200e-3;
FitParams.minBIR = 1;
FitParams.maxBIR = 1;
FitParams.minBeta = FitParams.initGlobalBeta;
FitParams.maxBeta = FitParams.initGlobalBeta;
FitParams.minBT2 = 1;
FitParams.maxBT2 = 1;


% -- Fitting option
% --- for T1T2 1FA measurements:
% --- 0 - fit T1/T2/BIR/BT2/Beta;  1 - constant BT2;  2: BT2 = BIR; 
% --- 3 - constant BT2 & Beta;     4 - BT2 = BIR, constant Beta
FitParams.ChooseConst = 0;
FitParams.flag2stepFitting = 1;
% -- Do fitting
FitResult = paramFit(FitParams);

%% 

% -- Check fitted curves
curve = plotSelectedVoxel(interp2dimage(ReconstructedImages.reconRecovery,FitResult.fitParams.interpFactor,'factor'),[],[],[],FitResult,1);   

% -- Display color maps and save as gif
T1map_color = displayMap(FitResult.T1map,warmmetal(256),[0 3],20,'T1map_const4.gif', FitResult.fitParams.filePath);
T2map_color = displayMap(FitResult.T2map,warmmetal(256),[0 0.15],20,'T2map_const4.gif',FitResult.fitParams.filePath);  


% -- Save as Dicom files
genDicomfromSiemens(abs(round(1000*FitResult.T1map)), TwixObj, ReconOptions.interpFactor, 0, 'T1map_const4');
genDicomfromSiemens(abs(round(1000*FitResult.T2map)), TwixObj, ReconOptions.interpFactor, 0, 'T2map_const4');

figure; imagesc(FitResult.T1map*1000,[0 3000]),axis image,colormap(warmmetal); colorbar; title('T1 (ms)'); hold on
f = gcf;
set(f,'Position',get(0,'screensize'));
exportgraphics(f,[erase(Params.fileString,'multitasking'),'T1map.tiff'],'Resolution',300);
close all hidden;
figure; imagesc(FitResult.T2map*1000,[0 150]),axis image,colormap(warmmetal); colorbar; title('T2 (ms)'); hold on
f = gcf;
set(f,'Position',get(0,'screensize'));
exportgraphics(f,[erase(Params.fileString,'multitasking'),'T2map.tiff'],'Resolution',300);
close all hidden;
figure; imagesc(FitResult.T2starmap*1000,[0 150]),axis image,colormap(warmmetal); colorbar; title('T2* (ms)'); hold on
f = gcf;
set(f,'Position',get(0,'screensize'));
exportgraphics(f,[erase(Params.fileString,'multitasking'),'T2starmap.tiff'],'Resolution',300);
close all hidden;
%% Save fit results

fitResultFileString = [FitParams.filePath '/multitasking_FitResult_sliceprofile.mat'];
save(fitResultFileString,'FitResult','-v7.3'); 
disp('FitResult Saved.')

dbclear all
