function [params, reconOptions, dataArray] = loadData(params, reconOptions, twixObj, dataArray)

extractVarFromStruct(params);

if ~isfield(params,'navShift')
    navShift = 1;
end

if ~isfield(params,'SGBlock')
    SGBlock = 1;
end

if SGBlock > 1
    isMultitasking = true;
else
    isMultitasking = false;
end

vec = @(x) x(:);

try

%% ScanType & moduleLength

if ~exist('isT2IR','var')
    isT2IR = 1;
end
if ~exist('isT1rhoIR','var')
    isT1rhoIR = 1;
end

if ~strcmp(params.ScanType,'CEST')    %~isfield(params,'ScanType')
    if numSRns == 0 && numIRns == 0 && numT2prep == 0 && numT1rhoPrep == 0
        params.ScanType = 'Cine';
    elseif numIRns > 0 && numT2prep == 0 && numT1rhoPrep == 0 && numFA == 1
        params.ScanType = 'IR';
    elseif numIRns > 0 && numT2prep == 0 && numT1rhoPrep == 0 && numFA > 1
        params.ScanType = 'IR_VFA';
    elseif numT2prep > 0 && numT1rhoPrep == 0 && numFA == 1
        if isT2IR
            params.ScanType = 'T2IR';
        else
            params.ScanType = 'T2prep';
        end
    elseif numT2prep > 0 && numT1rhoPrep == 0 && numFA > 1
        if isT2IR
            params.ScanType = 'T2IR_VFA';
        else
            params.ScanType = 'T2prep_VFA';
        end
    elseif numT2prep == 0 && numT1rhoPrep > 0 && numFA == 1
        if isT1rhoIR
            params.ScanType = 'T1rhoIR';
        else
            params.ScanType = 'T1rho';
        end
    elseif numT2prep == 0 && numT1rhoPrep > 0 && numFA > 1
        if isT1rhoIR
            params.ScanType = 'T1rhoIR_VFA';
        else
            params.ScanType = 'T1rho_VFA';
        end
    elseif numT2prep > 0 && numT1rhoPrep > 0
        params.ScanType = 'T1rhoT2IR';
    elseif numSRns > 0
        params.ScanType = 'SR';
    else
        fprintf('ScanType not supported\n');
    end 
end
fprintf('ScanType: %s. Loading k-space data... ', params.ScanType);

if strcmp(params.ScanType,'CEST')
    moduleLength = params.CESTNumRep;
else
    try
        moduleLength = lcm((numIRns+numSRns+numT2prep+numT1rhoPrep),numFA);
    catch
        moduleLength = numFA;
    end
end

%% re-arrange data array and encoding info

if isempty(linesPerShot) || linesPerShot == 0
    linesPerShot = params.lSegments;
end

if Necho == 1
    isVTR = false;
end

% Retrive raw data
kspaceData = twixObj.image.unsorted();        % NCol, NCha, Nsamples in acquisition order
%kspaceData = reshape(kspaceData, size(kspaceData,1), size(kspaceData,2), Necho, []);
kspaceData = permute(kspaceData,[3 1 4 2]);
Nsamples = size(kspaceData,1);

% physio recordings
try
    physio = twixObj.PMUdata.physio_all.';
catch
    physio = zeros(Nsamples,7);
end
if size(physio,2) < 7
    physio(:,end+1:7) = 0;
end

% Convert from hardware coordinate to LPS patient coordinate
kspaceData = conj(kspaceData);

% ACQ segments
if isMultitasking && isVTR
    SGBlockACQCount = 1 + (SGBlock-1)*Necho;
    ACQCount = 1;
    for n = 2:SGBlock
        ACQCount = [ACQCount; Necho];
    end
    ACQCount = vec(repmat(ACQCount,[1 ceil(Nsamples/SGBlockACQCount)]));
    ACQCount = cumsum(ACQCount);
    lineCount = find(ACQCount>=Nsamples,1);
    shotCount = floor(lineCount/linesPerShot);
else
    SGBlockACQCount = SGBlock*Necho;
    lineCount = Nsamples/Necho;
    ACQCount = (1:lineCount)*Necho;
    shotCount = floor(lineCount/linesPerShot);
end

% data cut-off 
if ~isMultitasking
    cutoff_shot = 0;
    cutoff = 0;
elseif strcmp(params.ScanType,'CEST') 
    cutoff_shot = 0;
    cutoff = 0; 
elseif strcmp(params.ScanType,'Cine')
    cutoff_shot = ceil(2/alTR_seconds);
    cutoff = ACQCount(linesPerShot*cutoff_shot); 
else
    cutoff_shot = 70;
    cutoff = ACQCount(linesPerShot*cutoff_shot); 
end

cutoff_end_shot = min(cutoff_shot + ceil(reconOptions.totalTime/params.alTR_seconds/moduleLength)*moduleLength, shotCount);
cutoff_end      = ACQCount(cutoff_end_shot*linesPerShot);
% cutoff_end = min(cutoff + ceil(reconOptions.totalTime/params.alTR_seconds/moduleLength)*moduleLength*linesPerShot, cutoff + floor((size(kspaceData,1)-cutoff)/linesPerShot/moduleLength)*moduleLength*linesPerShot);
cutoff_shift_line = SGBlock - mod(cutoff-1,SGBlock);
cutoff_shift      = (cutoff_shift_line-1)*Necho + navShift;   % Segments may not be a multiple of SGBlock
lineCount = (cutoff_end_shot - cutoff_shot)*linesPerShot;

% Crop k-space data
kspaceData = single(kspaceData((cutoff+1):cutoff_end,:,:,:,:));
physio = physio((cutoff+1):cutoff_end,:);

% SG line indices
navIndices = cutoff_shift_line:SGBlock:lineCount;
if isVTR
    navIndices_full = cutoff_shift:SGBlockACQCount:size(kspaceData,1);
else
    navIndices_full = zeros(Necho,numel(cutoff_shift_line:SGBlock:lineCount));
    for n = 1:Necho
        navIndices_full(n,:) = (cutoff_shift+(n-1)):SGBlockACQCount:size(kspaceData,1);
    end
    navIndices_full = vec(navIndices_full)';
end

% ACQ timing for navData and kspaceData 
if strcmp(params.ScanType,'CEST') 
    ACQTiming = [(CESTSatDuration+CESTSatSpoilerTotalDuration)*1e-3 ones(1,linesPerShot)*lEchoSpacing]; % Sat module + one shot
    ACQTiming = cumsum(repmat(ACQTiming,[1 ceil(size(kspaceData,1)/SGBlock)]));
    ACQTiming(1:(linesPerShot+1):end) = [];       % remove Sat module time points
    ACQTiming = ACQTiming((cutoff+1):cutoff_end) - ACQTiming(cutoff+1);
    navTiming = ACQTiming(navIndices);
elseif isMultitasking && isVTR
    navEchoSpacing = lEchoSpacing - alTE_seconds(end) + alTE_seconds(1);
    ACQTiming = [0 ones(1,numel(0:cutoff_shift_line-2))*lEchoSpacing];
    ACQTiming = cumsum([ACQTiming repmat([navEchoSpacing ones(1,SGBlock-1)*lEchoSpacing],[1 ceil(lineCount/SGBlock)])]);
    ACQTiming = ACQTiming(1:lineCount);
    navTiming = ACQTiming(navIndices);
%     tEnd = 0;
%     ACQTiming_pre = zeros((cutoff_shift_line-1)*Necho,1);
%     for n = 1:cutoff_shift_line-1
%         ACQTiming_pre((1:Necho)+(n-1)*Necho) = tEnd*ones(Necho,1);
%         tEnd = tEnd + lEchoSpacing;
%     end
% 
%     ACQTiming = zeros(SGBlockACQCount,ceil(lineCount/SGBlock));
%     for m = 1:ceil(lineCount/SGBlock)
%         ACQTiming(1,m) = tEnd;
%         tEnd = tEnd + navEchoSpacing;
%         for n = 2:SGBlock
%             ACQTiming((1:Necho)+Necho*(n-2)+1,m) = tEnd*ones(Necho,1);
%             tEnd = tEnd + lEchoSpacing;
%         end
%     end
% 
%     ACQTiming    = [vec(ACQTiming_pre);vec(ACQTiming)];
%     ACQTiming_full = ACQTiming(1:size(kspaceData,1));
%     ACQTiming(navIndices_full) = [];
%     navTiming      = ACQTiming(navIndices_full);
else
    ACQTiming = (1:size(kspaceData(1:Necho:end,:),1))*lEchoSpacing;
    navTiming = ACQTiming(navIndices);
%     ACQTiming_full(n,:) = vec(repmat(n:Necho:size(kspaceData,1),[Necho,1]))';
%     navTiming = ACQTiming_full(navIndices_full);
    
end

% Phase-encoding order
try
    if strcmp(twixObj.image.softwareVersion,'vb')
        sWipMemBlock = twixObj.hdr.MeasYaps.sWiPMemBlock;
    else
        sWipMemBlock = twixObj.hdr.MeasYaps.sWipMemBlock;
    end
    if isCartesian || strcmp(Trajectory,'Spiral') || sWipMemBlock.alFree{64} >= 230831
        linOrder_full = twixObj.image.Lin((cutoff+1):cutoff_end);
    else
        linOrder_full = zeros(1,cutoff_end-cutoff);
    end
catch
    if isCartesian
        linOrder_full = twixObj.image.Lin((cutoff+1):cutoff_end);
    else
        linOrder_full = zeros(1,cutoff_end-cutoff);
    end
end
parOrder_full = twixObj.image.Par((cutoff+1):cutoff_end);
physio_k = physio;

% Collect navData
if isMultitasking
    navData = single(kspaceData(navIndices_full,:,:,:,:));
    physio_nav = physio(navIndices_full,:);
    if ~isVTR
        navData = reshape(navData,Necho,[],size(navData,2),size(navData,3),size(navData,4));
        navData = permute(navData, [2 3 4 5 1]);
        physio_nav = permute(reshape(physio_nav,Necho,[],size(physio,2)),[2 3 1]);
    end
    % remove navData from k-spaceData
    kspaceData(navIndices_full,:,:,:,:) = [];
    linOrder_full(navIndices_full) = [];
    parOrder_full(navIndices_full) = [];
    physio_k(navIndices_full,:) = [];
else
    navData = [];
    physio_nav = physio(navIndices_full,:);
end

kspaceData = reshape(kspaceData,Necho,[],size(kspaceData,2),size(kspaceData,3),size(kspaceData,4));
kspaceData = permute(kspaceData,[2 3 4 5 1]);
linOrder = linOrder_full(1:Necho:end);
parOrder = parOrder_full(1:Necho:end);
physio_k = permute(reshape(physio_k,Necho,[],size(physio,2)),[2,3,1]);

% no need to flip even echoes, Siemens rawdata are already flipped
% if Necho > 1  
%     navData(:,:,:,:,2:2:end)     = -navData(:,:,:,:,2:2:end);
%     navData(:,2:end,:,:,2:2:end) = flip(navData(:,2:end,:,:,2:2:end),2);
%     kspaceData(:,:,:,:,2:2:end)     = -kspaceData(:,:,:,:,2:2:end);
%     kspaceData(:,2:end,:,:,2:2:end) = flip(kspaceData(:,2:end,:,:,2:2:end),2);
% end

% Data dimensions
Ntpoint     = lineCount;
Nread       = size(kspaceData,1);
Nkx         = size(kspaceData,2);
Ncoils      = size(kspaceData,4);
totalTime   = Ntpoint/linesPerShot*alTR_seconds;
fprintf('done. Data length = %.2f sec\n', totalTime);

Norig = lBaseResolution;
Nz    = scanParams.NoOfFourierPartitions;
% % If no oversampling on the Z direction, Nzorig = Nz
% try 
%     Nzorig = scanParams.NImagePar;
% catch
%     % Nzorig = scanParams.NoImagesPerSlab;
%     Nzorig = scanParams.NPar;
% end

% if SMS, reset partition encoding array and Nz
if MBfactor > 1
    for n = 1:MBfactor
        parOrder(n:MBfactor:end) = mod(mod((cutoff - ceil(cutoff/SGBlock)), MBfactor) + n - 1, MBfactor) + 1;
    end
    Nz = max(parOrder);
    parOrder = Nz - parOrder + 1;
    Nzorig = Nz;
end
Nzorig = min(Nzorig,Nz);

ovsSL = Nz - Nzorig;                                                    % Slice oversampling extension
%N     = scanParams.NImageCols*twixObj.hdr.Dicom.flReadoutOSFactor;     % oversample
N     = Norig*twixObj.hdr.Dicom.flReadoutOSFactor;                      % oversample
ovsRO = N - Norig;                                                      % FOV extension
DC    = Nkx - ceil(N/2) + 1;
DC_kx = floor(N/2) + 1;
DC_ky = twixObj.image.Lin(1);
    
if isMultitasking
    DC_kz = twixObj.image.Par(1);
else
    DC_kz = floor(Nz/2) + 1;
end

if strcmp(Trajectory,'Spiral')
    DC    = 1;
    DC_kx = 1;
end

if N ~= Nkx
    isAsymmEcho = true;
else
    isAsymmEcho = false;
end

if isCartesian
    Nx = N;
    %DC_ky = floor(scanParams.NoOfFourierLines/2) + 1;
    Ny = scanParams.NoOfFourierLines;
%     Nydisp = twixObj.hdr.Config.NImageLins;
%     Nxdisp = twixObj.hdr.Config.NImageCols;
    Nxdisp = Norig;
    Nydisp = round(Norig*dPhaseFOV_mm/dReadoutFOV_mm/2)*2;
else
    if params.Ntrajs
        Ntrajs = params.Ntrajs;
    else
        Ntrajs = 0;
    end
    if params.trajInc
        trajInc = params.trajInc;
    else
        trajInc = 0;
    end
    Ny = N;
    Nx = N;
    Nydisp = Norig;
    Nxdisp = Norig;
end

Nxdisp = min(Nxdisp,Nx);
Nydisp = min(Nydisp,Ny);
Nzdisp = Nzorig;

%% Output struct variables

params.cutoff       = cutoff;
params.cutoff_shot  = cutoff_shot;
params.cutoff_shift = cutoff_shift;
params.cutoff_shift_line = cutoff_shift_line;
params.N      = N;
params.Norig  = Norig;
params.Nzorig = Nzorig;
params.Nx = Nx;
params.Ny = Ny;
params.Nz = Nz;
params.Nxdisp = Nxdisp;
params.Nydisp = Nydisp;
params.Nzdisp = Nzdisp;
params.DC    = DC;
params.DC_kx = DC_kx;
params.DC_ky = DC_ky;
params.DC_kz = DC_kz;
params.ovsRO = ovsRO;
params.ovsSL = ovsSL;
params.Ntrajs  = Ntrajs;
params.trajInc = trajInc;
params.Ntpoint    = Ntpoint;
params.Nread      = Nread;
params.Nkx        = Nkx;
params.Ncoils     = Ncoils;
params.Necho      = Necho;
params.linesPerShot = linesPerShot;
params.moduleLength = moduleLength;
params.isAsymmEcho = isAsymmEcho;
params.isVTR       = isVTR;
params.isMultitasking = isMultitasking;

reconOptions.totalTime = totalTime;
reconOptions.interpFactor = [1 1];

dataArray.physio.all.EKG   = physio(:,1,:);
dataArray.physio.all.EKG2  = physio(:,2,:);
dataArray.physio.all.EKG3  = physio(:,3,:);
dataArray.physio.all.PULSE = physio(:,4,:);
dataArray.physio.all.EXT   = physio(:,5,:);
dataArray.physio.all.RESP  = physio(:,6,:);
dataArray.physio.all.TriggerDelay  = physio(:,7,:);
dataArray.physio.nav.EKG   = physio_nav(:,1,:);
dataArray.physio.nav.EKG2  = physio_nav(:,2,:);
dataArray.physio.nav.EKG3  = physio_nav(:,3,:);
dataArray.physio.nav.PULSE = physio_nav(:,4,:);
dataArray.physio.nav.EXT   = physio_nav(:,5,:);
dataArray.physio.nav.RESP  = physio_nav(:,6,:);
dataArray.physio.nav.TriggerDelay = physio_nav(:,7,:);
dataArray.physio.kspace.EKG   = physio_k(:,1,:);
dataArray.physio.kspace.EKG2  = physio_k(:,2,:);
dataArray.physio.kspace.EKG3  = physio_k(:,3,:);
dataArray.physio.kspace.PULSE = physio_k(:,4,:);
dataArray.physio.kspace.EXT   = physio_k(:,5,:);
dataArray.physio.kspace.RESP  = physio_k(:,6,:);
dataArray.physio.kspace.TriggerDelay  = physio_k(:,7,:);

dataArray.flagIsTrajectorySet  = false;
dataArray.flagIsPrewhitened    = false;
dataArray.flagIsCompressed     = false;
dataArray.flagIsDriftCorrected = false;
dataArray.navData    = navData;
dataArray.kspaceData = kspaceData;
dataArray.navIndices = navIndices;
dataArray.navIndices_full = navIndices_full;
dataArray.ACQTiming  = ACQTiming;
dataArray.navTiming  = navTiming;
dataArray.linOrder_full = linOrder_full(:);
dataArray.parOrder_full = parOrder_full(:);
dataArray.linOrder = linOrder(:);
dataArray.parOrder = parOrder(:);

if isfield(dataArray,'lrw')
   dataArray = rmfield(dataArray,'lrw');
end

catch ME
    rethrow(ME)
end


