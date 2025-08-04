function [reconOptions,dataArray,tempStructTensor] = setTensorParams(params,reconOptions,dataArray,temporalBasis)

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end

fprintf('Setup tensor subspace parameters...')
tic;

vec = @(x) x(:);
row = @(x) x(:).';

tempL = 72;
flagDataDriven = false;
flagUseFirstEcho = false;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);

if ~exist('mixer', 'var')
    mixer = eye(Ncoils);
end

if size(navData,5) == 1
    flagUseFirstEcho = true;
end
if flagUseFirstEcho
    navData = navData(:,:,:,:,1);
end
tempNecho = size(navData,5);

if strcmp(ScanType,'CEST') && size(curvePhi,1) == linesPerShot*CESTNumRep 
    flagDataDriven = true;
elseif strcmp(ScanType,'CEST') && size(curvePhi,1) == linesPerShot*CESTNumRep*CESTnumFreqOffsets 
    flagDataDriven = false;
elseif strcmp(ScanType,'CEST')
    curvePhi = [];
    fprintf(2,'Wrong curvePhi size.\n');
elseif moduleLength > 1 && size(curvePhi,1) == linesPerShot
    flagDataDriven = true;
elseif moduleLength > 1 && size(curvePhi,1) == linesPerShot*moduleLength 
    flagDataDriven = false;
elseif moduleLength > 1
    curvePhi = [];
    fprintf(2,'Wrong curvePhi size.\n');
end

rbins = max(Ridx);
cbins = max(Hidx);
reconOptions.rbins = rbins;
reconOptions.cbins = cbins;

% setup default lambda
switch ScanType
    case 'Cine'
        lr = 1e-4;  % how low-rank?
        sms = [0 0];
        smorders = [];
        circs = [];
    case 'SR'
        lr = 1e-4; 
        if rbins > 1 && cbins == 1
            sms = [5e-6 1e-7];    % [resp, wall-clock]
            smorders = [1 1];     % piecewise constant
            circs = [false false];
        elseif cbins > 1
            sms = [5e-6 5e-6 5e-6];   % [cardiac, resp, wall-clock]
            smorders = [1 1 1];       % piecewise constant
            circs = [true false false];
            if rbins == 1; sms(2) = 0; end
        else
            sms = 1e-9;
            smorders = 1;
            circs = false;
        end
    case 'IR'
        lr = 1e-5;%1e-4;
        if rbins > 1 && cbins == 1
            sms = [3e-6 1e-6];  % [cardiac, resp], high->smooth, low->more contrast change
            smorders = [1 1];
            circs = [true false];
            if rbins == 1; sms(2) = 0; end
            if cbins == 1; sms(1) = 0; end
        elseif cbins > 1
            sms = [1e-6, 1e-6, 1e-5];%[1e-6 1e-6 1e-5]; % [1e-10 1e-10 1e-4] time-domain notok%[5e-6 5e-6 5e-6];   % [cardiac, resp, wall-clock]
            smorders = [1 1 1];       % piecewise constant
            circs = [true false false];
            if rbins == 1; sms(2) = 0; end
        else
            sms = 1e-9;
            smorders = 1;
            circs = false;
        end
    case 'T2prep'
        lr = 1e-4;
        sms = [3e-6 1e-6 3e-6];
        smorders = [1 1 1];
        circs = [true false false];
        if rbins == 1; sms(2) = 0; end
        if cbins == 1; sms(1) = 0; end
    otherwise
        lr = 1e-4; 
        sms = [1e-6 1e-6];  % [cardiac, resp], high->smooth, low->more contrast change
        smorders = [1 1];
        circs = [true false];
        if rbins == 1; sms(2) = 0; end
        if cbins == 1; sms(1) = 0; end
end

% setup wall clock
switch ScanType
    case 'CEST'
        NDCEAve = 1;
        [dataArray.CESTSatFreqOffsetppmListUnique,dataArray.CESTSatFreqIdx,wallClock] = unique(CESTSatFreqOffsetppmList);
%         wallClock(curveMaskNav==0) = 0;
        wallClock = row(repmat(row(wallClock),[linesPerShot*CESTNumRep/SGBlock 1]));
%         wallClock_full = wallClock;
        isWallClockSet = true;
        navData_temp = reshape(reshape(permute(navData,[5 1 2 3 4]),[],Ncoils)*mixer(:,1:max(ceil(Ncoils/3),4)),size(navData,1)*tempNecho,[]);
    case 'SR'
        if ~exist('NDCEAve','var')
            NDCEAve = 2;
        end
        wallClock = ceil((1:Ntpoint)/linesPerShot/NDCEAve);
        wallClock = wallClock(navIndices);
        if ~flagUseFirstEcho
            wallClock = row(repmat(wallClock,[tempNecho 1]));
        end
%         wallClock_full = wallClock;
        isWallClockSet = true;
        navData_temp  = reshape(reshape(permute(navData,[5 1 2 3 4]),[],Ncoils)*mixer(:,1:Ncoils),size(navData,1),[]);
    case 'T2IR'
        NDCEAve = 1;
        if flagDataDriven
%             wallClock_full = row(mod(ceil(navIndices_full/linesPerShot/Necho)-1,moduleLength)+1);
            wallClock = row(mod(ceil(navIndices/linesPerShot)-1,moduleLength)+1);
        else
%             wallClock_full = ones(1,numel(navIndices_full));
            wallClock = ones(1,numel(navIndices));
        end
        if ~flagUseFirstEcho
            wallClock = row(repmat(wallClock,[tempNecho 1]));
        end
        isWallClockSet = true;
        navData_temp = reshape(reshape(permute(navData,[5 1 2 3 4]),[],Ncoils)*mixer(:,1:max(ceil(Ncoils/3),4)),size(navData,1)*tempNecho,[]);
    case 'IR'
        if ~exist('Ngd_dur','var')
            Ngd_dur = 1;
        end

        NDCEAve = 1;
        wallClock = ceil((1:Ntpoint)/linesPerShot/Ngd_dur/NDCEAve);
        wallClock = wallClock(navIndices);

        isWallClockSet = true;
        navData_temp   = reshape(reshape(permute(navData,[5 1 2 3 4]),[],Ncoils)*mixer(:,1:ceil(Ncoils/3)),size(navData,1)*tempNecho,[]);
    otherwise
        NDCEAve = 1;
        wallClock = ones(1,numel(navIndices));
%         wallClock_full = row(repmat(wallClock,[tempNecho 1]));
        if ~flagUseFirstEcho
            wallClock = row(repmat(wallClock,[tempNecho 1]));
        end
        
        isWallClockSet = false;
        navData_temp   = reshape(reshape(permute(navData,[5 1 2 3 4]),[],Ncoils)*mixer(:,1:ceil(Ncoils/3)),size(navData,1)*tempNecho,[]);
end

[UU,Ss,ts_proj] = svde(navData_temp(:,:));
ts_proj = ts_proj(:,1:tempL);
if flagDataDriven
    navData_temp = UU(:,1:tempL)*Ss(1:tempL, 1:tempL);
else
    navData_temp = navData_temp*ts_proj;
end

if strcmp(ScanType,'CEST')
    curvePhi_orig = curvePhi;
    cL = size(curvePhi,2);
    if flagDataDriven
        Nseg_orig = size(curvePhi_orig,1);
    else
        Nseg_orig = size(curvePhi_orig,1)/CESTnumFreqOffsets;
    end
elseif strcmp(ScanType,'Cine')
    curvePhi_orig = [];
    curvePhi = [];
    cL = 1;
    Nseg_orig   = 1;
elseif isempty(curvePhi)
    curvePhi_orig = [];
    Nseg_orig = linesPerShot;
else
    curvePhi_orig = curvePhi;
    cL = size(curvePhi,2);
    Nseg_orig = size(curvePhi_orig,1);
end
Segidx_orig = mod(navIndices-1, Nseg_orig) + 1;

% Hard coded here XZ 11/13/2024
%Segidx_temp = repmat([ceil([1:80]/2), ceil([0:99]/4)+40, ceil([0:159]/8)+65],[1 319]);
%Segidx_orig = (2*Segidx_temp-1);
%navidx_temp = (2*Segidx_temp);

if  lowmem == 0
    Nseg   = Nseg_orig;
    Segidx = Segidx_orig;
    [uni,~,~] = unique(Segidx_orig);
else
    % XZ 11/14/2024 
    [uni1,~,~] = unique(Segidx_orig);
    [uni2,~,~] = unique(navidx_temp);
    uni = [uni1; uni2];
    uni = 4*uni(:);
    Segidx = Segidx_orig;
    Nseg = numel(uni);


    if strcmp(ScanType,'CEST') && ~flagDataDriven
        curvePhi_orig = reshape(curvePhi_orig,Nseg_orig,[],cL);
        curvePhi_orig = curvePhi_orig(:,CESTfreqStart:CESTfreqEnd,:);
        curvePhi = reshape(curvePhi_orig(uni,:,:),[],cL);
        curvePhi_orig = reshape(curvePhi_orig,[],cL);
    elseif ~isempty(curvePhi_orig)
        curvePhi = curvePhi_orig(uni,:);
    end
end

echoIdx = repmat(1:tempNecho,[1 numel(Segidx)]);
Segidx  = row(repmat(row(Segidx),[tempNecho 1]));
HidxNecho = row(repmat(row(Hidx),[tempNecho 1]));
RidxNecho = row(repmat(row(Ridx),[tempNecho 1]));

if ~strcmp(ScanType,'CEST')
    curveMaskNav = ones(size(wallClock));
end
if strcmp(ScanType,'SR') || strcmp(ScanType, 'IR') %XZ 11/20/2024
    % Regular (non-weighted) recon:
    [navData_tensor,mask] = regate(navData_temp,Segidx,HidxNecho,RidxNecho,wallClock,echoIdx);
else
    if max(Hidx) == 1 && max(Ridx) == 1
        [navData_tensor,mask] = regate(navData_temp,Segidx,HidxNecho,RidxNecho,wallClock(:).*curveMaskNav(:),echoIdx);
    else
        % Weighted recon:
        % find bestres for each imaging block
        bestresi = 1000*ones(ceil(Ntpoint/(linesPerShot*moduleLength))*(linesPerShot*moduleLength),1); %can change interpolation method
        if strcmp(ScanType,'Cine')
            bestresi(1:Ntpoint) = interp1(navIndices,bestres(:),1:Ntpoint,'nearest','extrap'); %can change interpolation method
        else
            bestresi(1:Ntpoint) = interp1Segmented(bestres(:),navIndices,linesPerShot,'cols','nearest'); %can change interpolation method
        end
        motionw = reshape(bestresi, Nseg_orig, []);
        motionw = vec(bsxfun(@ldivide,motionw,median(motionw,2)));
        motionw = min(motionw,sqrt(3));
        motionw = motionw(navIndices); %get back weightings at nav indices
        if ~flagUseFirstEcho
            motionw = vec(repmat(row(motionw),[tempNecho 1]));
        end

        tensorw2 = regate(motionw.^2,Segidx,HidxNecho,RidxNecho,wallClock(:).*curveMaskNav(:),echoIdx);

        navData_tensor = regate(bsxfun(@times,motionw.^2,navData_temp),Segidx,HidxNecho,RidxNecho,wallClock(:).*curveMaskNav(:),echoIdx);
        navData_tensor = bsxfun(@rdivide,navData_tensor,tensorw2);
        navData_tensor(~isfinite(navData_tensor)) = 0;
        mask = repmat(sqrt(tensorw2),[size(navData_tensor,1) 1 1 1 1]);
    end
end

if size(navData_tensor,2) < Nseg
    navData_tensor(:,Nseg,:,:,:) = 0;
    mask(:,Nseg,:,:,:) = 0;
end

sizes_orig = [1, 1, 1, 1, 1, 1];
sizes_orig(1:ndims(navData_tensor)) = size(navData_tensor);
sizes_lm = sizes_orig;
sizes_lm(2) = cL;

if strcmp(ScanType,'CEST') 
    if flagDataDriven && sizes_orig(3)*sizes_orig(4)>1 && lowmem == 0
        curvePhitemp = reshape(curvePhi_alt,Nseg_orig,[],cL);
        curvePhitemp = reshape(curvePhitemp(:,CESTfreqStart:CESTfreqEnd,:),[],cL);
        navData_tensor = reshape(permute(navData_tensor,[1 2 5 3 4 6]),sizes_orig(1),sizes_orig(2)*sizes_orig(5),sizes_orig(3),sizes_orig(4),1,[]);
        mask = reshape(permute(mask,[1 2 5 3 4 6]),size(navData_tensor));
    elseif flagDataDriven && sizes_orig(3)*sizes_orig(4)>1 
        curvePhitemp = reshape(curvePhi_alt,Nseg_orig,[],cL);
        curvePhitemp = reshape(curvePhitemp(uni,CESTfreqStart:CESTfreqEnd,:),[],cL);
        navData_tensor = reshape(permute(navData_tensor,[1 2 5 3 4 6]),sizes_orig(1),sizes_orig(2)*sizes_orig(5),sizes_orig(3),sizes_orig(4),1,[]);
        mask = reshape(permute(mask,[1 2 5 3 4 6]),size(navData_tensor));
    elseif flagDataDriven
        curvePhitemp = curvePhi;
    else
        navData_tensor = reshape(permute(navData_tensor,[1 2 5 3 4 6]),sizes_orig(1),sizes_orig(2)*sizes_orig(5),sizes_orig(3),sizes_orig(4),1,[]);
        mask = reshape(permute(mask,[1 2 5 3 4 6]),size(navData_tensor));
        curvePhitemp = curvePhi;
    end
    Nseg = size(mask,2);
else
    curvePhitemp = curvePhi;
end
sizes = [1, 1, 1, 1, 1, 1];
sizes(1:ndims(navData_tensor)) = size(navData_tensor);
sizes_orig(2) = Nseg_orig;
sizes_lm = sizes;
sizes_lm(2) = cL;

morozov = dataArray.msdev^2*norm(double(mask(:)))^2;

doBloch = ~isempty(curvePhi);
if doBloch
    fprintf('generating navData_bloch...')
    
    if strcmp(ScanType,'CEST') 
        navData_bloch = zeros(sizes);
        navData_bloch_cL = zeros(sizes_lm);
        for nW = 1:size(navData_tensor,5)
            for nR = 1:size(navData_tensor,4)
                for nH = 1:size(navData_tensor,3)
                    tempmask = mask(1,:,nH,nR,nW,1);
                    tempnav = reshape(permute(tempmask.^2.*navData_tensor(:,:,nH,nR,nW,:),[1 3 4 5 6 2]),[],Nseg);
                    tempcurvePhi = curvePhitemp.*tempmask.^2';
                    temp = tempnav*pinv(tempcurvePhi).';
                    navData_bloch(:,:,nH,nR,nW,:) = permute(reshape(temp*(curvePhitemp).',sizes([1 6 2])),[1 3 4 5 6 2]);
                    navData_bloch_cL(:,:,nH,nR,nW,:) = permute(reshape(temp,sizes_lm([1 6 2])),[1 3 4 5 6 2]);
                end
            end
        end
    else
        navData_bloch = pcg(@(x)vec((reshape(permute(mask.^2,[1 3 4 5 6 2]),[],Nseg).*(reshape(x,[],cL)*curvePhitemp.'))*curvePhitemp),...
                    vec(reshape(permute(mask.^2.*navData_tensor,[1 3 4 5 6 2]),[],Nseg)*curvePhitemp),[],200);
        if lowmem == 2
            navData_bloch_cL = ipermute(reshape(navData_bloch,sizes_lm([1 3 4 2 5])),[1 3 4 2 5]);
        end
        navData_bloch = ipermute(reshape(reshape(navData_bloch,[],cL)*curvePhitemp.',sizes([1 3 4 5 6 2])),[1 3 4 5 6 2]);
    end
    % morozov_new=morozov-norm((navData_bloch(:)-navData_tensor(:)).*mask(:))^2
else
    navData_bloch = [];
end

if strcmp(ScanType,'CEST') && flagDataDriven && sizes_orig(3)*sizes_orig(4)>1
    navData_tensor = ipermute(reshape(navData_tensor,sizes_orig(1),[],sizes_orig(5),sizes_orig(3),sizes_orig(4),sizes_orig(6)),[1 2 5 3 4 6]);
    navData_bloch  = ipermute(reshape(navData_bloch,sizes_orig(1),[],sizes_orig(5),sizes_orig(3),sizes_orig(4),sizes_orig(6)),[1 2 5 3 4 6]);
    mask = ipermute(reshape(mask,sizes_orig(1),[],sizes_orig(5),sizes_orig(3),sizes_orig(4),sizes_orig(6)),[1 2 5 3 4 6]);
    Nseg = size(mask,2);
    sizes = [1, 1, 1, 1, 1, 1];
    sizes(1:ndims(navData_tensor)) = size(navData_tensor);
end

try
if flagAutoLambdaTensor && ~strcmp(ScanType,'Cine') && ~isempty(navData_bloch)
    if lowmem == 0
        [lr,sms] = autolambdaTensor(navData_bloch(:,uni,:,:,:),ts_proj,dataArray.msdev);
    else
        [lr,sms] = autolambdaTensor(navData_bloch,ts_proj,dataArray.msdev);
    end
end
catch
end

t = toc;
strOutput = '\nRegularization parameters:';
strOutput = [strOutput '\n - lr    = ' num2str(lr)];
strOutput = [strOutput '\n - sms  = [' num2str(sms) ']'];
strOutput = [strOutput '\n Calculation time = ' num2str(t) ' sec\n'];
fprintf(strOutput);

reconOptions.tensor.lr = lr;
reconOptions.tensor.sms = sms;
reconOptions.tensor.smorders = smorders;
reconOptions.tensor.circs = circs;
reconOptions.NDCEAve = NDCEAve;
reconOptions.flagDataDriven = flagDataDriven;
reconOptions.flagUseFirstEcho = flagUseFirstEcho;

dataArray.wallClock = wallClock;
% dataArray.wallClock_full = wallClock_full;

tempStructTensor.cL        = cL;
tempStructTensor.NDCEAve   = NDCEAve;
tempStructTensor.Nseg      = Nseg;
tempStructTensor.Nseg_orig = Nseg_orig;
tempStructTensor.Segidx      = Segidx;
tempStructTensor.Segidx_orig = Segidx_orig;
tempStructTensor.uni         = uni;
tempStructTensor.isWallClockSet = isWallClockSet;
tempStructTensor.doBloch        = doBloch;
tempStructTensor.doSpline       = sum(abs(sms)) ~= 0;
tempStructTensor.isUndersampled = ~exp(sum(log(double(mask(:)))));  %~prod(mask(:));
tempStructTensor.navData_tensor = navData_tensor;
tempStructTensor.navData_bloch  = navData_bloch;
tempStructTensor.mask           = mask;
tempStructTensor.morozov = morozov;
tempStructTensor.sizes      = sizes;
tempStructTensor.sizes_orig = sizes_orig;
tempStructTensor.curvePhi      = curvePhi;
tempStructTensor.curvePhi_orig = curvePhi_orig;
tempStructTensor.ts_proj  = ts_proj;
if lowmem == 2
    tempStructTensor.navData_bloch_cL = navData_bloch_cL;
    tempStructTensor.sizes_lm         = sizes_lm;
end



%% License check
function isValid = checkLicense(reconOptions)

% get machine ID
[~,strAdd] = genID;

md = java.security.MessageDigest.getInstance('MD5');    
try
    % load license file
    licenseID = loadLicenseFile(reconOptions);
    [hashes,dateNr] = getHash(licenseID);
    
    % check hash
    for n = 1:length(strAdd)
        ID   = double(strAdd{n})*sum(dateNr);
        hash = dec2hex(uint8(double(md.digest(ID))+128));
        hash = hash(:).';
        isValid = contains(hashes,hash);
        if isValid; break; end
    end
    if ~isValid
        fprintf(2,'License check error: invalid hash.\n');
        disp(['Current date: ' date]);
        disp('Current system info:');
        for n = 1:length(strAdd)
            disp(['    ' strAdd{n}]);
        end
    end
catch errormsg
    fprintf(2,'License check error: %s\n', errormsg.message);
    isValid = false;
end