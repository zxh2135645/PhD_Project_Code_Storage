function dataArray = binningResp(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI,usePMUdata)
%v2.0

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end

if nargin < 7
    usePMUdata = false;
end
  
if nargin < 6
    selectROI = false;
end

% initialize breathing rate filter range (cycle/min)
% will be overwritten if set in reconOptions
BRlow  = 3;
BRhigh = 50;

% load variables
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

if ~exist('physio','var') || sum(diff(physio.nav.RESP(:,:,1))) == 0
    usePMUdata = false;
end

% if strcmp(params.ScanType,'SR')  % respiratory dimension for perfusion
%     rbins = 1;
% end

if strcmp(ScanType,'CEST')
    curvePhi = curvePhi_binning;
end

L_init = size(Phi_rt_small_init,1);
cL = size(curvePhi,2);
Segidx = mod(navIndices-1,size(curvePhi,1)) + 1;

vec = @(x) x(:);
row = @(x) x(:).';

if rbins > 1

    %% Prepare weighted images
    Nzshift = floor(Nz/2) - floor(Nzorig/2);
    Nyshift = floor(Ny/2) - floor(Nydisp/2);
    Nxshift = floor(Nx/2) - floor(Nxdisp/2);
    roi_weighting = zeros(Ny,Nx,Nz);
    roi_weighting(Nyshift+(1:2*floor(Nydisp/2)),Nxshift+(1:2*floor(Nxdisp/2)),Nzshift+(1:Nzorig)) = repmat(hanning(2*floor(Nydisp/2))*hanning(2*floor(Nxdisp/2)).',[1 1 Nzorig]); %ROI is center region of center slices
    roi_weighting = fftshift(roi_weighting,3);
    windowfun = roi_weighting(:,:,Nzshift+1);

    Wti = reshape(U_init,[],L_init);
    Wti = bsxfun(@times,Wti,roi_weighting(:));
    Wti = Wti'*Wti;
    Wti = sqrtm(Wti); %actually inverse of Wt

    Phi_rt       = Wti*Phi_rt_init;
    Phi_rt_full  = Wti*Phi_rt_full_init;
    Phi_rt_small = double(Wti*Phi_rt_small_init);
    
    %% Respiratory window
    newTiming = 0:floor(lEchoSpacing*SGBlock/2*1e3)/1e3:max(navTiming);
    Phi_rt_small_orig = Phi_rt_small;
    physioNavResp = row(interp1(navTiming,physio.nav.RESP,newTiming,'pchip','extrap'));
    Phi_rt_small = interp1(navTiming,realify(Phi_rt_small.','cols'),newTiming,'pchip','extrap').';
    curvePhi = interp1(navTiming,realify(curvePhi(Segidx,:),'cols'),newTiming,'pchip','extrap');
    temp = [0 diff(mod(navIndices-1,linesPerShot)+1)];
    idx = find(temp<0);
    for n = 1:numel(idx)
        zeroidx = intersect(find(newTiming>navTiming(idx(n)-1)),find(newTiming<navTiming(idx(n))));
        physioNavResp(zeroidx) = 0;
        Phi_rt_small(:,zeroidx) = 0;
        curvePhi(zeroidx,:) = 0;
    end
    Segidx = 1:size(Phi_rt_small,2);
    dt = floor(lEchoSpacing*SGBlock/2*1e3)/1e3;
    fs  = 1/dt;
    fsm = 1/ACQTiming(linesPerShot+1);
    df  = fs/size(Phi_rt_small,2);


    %% run binning
    Z = physioNavResp;  
    Z = (Z-min(Z))/range(Z)*(rbins-1)+1;
    [Zn,~] = hist(Z*10,1:round(max(Z)*10));
    Zn = cumsum(Zn)/sum(Zn);
    Z = Zn(round(Z*10));
    Z = ceil(Z*rbins);
    Z(Z<1) = 1;
    Z(Z>rbins) = rbins;
    Z = Z(:).';

    Ridx = Z;

    Phiresp = zeros(rbins,L_init,cL);
    bestresnorm = inf;
    resp_pct_idx = 1;
    outeritLoop = 1;
    outerit = outeritLoop;

    for j = 1:rbins
        Phiresp(j,:,:) = Phi_rt_small(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
    end

    res = sqrt(sum(abs(Phi_rt_small.' - sum(Phiresp(Z,:,:).*permute(curvePhi(Segidx,:),[1 3 2]),3)).^2,2));
    resnorm = norm(res);
    fprintf('resnorm = %f\n', resnorm);

    if resnorm < bestresnorm
        Ridx = Z;
        bestres = res;
        bestresnorm = resnorm;
        bestPhi = resp_pct_idx(outerit);
    end

    bestresResp = accumarray(Ridx(:),bestres(:),[],@(x)norm(x)/sqrt(numel(x)));

    fprintf('Use Phi #%d\n', bestPhi);

    if bestresResp(1) > bestresResp(end)
        Ridx = rbins - Ridx + 1;
        bestresResp = flip(bestresResp);
    end
    dataArray.bestresResp = bestresResp;

    temp = fft(Ridx);
    temp(1) = 0;
    [~,idx] = max(abs(temp));
    dataArray.meanRPeriod = 1/(df*(idx-1));

    dataArray.Ridx_full   = row(Ridx);
    Ridx = interp1(newTiming,Ridx,navTiming,'nearest','extrap');
    Segidx = mod(navIndices-1,size(curvePhi,1)) + 1;

    dataArray.Ridx   = row(Ridx);
    dataArray.Segidx = row(Segidx);


%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end
% if rbins > 1
%     Utemp = Utemp.*windowfun;
% end

temp = zeros(Nydisp,Nxdisp,numel(dispSlice),rbins);
for j = 1:rbins
    temp(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                    *mean(Phi_rt_small_init(:,Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[])); 
    % temp(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
    %                 *mean(Phi_rt_small(:,Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[])); % XZ 11/19/2024
end
cw = prctile(temp(:),99.9);
dataArray.binsRespMean = abs(temp)/cw;

% for 3D volume acquired in transverse orientation, also display resp bin images in coronal view
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);
if Nz > 1 && abs(zDir) == 3 && flagCommandLine && rbins > 1
    bins = zeros(Nydisp,Nxdisp,Nz,rbins);      
    for j = 1:rbins
        bins(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)*mean(Phi_rt_small_init(:,Ridx==j),2),Nydisp,Nxdisp,Nz,1));
        % bins(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)*mean(Phi_rt_small(:,Ridx==j),2),Nydisp,Nxdisp,Nz,1)); % XZ 11/19/2024
    end
    if abs(xDir) == 1
        bins_slice = permute(squeeze(bins(:,floor(Nxdisp/2),:,:)),[2 1 3]);  % select the slice here
        bins_slice = imresize(bins_slice,[floor(Nz*voxelSpacing(3)/voxelSpacing(1)) Ny]);
    else
        bins_slice = permute(squeeze(bins(floor(Nydisp/2),:,:,:)),[2 1 3]);  % select the slice here
        bins_slice = imresize(bins_slice,[floor(Nz*voxelSpacing(3)/voxelSpacing(2)) Nx]);
    end
    if xDir < 0
        bins_slice = flip(bins_slice,2);
    end
    if zDir > 0
        bins_slice = flip(bins_slice,1);
    end
    cw = prctile(bins_slice(:),99.9);
    bins_slice = bins_slice/cw;
    implayZoom(bins_slice,2);
    dataArray.bins_slice = bins_slice;
end

% dataArray.binsResp = cell(rbins,1);
% for j = 1:rbins
%     temp = abs(reshape(reshape(dispim(U_init),[],L_init)...
%                     *Phi_rt_small_init(:,Ridx==j),Nydisp,Nxdisp,numel(dispSlice),[]));
%     cw = prctile(temp(:),99);
%     dataArray.binsResp{j} = abs(temp)/cw;
% end

dt = lEchoSpacing*SGBlock;

if flagCommandLine && rbins > 1
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(dt:dt:dt*numel(Ridx),Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
    
    implayZoom(imageOrientLPS(dataArray.binsRespMean(:,:,1,:),params),2);
end
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