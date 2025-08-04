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
    
    if selectROI
        h = figure;imagesc(abs(fbpComposite(:,:,floor(Nz/2)+1,1)));axis equal tight;colormap('gray');title('Draw Resp motion ROI')
        roiResp = imellipse;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);

        Nxshift = roiPosition(1);
        Nyshift = roiPosition(2);
        ROINx   = roiPosition(3);
        ROINy   = roiPosition(4);

        roi_weighting = zeros(Ny,Nx,Nz);
        roi_weighting(Nyshift+(1:ROINy),Nxshift+(1:ROINx),Nzshift+(1:ceil(Nzorig))) = repmat(hanning(ROINy)*hanning(ROINx).',[1 1 Nzorig]); %ROI is center region of center slices
        roi_weighting = roi_weighting(1:Ny,1:Nx,1:Nz);
        windowfun = roi_weighting(:,:,floor(Nz/2)+1);
        roi_weighting = fftshift(roi_weighting,3);
    end

    Wti = reshape(U_init,[],L_init);
    Wti = bsxfun(@times,Wti,roi_weighting(:));
    Wti = Wti'*Wti;
    Wti = sqrtm(Wti); %actually inverse of Wt
%     U = vec(reshape(U_init,[],L_init)/Wti);

    Phi_rt       = Wti*Phi_rt_init;
    Phi_rt_full  = Wti*Phi_rt_full_init;
    Phi_rt_small = double(Wti*Phi_rt_small_init);

    temporalBasis.Phi_rt       = Phi_rt;
    temporalBasis.Phi_rt_full  = Phi_rt_full;
    temporalBasis.Phi_rt_small = Phi_rt_small;

    %% Respiratory window

    if sum(abs(diff(diff(navTiming)))) >= 1e-6 || usePMUdata% Changed condition for PMU data usage - XZ 11/19/2024
    % if sum(abs(diff(diff(navTiming)))) < 1e-6
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
        
    else
        dt  = navTiming(2) - navTiming(1);    
        fs  = 1/dt;
        fsm = 1/ACQTiming(linesPerShot+1);
        df  = fs/size(Phi_rt_small,2);
    end
    
    winlp = 2*floor((BRhigh/60)/df);
    windowlp = zeros(1,size(Phi_rt_small,2));
    windowlp(1:winlp) = hamming(winlp,'periodic');
    hwindow = circshift(windowlp,-winlp/2);

    if strcmp(ScanType, 'Cine')
        winhp = 2*floor((BRlow/60)/df); 
        windowhp = zeros(1,size(Phi_rt_small,2));
        windowhp(1:winhp) = 1;
        hwindow = hwindow.*(1-circshift(windowhp,-winhp/2));
    end

    winwidth = 0.1/dt; %100 ms window width
    winwidth = ceil((winwidth-1)/2)*2 + 1; %make odd

    %% Select candidate initial guesses
    if usePMUdata
        its = 1;
        outeritLoop = 1;
        Phi_rt_motion = ifft(fft(realify(physioNavResp)).*hwindow/sqrt(numel(physioNavResp)));
%         [B,TFrm] = rmoutliers(abs(Phi_rt_motion),"SamplePoints",1:numel(Phi_rt_motion));
%         physioNavResp(TFrm) = physioNavResp(TFrm)/max(physioNavResp(TFrm))*mean(physioNavResp(~TFrm));
%         Phi_rt_motion = ifft(fft(realify(physioNavResp)).*hwindow/sqrt(numel(physioNavResp)));
        resp_pct_idx = 1;
    else
        its = 35;
        outeritLoop = 1:min(3, L_init);
        Phi_rt_motion = Phi_rt_small - Phi_rt_small*pinv(curvePhi(Segidx,:).')*curvePhi(Segidx,:).'; %orthogonalize to Bloch subspace
        resp_pct = zeros(L_init,1);
        Phi_rt_motion2 = zeros(size(Phi_rt_motion));
        for outerit = 2:L_init
            Phi_rt_motion2(outerit,:) = ifft(fft(realify(Phi_rt_motion(outerit,:))).*hwindow/sqrt(size(Phi_rt_small,2)));
            resp_pct(outerit) = norm(fft(realify(Phi_rt_motion(outerit,:))).*hwindow/sqrt(size(Phi_rt_small,2)))/norm(Phi_rt_small(1,:));
        end
    
        [~,resp_pct_idx] = sort(resp_pct,'descend');
    end
    % plotPhi(Phi_rt_motion2,dt,1:8,'Phi_rt_motion2');    
%     if flagCommandLine
%         figure,plot(resp_pct,'.-'),drawnow;
%     end

    %% run binning

    Phiresp = zeros(rbins,L_init,cL);
    bestresnorm = inf;
    for outerit = outeritLoop %how many guesses?
        fprintf(' -- Respiratory binning iteration = %d: Phi_motion component %d, ', outerit, resp_pct_idx(outerit));

%         if flagUseEMD
%             Z = realify(Phi_rt_small(resp_pct_idx(outerit),:)).';
%             IMF = emd(Z);
%             temp0 = 0;
%             for n = 1:size(IMF,2)
%                 temp = norm(fft(IMF(:,n)).*hwindow);
%                 if temp > temp0
%                     temp0 = temp;
%                     Z = IMF(:,n);
%                 end
%             end
%         else
            Z = realify(Phi_rt_motion(resp_pct_idx(outerit),:));
            Z = real(ifft(fft(Z).*hwindow));
%         end
        
        Z = (Z-min(Z))/range(Z)*(rbins-1)+1;
        [Zn,~] = hist(Z*10,1:round(max(Z)*10));
        Zn = cumsum(Zn)/sum(Zn);
        Z = Zn(round(Z*10));
        Z = ceil(Z*rbins);
        Z(Z<1) = 1;
        Z(Z>rbins) = rbins;
        Z = Z(:).';

    %     resmat = zeros(rbins,numel(Z));
        for it = 1:its
            for j = 1:rbins
                Phiresp(j,:,:) = Phi_rt_small(:,Z==j) * pinv(curvePhi(Segidx(Z==j),:).');
            end

    %          parfor j=1:size(Phi_rt_small,2)
    %              resmat(:,j) = sum(abs(bsxfun(@minus,Phi_rt_small(:,j),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
    %          end

            temp   = reshape(reshape(Phiresp,[],cL)*curvePhi(:,1:cL).',rbins, L_init, []);
            resmat = squeeze(sum(abs(bsxfun(@minus, reshape(Phi_rt_small,1,L_init,[]), temp(:, :, Segidx))).^2,2));

            resmat = bsxfun(@minus,resmat,mean(resmat));
            resmat_filt = sgolayfilt(resmat.',0,winwidth).';
            [~,Z] = min(resmat_filt);
            Z = Z(:).';

            Z = real(ifft(fft(Z).*hwindow));
            Z = Z - min(Z) + 1;

            [Zn,~] = hist(Z*10,1:round(max(Z)*10));
            Zn = cumsum(Zn)/sum(Zn);
            Z  = Zn(round(Z*10));
            Z  = ceil(Z*rbins);
            Z(Z<1) = 1;
            Z(Z>rbins) = rbins;
            Z  = Z(:).';

        end

        for j = 1:rbins
            Phiresp(j,:,:) = Phi_rt_small(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
        end

%         res = zeros(size(Z));
%         parfor j = 1:size(Phi_rt_small,2)
%             res(j) = norm(Phi_rt_small(:,j) - squeeze(Phiresp(Z(j),:,:))*curvePhi(Segidx(j),:).');
%         end
        res = sqrt(sum(abs(Phi_rt_small.' - sum(Phiresp(Z,:,:).*permute(curvePhi(Segidx,:),[1 3 2]),3)).^2,2));
        resnorm = norm(res);
        fprintf('resnorm = %f\n', resnorm);

        if resnorm < bestresnorm
            Ridx = Z;
            bestres = res;
            bestresnorm = resnorm;
            bestPhi = resp_pct_idx(outerit);
        end
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

    if sum(abs(diff(diff(navTiming)))) >= 1e-6 || usePMUdata
        dataArray.Ridx_full   = row(Ridx);
        Ridx = interp1(newTiming,Ridx,navTiming,'nearest','extrap');
        Segidx = mod(navIndices-1,size(curvePhi,1)) + 1;
    end
else
    Ridx = ones(size(Phi_rt_small_init,2),1);
end

dataArray.Ridx   = (Ridx)';
dataArray.Segidx = (Segidx)';


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
    %temp(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
    %                *mean(Phi_rt_small(:,Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[])); % XZ 11/19/2024
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