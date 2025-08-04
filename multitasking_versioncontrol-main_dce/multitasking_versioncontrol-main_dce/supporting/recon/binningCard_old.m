function [dataArray,temporalBasis] = binningCard_old(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI)
% v2.0

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Invalid hash: Multitasking "license check" failed');
%     waitfor(dlg);
%     return;
% end

if nargin < 6
    selectROI = false;
end

% initialize heart rate filter range (cycle/min)
% will be overwritten if set in reconOptions
% HRlow  = 20;
% HRhigh = 70;

% load variables
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L_init = size(Phi_rt_small_init,1);
ccL = min(10,size(curvePhi,2)); 
Segidx = mod(navIndices-1,linesPerShot*moduleLength)+1;

vec = @(x) x(:);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end

if cbins > 1

    % Set looping params
    % If binning quality is not good, try:
    % 1) change max number of outerit loop
    % 2) change number of iterations its
    % 3) change alpha
    its = 35;
    alpha = cbins*0.5;  
    outeritLoop = 1:min(2, L_init);


    %% weighted images

    if selectROI
        h = figure;imagesc(abs(fbpComposite(:,:,floor(Nz/2)+1,1)));axis equal tight;colormap('gray');title('Draw cardiac motion ROI')
        roiResp = imellipse;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);

        Nzshift = floor(Nz/2)- floor(Nzorig/2);
        Nxshift = roiPosition(1);
        Nyshift = roiPosition(2);
        ROINx = roiPosition(3);
        ROINy = roiPosition(4);

        windowfun = zeros(Ny,Nx);
        windowfun(Nyshift+(1:ROINy),Nxshift+(1:ROINx)) = hanning(ROINy)*hanning(ROINx).';  
        roi_weighting = zeros(Ny,Nx,Nz);
        roi_weighting(:,:,Nzshift+(1:Nzorig)) = repmat(windowfun,[1 1 Nzorig]); 
    else
%         tempLength = min(size(Phi_rt_small_init,2),4000);
%         fs = 1/(lEchoSpacing*SGBlock);
%         df = fs/tempLength;
%         winlp = 2*floor((HRhigh/60)/df); %highest HR: 130 bpm
%         windowlp = zeros(tempLength,1);
%         windowlp(1:winlp) = 1;
%         windowlp = circshift(windowlp,[-winlp/2, 0]);
% 
%         winhp = 2*floor((HRlow/60)/df); %lowest HR: 50 bpm
%         windowhp = zeros(tempLength,1);
%         windowhp(1:winhp) = 1;
%         hwindow = windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
%     
%         imgs = reshape(Utemp(:,:,floor(Nz/2)+1,:),[],L_init)*Phi_rt_small_init(:,1:tempLength);
%         imgs = abs(reshape(imgs,Ny,Nx,[]));
%         imgs_fft = fft(ifftshift(imgs,3),[],3);
%         imgs_h = sum(abs(ifft(bsxfun(@times,imgs_fft,reshape(hwindow,1,1,[])),[],3)).^2,3);
%         [x,y] = meshgrid(1:size(imgs_h, 2), 1:size(imgs_h, 1));
%         weightedx = x .* imgs_h;
%         weightedy = y .* imgs_h;
%         xcenter = sum(weightedx(:)) / sum(imgs_h(:));
%         ycenter = sum(weightedy(:)) / sum(imgs_h(:));

        % Define window
        Ncard = 2*floor(3*Norig/8);
        Nyshift = floor(Ny/2) - Ncard/2;
        Nxshift = floor(Nx/2) - Ncard/2;
        Nzshift = floor(Nz/2)- floor(Nzorig/2);
        windowfun = zeros(Ny,Nx);
        windowfun(Nyshift+(1:Ncard),Nxshift+(1:Ncard)) = hanning(Ncard)*hanning(Ncard).';   %ROI is center region of center slices
        %windowfun = circshift(windowfun,floor([ycenter-floor(Ny/2)-1,xcenter-floor(Nx/2)-1]));
        roi_weighting = zeros(Ny,Nx,Nz);
        roi_weighting(:,:,Nzshift+(1:Nzorig)) = repmat(windowfun,[1 1 Nzorig]);
    end

    if MBfactor == 1
        roi_weighting = ifftshift(roi_weighting,3);
    end
    if isCartesian
        roi_weighting = ifftshift(roi_weighting,1);
    end

    Wti = reshape(U_init,[],L_init);
    Wti = bsxfun(@times,Wti,roi_weighting(:));
    Wti = Wti'*Wti;
    Wti = sqrtm(Wti); %actually inverse of Wt
    U = vec(reshape(U_init,[],L_init)/Wti);

    Phi_rt       = Wti*Phi_rt_init;
    Phi_rt_full  = Wti*Phi_rt_full_init;
    Phi_rt_small = Wti*Phi_rt_small_init;

    temporalBasis.Phi_rt = Phi_rt;
    temporalBasis.Phi_rt_full  = Phi_rt_full;
    temporalBasis.Phi_rt_small = Phi_rt_small;

    %% Cardiac window
        
    fs = 1/(lEchoSpacing*SGBlock);
    df = fs/size(Phi_rt_small,2);
        
    winlp = 2*floor((HRhigh/60)/df); %highest HR: 130 bpm
    windowlp = zeros(size(Phi_rt_small,2),1);
    windowlp(1:winlp) = 1;
    windowlp = circshift(windowlp,[-winlp/2, 0]);

    winhp = 2*floor((HRlow/60)/df); %lowest HR: 40 bpm
    windowhp = zeros(size(Phi_rt_small,2),1);
    windowhp(1:winhp) = 1;
    hwindow = windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
    
    winwidth = 0.035/(lEchoSpacing*SGBlock); %35 ms window width
    winwidth = ceil((winwidth-1)/2)*2 + 1; %make odd

    cbins_half = ceil(cbins/2);
    cbins_full = cbins_half*2;
    
    %% Select candidate initial guesses

    Phi_rt_motion = Phi_rt_small - Phi_rt_small*pinv(curvePhi(Segidx,:).')*curvePhi(Segidx,:).';  %orthogonalize to Bloch subspace
%     Phi_rt_motion = Phi_rt_small_init;
%     Phi_rt_motion(1,:) = repmat(Phi_rt_small_init(1,340),1,size(Phi_rt_small_init,2));
    card_pct = zeros(L_init,1);
    for outerit = 2:L_init
        card_pct(outerit) = norm(fft(realify(Phi_rt_motion(outerit,:))).'.*hwindow/sqrt(size(Phi_rt_small,2)))/norm(Phi_rt_small(1,:));
    end
    [~,card_pct_idx] = sort(card_pct,'descend');

    % if flagCommandLine
    %     figure,plot(card_pct,'.-'),drawnow
    % end

    %%
    bestresnorm = inf;

    for outerit = outeritLoop
        fprintf(' -- Cardiac binning iteration %d: Phi_motion component %d, ', outerit, card_pct_idx(outerit));

        Phicard = zeros(rbins,cbins_half,L_init,ccL);
        [Z,Phicard,res,resnorm] = runBinningCardLoop(alpha,its,Ridx,Phi_rt_motion(card_pct_idx(outerit),:),Phi_rt_small,curvePhi,Phicard,Segidx,hwindow,winwidth);
        
        if resnorm < bestresnorm
            Phicardbest = Phicard;
            Hidx = Z;
            bestres = res;
            bestresnorm = resnorm;
            bestPhi = card_pct_idx(outerit);
        end
    end
    fprintf('Use Phi #%d\n', bestPhi);
else
    cbins_full = 1;
    Hidx    = ones(size(Ridx));
    bestres = ones(size(Ridx))/numel(Ridx);
end

%%
% cbins   = cbins_full;
% Phicard = Phicardbest;
% bestresnorm
% try
%     temp = reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL)*curvePhi(120,1:ccL).',[],L_init).',Nydisp,Nxdisp,[]);
%     implay(2*abs(temp)/max(abs(temp(:))))
% catch
%     temp = pinv(Phi_rt_small.')*nav_data(:,:);
%     temp = temp.'*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL)*curvePhi(end,1:ccL).',[],L_init).';
%     imagesc(abs(temp))
% end
  
% figure,hist(diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*SGblock*1000,10)
% xlabel('RR interval (ms)')

bestresCard = accumarray(Hidx(Ridx==1),bestres(Ridx==1),[],@(x)norm(x)/sqrt(numel(x)));

%% wall cloack
switch ScanType
    case 'Cine'
        wall_clock = 1 + cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1;
    case {'T2prep','T2IR'}
        wall_clock = 1 + cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock = ceil(wall_clock/HRinterv);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1;
    otherwise
        wall_clock = 1+cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1; %the last heartbeat may be incomplete, so group it with the next-to-last
        wall_clock = max(wall_clock-1,1); %the first heartbeat may be incomplete, so group it with the second
end

dataArray.Hidx = Hidx;
dataArray.bestresCard = bestresCard;
dataArray.bestres     = bestres;
dataArray.wallClock   = wall_clock;

temp = fft(Hidx);
temp(1) = 0;
[~,idx] = max(abs(temp(1:floor(numel(Hidx)/2))));
dataArray.meanHBPM = round(60*idx/(lEchoSpacing*SGBlock*numel(Hidx)));

HRinterv_temp=1;
HRintv=diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*2000;
HRintv=movmean(HRintv,4);
figure;plot(HRintv(1:HRinterv_temp:end)/1000);
ylabel('RR interval (s)');
dataArray.RRinterval = HRintv(1:HRinterv_temp:end)/1000;

% dataArray.diastoleIdx = 1;
% dataArray.Hidx = mod(Hidx - dataArray.diastoleIdx, max(Hidx)) + 1;


%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

if cbins > 1
%     Utemp = Utemp.*windowfun;

    temp = zeros(Nydisp,Nxdisp,numel(dispSlice),cbins_full,rbins);
    for j = 1:rbins
        for k = 1:cbins_full
            temp(:,:,:,k,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                              *mean(Phi_rt_small_init(:,Hidx==k & Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[]));
        end
    end
    cw = prctile(temp(:),99.5);
    dataArray.binsCardMean = abs(temp)/cw;

    dataArray.binsCard = cell(cbins,1);
    for k = 1:cbins_full
            temp = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                              *Phi_rt_small_init(:,Hidx==k),Nydisp,Nxdisp,numel(dispSlice),[]));
            dataArray.binsCard{k} = abs(temp)/prctile(temp(:),99.5);
    end

    dataArray.diastoleIdx = testLVvolume(squeeze(abs(temp)));
    dataArray = shiftHidx(params,reconOptions,dataArray,temporalBasis,spatialCoeff);
else
    dataArray.diastoleIdx = 1;
end

% if flagCommandLine && cbins > 1
%     h = findall(groot,'Type','figure','Name','Binning Results');
%     if isempty(h)
%         figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
%     else
%         figure(h);
%     end
%     subplot(2,2,3),plot(Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
%     subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
%     subplot(2,2,1),plot(Hidx,'.-');axis([-inf inf 0 cbins_full+1]);title(sprintf('Hidx, mean BPM %d',dataArray.meanHBPM));
%     subplot(2,2,2);
%     for rphase = 1:rbins
%         plot(Segidx((Ridx==rphase)), Hidx((Ridx==rphase)),'.');hold on;
%     end
%     hold off; axis([0 linesPerShot*moduleLength 0 cbins_full+1]);title('Hidx / Segment Index');
%     
%     implayZoom(squeeze(abs(dataArray.binsCard(:,:,1,:,1))),cbins_full);
% end


function  [Z,Phicard,res,resnorm] = runBinningCardLoop(alpha,its,Ridx,Phi_rt_motion,Phi_rt_small,curvePhi,Phicard,Segidx,hwindow,winwidth)

    [rbins,cbins_half,L_init,ccL] = size(Phicard);

%         if flagUseEMD  % try EMD
%             Z = realify(Phi_rt_motion).';   
%             IMF = emd(Z/sqrt(numel(Z)));
%             temp0 = 0;
%             for n = 1:size(IMF,2)
%                 temp = norm(fft(IMF(:,n)).*hwindow);
%                 if temp > temp0
%                     temp0 = temp;
%                     Z = IMF(:,n);
%                 end
%             end
%         else
         Z = realify(Phi_rt_motion).';   
         Z = fft(Z/sqrt(numel(Z))).*hwindow;
         Z = real(ifft(Z));
%         end

    Z = Z./abs(hilbert(Z));
    Z = ceil((Z+1)*cbins_half/2);
    Z(Z<1) = 1;
    Z(Z>cbins_half) = cbins_half;
    %resmat = zeros(cbins_half,numel(Z));
    for it = 1:its
        for j = 1:rbins
            for k = 1:cbins_half
                Phicard(j,k,:,:) = Phi_rt_small(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
            end
        end
        %switch systole/diastole?
        nucnorm = zeros(2^(rbins-1),1);
        temp = rbins;
        for j = 1:2^(temp-1) %never switch first one
            reverse_these = (dec2bin(j-1,temp)=='1');
            Phicardtemp = Phicard;
            Phicardtemp(reverse_these,:,:,:) = flip(Phicardtemp(reverse_these,:,:,:),2);
            nucnorm(j)  = sum(svde(Phicardtemp(:,:)));
            Phicardtemp = permute(Phicardtemp,[2 1 3 4]);
            nucnorm(j)  = nucnorm(j) + sum(svde(Phicardtemp(:,:)));
        end
        [~,j] = min(nucnorm);
        reverse_these = (dec2bin(j-1,rbins)=='1');
        Phicard(reverse_these,:,:,:) = flip(Phicard(reverse_these,:,:,:),2);

%         temp = Ridx;
%         parfor j = 1:size(Phi_rt_small,2)
%             resmat(:,j) = sum(abs(bsxfun(@minus,Phi_rt_small(:,j),reshape(reshape(Phicard(temp(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',cbins_half,[]).')).^2);
%         end
        
        resmat = squeeze(sum(abs(reshape(Phi_rt_small.',[],1,L_init) - sum(Phicard(Ridx,:,:,:).*reshape(curvePhi(Segidx,1:ccL),[],1,1,ccL),4)).^2,3)).';
        
        resmat = bsxfun(@minus,resmat,mean(resmat));
        resmat_filt = sgolayfilt(double(resmat).', 0, winwidth).';
        [~,Z] = min(resmat_filt);
        Z = Z(:);

        Z = fft(Z/sqrt(numel(Z))).*hwindow;
        if (alpha > 0) && (it > 10) && (max(abs(Z))>alpha)
            Z = sign(Z).*max(abs(Z)-alpha,0);
        end
        Z = real(ifft(Z));
        if it < its
            Z = Z./abs(hilbert(Z));
            Z = ceil((Z+1)*cbins_half/2);
            Z(Z<1) = 1;
            Z(Z>cbins_half) = cbins_half;
        else
            Z = angle(hilbert(Z));
            Z = ceil((Z/pi+1)*cbins_half);
            Z(Z==0) = 1;
        end
    end

    Phicard = zeros(rbins,cbins_half*2,L_init,ccL);
    for j = 1:rbins
        for k = 1:cbins_half*2
            Phicard(j,k,:,:) = Phi_rt_small(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
        end
    end

%     res = zeros(size(Z));
%     temp = Ridx;
%     parfor j = 1:size(Phi_rt_small,2)
%         res(j) = norm(Phi_rt_small(:,j)-squeeze(Phicard(temp(j),Z(j),:,:))*curvePhi(Segidx(j),1:ccL).');
%     end
    ind = sub2ind([rbins,cbins_half*2],Ridx,Z);
    Phicard = reshape(Phicard,[],L_init,ccL);
    res = sqrt(sum(abs(Phi_rt_small.' - sum(Phicard(ind,:,:).*reshape(curvePhi(Segidx,1:ccL),[],1,ccL),3)).^2,2));
    
    resnorm = norm(res);
    fprintf('resnorm = %f\n', resnorm);


%% License check
function isValid = checkLicense(reconOptions)

% get machine ID
if ispc
    [~,info] = system('getmac');
elseif ismac
    [~,info] = system('ifconfig en0 | grep ether');
elseif isunix
    [~,info] = system('ip addr | grep ether');
else
    error('OS not recognized.');
end

md = java.security.MessageDigest.getInstance('MD5');    
try
    % load license file
    [licenseID,flag] = loadLicenseFile(reconOptions);
    licenseID = licenseID/sum(double('multitasking'));
    if flag == 1        % offline license check
        datestr = licenseID(end);
        licenseID(end) = [];
        licenseDays = datestr - datenum(date);
        if licenseDays > 0 
            ID = double(info) * datestr;
            hash = dec2hex(uint8(double(md.digest(ID))+128));
            hash = hash(:).';
            isValid = contains(char(licenseID),hash);
        else
            error('License expired.');
        end
    else                % online license check
        % get Hash
        ID   = double(info)*sum(double(date));
        hash = dec2hex(uint8(double(md.digest(ID))+128));
        hash = hash(:).';

        %hashes  = webread(['https://agchristodoulou.github.io/MTcheck/' htmlID '.html']);
        URL = char(licenseID/sum(double('multitasking')));
        hashes  = webread(URL);
        isValid = contains(hashes,hash);
    end
catch
    isValid = false;
end
