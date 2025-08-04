function [dataArray,temporalBasis] = binningCard_SA(params,reconOptions,dataArray,temporalBasis,spatialCoeff)
% v2.0

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

if nargin < 6
    selectROI = false;
end

% initialize heart rate filter range (cycle/min)
% will be overwritten if set in reconOptions
HRlow  = 40;
HRhigh = 130;

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
    dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

    h = figure;imagesc(abs(dispim(fbpComposite(:,:,:,1))));axis equal tight;colormap('gray');title('Draw cardiac motion ROI')
    roiResp = imrect;
    roiPosition = roiResp.getPosition();
    roiPosition = round(roiPosition);
    close(h);

    Nzshift = floor(Nz/2)- floor(Nzorig/2);
    Nxshift = floor(Nx/2)-floor(Nxdisp/2)+roiPosition(1);
    Nyshift = floor(Ny/2)-floor(Nydisp/2)+roiPosition(2);
    ROINx = roiPosition(3);
    ROINy = roiPosition(4);

    windowfun = hanning(ROINy)*hanning(ROINx).';  
    windowfun(windowfun<0.02) = 0;
    windowfun(windowfun>0) = 1; 
%     roi_weighting = zeros(Ny,Nx);
%     roi_weighting(Nyshift+(1:ROINy),Nxshift+(1:ROINx)) = windowfun; 
%     se = strel('disk',3);
%     roi_weighting = imdilate(roi_weighting,se);
%     roi_weighting = repmat(roi_weighting,[1 1 Nz]); 
%     roi_weighting(:,:,1:Nzshift) = 0;
%     roi_weighting(:,:,Nzshift+1+Nzorig:end) = 0;
    
    % image of the first 20 seconds
    Utemp2 = Utemp(Nyshift+(1:ROINy),Nxshift+(1:ROINx),:);
    Phitemp = Phi_rt_small_init(:,1:floor(20/lEchoSpacing/SGBlock));
    recon = reshape(reshape(Utemp2,[],L_init)*Phitemp,size(Utemp2,1),size(Utemp2,2),[]);
    recon = recon ./max(abs(recon(:)));
    temp  = mean(abs(recon.*windowfun),3);
   
    h = figure;imagesc(mean(abs(recon),3));axis equal tight;colormap('gray');title('Select an ROI inside LV blood pool');
    roiResp = imellipse;
    roiPosition = roiResp.getPosition();
    roiPosition = round(roiPosition);
    close(h);

    Nxshiftb = roiPosition(1);
    Nyshiftb = roiPosition(2);
    ROINxb = roiPosition(3);
    ROINyb = roiPosition(4);

    mask = zeros(size(recon,1),size(recon,2));
    mask(Nyshiftb+(1:ROINyb),Nxshiftb+(1:ROINxb)) = hanning(ROINyb)*hanning(ROINxb).'; 
    mask(mask>0.05) = 1;
    mask(mask<=0.05) = 0;
    %figure;imshow(mask);
 
    maskh = grayconnected(temp,floor(Nyshiftb+ROINyb/2),floor(Nxshiftb+ROINxb/2));
    
    mask1 = mask(:);
    recon1 = reshape(recon,numel(mask),[]);
    recon1 = realify(recon1,'rows');
    curve_blood = (mean(recon1(mask1>0,:))).';
    curve_blood(:,2) = 1;
    %figure;plot(curve_blood);

    res = zeros(size(recon1,1),size(recon1,2));
    beta = zeros(size(recon1,1),2);
    resnorm = zeros(size(recon1,1),1);
    for n = 1:size(recon1,1)
        b = curve_blood\recon1(n,:).';
        res(n,:)  = recon1(n,:).' - curve_blood*b;
        beta(n,:) = b.';
        resnorm(n) = norm(res(n,:))/norm(recon1(n,:));
        res(n,:) = res(n,:)/resnorm(n,1);
    end

    temp1 = reshape(resnorm,size(mask,1),size(mask,2));
    temp1(temp1>mean(temp1(:)))  = 0;
    temp1(temp1>mean(temp1(maskh>0))) = 0;
    %figure;imagesc(temp1);

    se = strel('disk',3);
    temp2 = imclose(temp1,se);
    temp2 = imopen(temp2,se);
    %figure;imagesc(temp2);
    %figure;imagesc(imfill(temp2));
    mask1 = temp2;
    mask1(mask1>0) = 1;
    %figure;imagesc(mask1);

    mask2 = imfill(temp2);
    mask2(mask2>0) = 1;
    se = strel('disk',3);
    mask2 = imclose(mask2,se);
    mask2 = imdilate(mask2,se);
    %figure;imagesc(mask2);
    
    roi_weighting1 = zeros(Ny,Nx,Nz);
    roi_weighting1(Nyshift+(1:ROINy),Nxshift+(1:ROINx),:) = repmat(mask2,[1 1 Nz]);  
    
    if MBfactor == 1
        roi_weighting1 = ifftshift(roi_weighting1,3);
    end
    if isCartesian
        roi_weighting1 = ifftshift(roi_weighting1,1);
    end

    Wti = reshape(U_init,[],L_init);
    Wti = bsxfun(@times,Wti,roi_weighting1(:));
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
        
    dt  = lEchoSpacing*SGBlock;    
    fs  = 1/dt;
    fsm = 1/(lEchoSpacing*linesPerShot);
    df  = fs/size(Phi_rt_small,2);
        
    winlp = 2*floor((HRhigh/60)/df); 
    windowlp = zeros(1,size(Phi_rt_small,2));
    windowlp(1:winlp) = 1;
    windowlp = circshift(windowlp,-winlp/2);

    winhp = 2*floor((HRlow/60)/df); 
    windowhp = zeros(1,size(Phi_rt_small,2));
    windowhp(1:winhp) = 1;
    hwindow = windowlp.*(1-circshift(windowhp,-winhp/2));
    
    winwidth = 0.035/(lEchoSpacing*SGBlock); %35 ms window width
    winwidth = ceil((winwidth-1)/2)*2 + 1; % make odd

    cbins_half = ceil(cbins/2);
    cbins_full = cbins_half*2;
    
    %% Select candidate initial guesses

    n = size(Phi_rt_small,2);
    f(1:ceil(n/2)) = (0:ceil(n/2)-1)*(fs/n);            % frequency range
    f(ceil(n/2)+1:n) = (-floor(n/2):1:-1)*(fs/n);       % frequency range
 
    Phi_rt_motion = Phi_rt_small;
    %plotPhi(Phi_rt_small,lEchoSpacing*SGBlock,1:q,'Phi_rt_small');
    Phi_rt_motion = Phi_rt_motion - Phi_rt_motion*pinv(curvePhi(Segidx,:).')*curvePhi(Segidx,:).';  %orthogonalize to Bloch subspace
    %plotPhi(Phi_rt_motion,lEchoSpacing*SGBlock,1:q,'Phi_rt_motion1');
    
%     % notch filter foro prep pulse frequency and harmonics
%     Phi_rt_motion = double(Phi_rt_motion);
%     for fn = 1:floor(max(f)/fsm)
%         d(fn) = designfilt('bandstopiir','FilterOrder',2, ...
%             'HalfPowerFrequency1',fsm*fn-df/0.5,'HalfPowerFrequency2',fsm*fn+df/0.5, ...
%             'DesignMethod','butter','SampleRate',fs);
%         for i = 2:L_init
%             Phi_rt_motion(i,:)= filtfilt(d(fn),Phi_rt_motion(i,:));
%         end
%     end
    
    card_pct = zeros(L_init,1);
    Phi_rt_motion2 = zeros(size(Phi_rt_motion));
    for outerit = 2:L_init
        Phi_rt_motion2(outerit,:) = ifft(fft(realify(Phi_rt_motion(outerit,:))).*hwindow/sqrt(size(Phi_rt_small,2)));
        card_pct(outerit) = norm(fft(realify(Phi_rt_motion(outerit,:))).*hwindow/sqrt(size(Phi_rt_small,2)))/norm(Phi_rt_small(1,:));
    end
    [~,card_pct_idx] = sort(card_pct,'descend');

    % plotPhi(Phi_rt_motion2,dt,1:8,'Phi_rt_motion2');
    
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

bestresCard = accumarray(vec(Hidx(Ridx==1)),vec(bestres(Ridx==1)),[],@(x)norm(x)/sqrt(numel(x)));

%% wall cloack
switch ScanType
    case 'Cine'
        wall_clock = 1 + cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1;
    case 'T2prep'
        wall_clock = 1 + cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock = ceil(wall_clock/5);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1;
    otherwise
        wall_clock = 1+cumsum(diff(Hidx)==(1-cbins));
        wall_clock(end+1) = wall_clock(end);
        wall_clock(wall_clock == wall_clock(end)) = wall_clock(end)-1; %the last heartbeat may be incomplete, so group it with the next-to-last
        wall_clock = max(wall_clock-1,1); %the first heartbeat may be incomplete, so group it with the second
end

dataArray.Hidx = Hidx(:).';
dataArray.bestresCard = bestresCard;
dataArray.bestres     = bestres(:).';
dataArray.wallClock   = wall_clock;
dataArray.meanHBPM    = round(numel(find(diff(Hidx)<0))*60/(lEchoSpacing*SGBlock*numel(Hidx)));

% dataArray.diastoleIdx = 1;
% dataArray.Hidx = mod(Hidx - dataArray.diastoleIdx, max(Hidx)) + 1;


%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

if cbins > 1
    %Utemp = Utemp.*roi_weighting;

    temp = zeros(Nydisp,Nxdisp,numel(dispSlice),cbins_full,rbins);
    for j = 1:rbins
        for k = 1:cbins_full
            temp(:,:,:,k,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                              *mean(Phi_rt_small_init(:,Hidx==k & Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[]));
        end
    end
    cw = prctile(temp(:),99.5);
    dataArray.binsCard = abs(temp)/cw;
  
    mask1 = zeros(Ny,Nx,Nz);
    mask1(Nyshift+(1:ROINy),Nxshift+(1:ROINx)) = mask;
    mask = dispim(mask1);
    dataArray.diastoleIdx = testLVvolume(squeeze(abs(temp(:,:,1,:,1))),mask);
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
         Z = realify(Phi_rt_motion);   
         Z = fft(Z/sqrt(numel(Z))).*hwindow(:).';
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
                Phicard(j,k,:,:) = Phi_rt_small(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)).',1:ccL).');
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
        Z = Z(:).';

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