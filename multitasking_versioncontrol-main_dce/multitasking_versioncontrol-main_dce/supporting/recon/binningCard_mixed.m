function [dataArray,temporalBasis] = binningCard(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI)
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
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice,:);
disp2d = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:,:);

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
        
        % for testLVvolume, experimental
        h = figure;imagesc(abs(fbpComposite(:,:,floor(Nz/2)+1,1)));axis equal tight;colormap('gray');title('Select an ROI inside LV');
        roiResp = imellipse;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        
        Nxshift = roiPosition(1);
        Nyshift = roiPosition(2);
        ROINx = roiPosition(3);
        ROINy = roiPosition(4);

        mask = zeros(Ny,Nx);
        mask(Nyshift+(1:ROINy),Nxshift+(1:ROINx)) = hanning(ROINy)*hanning(ROINx).'; 
        mask(mask>0.1) = 1;
        mask(mask<=0.1) = 0;
        mask = repmat(mask,[1 1 Nz]);
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
        
        % for testLVvolume, experimental
        mask = zeros(Ny,Nx);
        mask(floor(Ny/2),floor(Nx/2)) = 1;
        mask = repmat(mask,[1 1 Nz]);
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

    fprintf(' -- Cardiac binning ICA + clustering, ');
    
    cbins_half = ceil(cbins/2);
    cbins_full = cbins_half*2;
    cbins = cbins_full;
    
    %% Cardiac window
    dt  = lEchoSpacing*SGBlock;    
    fs  = 1/dt;
    fsm = 1/(lEchoSpacing*linesPerShot);
    df  = fs/size(Phi_rt_small,2);
    
    winlp = 2*floor((HRhigh/60)/df); %highest HR: 130 bpm
    windowlp = zeros(size(Phi_rt_small,2),1);
    windowlp(1:winlp) = 1;
    windowlp = circshift(windowlp,[-winlp/2, 0]);

    winhp = 2*floor((HRlow/60)/df); %lowest HR: 40 bpm
    windowhp = zeros(size(Phi_rt_small,2),1);
    windowhp(1:winhp) = 1;
    hwindow = windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
    
    winwidth = 0.05/(lEchoSpacing*SGBlock); % 35 ms window width
    winwidth = ceil((winwidth-1)/2)*2 + 1;   % make odd
    
    %Phi_rt_motion = Phi_rt_small - Phi_rt_small*pinv(curvePhi(Segidx,:).')*curvePhi(Segidx,:).';  %orthogonalize to Bloch subspace
    Phi_rt_motion = Phi_rt_small;
    roi_mask  = disp2d(logical(roi_weighting>0.2));
    tempmask  = repmat(roi_mask,[1,1,size(Utemp,3),size(Utemp,4)]);
    Utempmask = disp2d(Utemp);
    Utempmask = reshape(Utempmask(tempmask),sum(roi_mask(:)),L_init);
    tempmasksig = Utempmask*Phi_rt_motion;
    masksig = reshape(tempmasksig,sum(roi_mask(:)),[])';
    clear recon tempmask tempmasking

    input = masksig;
    plotrange = 1:size(masksig,1);
    input = abs(input);
    input = input/max(input(:));
    q = 5;
    Mdl = rica(input,q,'NonGaussianityIndicator',ones(q,1));
    unmixed = transform(Mdl,input);
    
    plotPhi(unmixed',lEchoSpacing*SGBlock,1:min(q,6),'ICA results');
    
    n   = size(Phi_rt_small,2);
    f(1:ceil(n/2)) = (0:ceil(n/2)-1)*(fs/n);            % frequency range
    f(ceil(n/2)+1:n) = (-floor(n/2):1:-1)*(fs/n);      % frequency range
 
    temp = unmixed;
    for fn = 1:floor(max(f)/fsm)
        d(fn) = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',fsm*fn-df/0.5,'HalfPowerFrequency2',fsm*fn+df/0.5, ...
            'DesignMethod','butter','SampleRate',fs);
        for i = 1:q
            temp(:,i)= filtfilt(d(fn),temp(:,i));
        end
    end
    plotPhi(temp',lEchoSpacing*SGBlock,1:min(q,6),'ICA filtered results');
    unmixed_decoupl = temp;
    clear temp
    
    for p = 1:q
        x = abs(unmixed_decoupl(:,p));
        y = fft(x);
        n = length(x);          % number of samples
        harmonicsm = (abs(mod(f,fsm))<0.05)|(abs(mod(f,fsm))>(fsm-0.05)); %harmocin frequency for prep modulation
        Power(p,:) = abs(y).^2/n;    % power of the DFT
        frange = (f>0.001)&(f<5);
%         figure(103);subplot(q,1,p); plot(f(frange),Power(p,frange));%/sum(power(:,p)));

        Cardiacrange = (f>(HRlow/60))&(f<(HRhigh/60)).*~harmonicsm;
        AcceptRange  =( f>(1/100))&(f<(200/60)).*~harmonicsm;
%         figure(104);
%         subplot(q,1,p); plot(f(frange),Power(p,frange).*Cardiacrange(frange));
        prepmodscore(p) = sum(Power(p,frange).*harmonicsm(frange))/sum(Power(p,AcceptRange));
        cscore(p) = sum(Power(p,frange).*Cardiacrange(frange))/sum(Power(p,frange).*AcceptRange(frange));
        %rscore(p) = sum(Power(p,frange).*respiratoryrange(frange))/sum(Power(p,frange).*AcceptRange(frange));
    end
    
    Cardiacp = find(cscore==max(cscore));
    %Respp=find(rscore==max(rscore));
    clear f Power;

    hf = designfilt('bandpassiir','FilterOrder',2, ...
        'HalfPowerFrequency1',HRlow/60,'HalfPowerFrequency2',HRhigh/60, ...
        'SampleRate',fs);
    Z = filtfilt(hf,unmixed_decoupl(:,Cardiacp));
    
    Phicard = zeros(rbins,cbins_half,L_init,ccL);
    [Hidx,Phicard,bestres,bestresnorm] = runBinningCardLoop(alpha,its,Ridx(:),Z(:),Phi_rt_small,curvePhi,Phicard,Segidx,hwindow,winwidth);
    
    Hidx = Hidx(:);
    Ridx = Ridx(:);

    bestresnorm = norm(bestres);
else
    cbins = 1;
    Hidx    = ones(size(Ridx));
    bestres = ones(size(Ridx))/numel(Ridx);
end

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

dataArray.Hidx = Hidx;
dataArray.bestresCard = bestresCard;
dataArray.bestres     = bestres;
dataArray.wallClock   = wall_clock;
dataArray.meanHBPM    = round(numel(find(diff(Hidx)<0))*60/(lEchoSpacing*SGBlock*numel(Hidx)));

% dataArray.diastoleIdx = 1;
% dataArray.Hidx = mod(Hidx - dataArray.diastoleIdx, max(Hidx)) + 1;


%% Images for individual bins

if cbins > 1
    Utemp = Utemp.*windowfun;

    temp = zeros(Nydisp,Nxdisp,numel(dispSlice),cbins,rbins);
    for j = 1:rbins
        for k = 1:cbins
            temp(:,:,:,k,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                              *mean(Phi_rt_small_init(:,Hidx==k & Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[]));
        end
    end
    cw = prctile(temp(:),99.5);
    dataArray.binsCard = abs(temp)/cw;

    mask = dispim(mask);
    dataArray.diastoleIdx = testLVvolume(squeeze(abs(temp(:,:,1,:,1))),mask);
    dataArray = shiftHidx(params,reconOptions,dataArray,temporalBasis,spatialCoeff);
else
    dataArray.diastoleIdx = 1;
end


function  [Z,Phicard,res,resnorm] = runBinningCardLoop(alpha,its,Ridx,Z,Phi_rt_small,curvePhi,Phicard,Segidx,hwindow,winwidth)

    [rbins,cbins_half,L_init,ccL] = size(Phicard);
    Z = realify(Z);   
    Z = Z./abs(hilbert(Z));
    Z = ceil((Z+1)*cbins_half/2);
    Z(Z<1) = 1;
    Z(Z>cbins_half) = cbins_half;
    
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
        Z = Z(:);

        Z = fft(Z/sqrt(numel(Z))).*hwindow;
        if (alpha > 0) && (it > 10) && (max(abs(Z))>alpha)
            Z = sign(Z).*max(abs(Z)-alpha,0);
        end
        Z = realify(ifft(Z));
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