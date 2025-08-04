function [dataArray,temporalBasis] = binningCardSelfGate(params,reconOptions,dataArray,temporalBasis,spatialCoeff)
%v0.2.1

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

if strcmp(ScanType,'CEST')
    curvePhi = curvePhi_binning;
end

L_init = size(Phi_rt_small_init,1);
%U_init = spatialCoeff.U_init;
ccL = min(10,size(curvePhi,2)); % increase it to improve the binning quality, 8 etc. Vincent X. Mao 02/13/20
Segidx = mod(navIndices-1,size(curvePhi,1))+1;

if sum(abs(diff(diff(navTiming)))) < 1e-6
    dt = navTiming(2) - navTiming(1);    
    frame_rate  = 1/dt;
else
    newTiming = 0:floor(lEchoSpacing*SGBlock/2*1e3)/1e3:max(navTiming);
    Phi_rt_small_init_orig = Phi_rt_small_init;
    Phi_rt_small_init = interp1(navTiming,realify(Phi_rt_small_init.','cols'),newTiming,'pchip','extrap').';
    curvePhi = interp1(navTiming,realify(curvePhi(Segidx,:),'cols'),newTiming,'pchip','extrap');
    temp = [0 diff(mod(navIndices-1,linesPerShot)+1)];
    idx = find(temp<0);
    for n = 1:numel(idx)
        zeroidx = intersect(find(newTiming>navTiming(idx(n)-1)),find(newTiming<navTiming(idx(n))));
        Phi_rt_small_init(:,zeroidx) = 0;
        curvePhi(zeroidx,:) = 0;
    end
    Segidx = 1:size(Phi_rt_small_init,2);
    Ridx   = Ridx_full;
    dt = floor(lEchoSpacing*SGBlock/2*1e3)/1e3;
    frame_rate = 1/dt;
end

vec = @(x) x(:);
row = @(x) x(:).';

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)+1, :);
ovs_1 = floor(Ny/2)-floor(Nydisp/2);
ovs_2 = ceil(Ny/2)-ceil(Nydisp/2);

if cbins > 1
    if ~strcmp(ScanType,'Cine')
        h = figure;
        imshow(dispim(sqrt(abs(fbpComposite))),[]);
        title('Select Cardiac ROI')
        card_roi = imrect;
        card_roi = padarray(padarray(createMask(card_roi),[ovs_1 ovs_1],'pre'),[ovs_2 ovs_2],'post');
        close(h);
        drawnow;

        Nseg = linesPerShot;
        temp = reshape(U_init,[],size(Phi_rt_small_init,1));

        %temp = Phi_rt_small_init.'*temp(card_roi,:).';
        Phi_rt_motion = Phi_rt_small_init - Phi_rt_small_init * pinv(curvePhi(Segidx,:).') * curvePhi(Segidx,:).'; %orthogonalize to Bloch subspace
        temp = Phi_rt_motion.'*temp(card_roi,:).';

        fs = linspace(-frame_rate/2,frame_rate/2,floor(size(temp,1)/2)*2+1);
        if mod(size(temp,1),2)==0 %iseven
            fs = fs(1:(end-1));
        end
    %     notchfilt = mod(fs,1/(Nseg*lEchoSpacing));
    %     notchfilt = notchfilt < (fs(2)-fs(1)) | notchfilt > (1/(Nseg*lEchoSpacing)-(fs(2)-fs(1)));
    %     notchfilt = ifftshift(~notchfilt);
    %     temp = ifft(fft(temp).*repmat(notchfilt(:),[1 size(temp,2)]));
        [~,~,Hidx,~,wallClock] = binningSelfGate(temp,[HRlow HRhigh]/60,cbins,[BRlow BRhigh]/60,1,frame_rate,Nseg/2);
    else
        Nseg = 2;
        navData_temp = navData(:,:);
    %     navdata_temp=reshape(nav_data(:,:,:,1),size(nav_data,1),[]); %assumes 1st "coil" is virtual coil covering cardiac region
        [~,~,Hidx,~,wallClock] = binningSelfGate(navData_temp,[HRlow HRhigh]/60,cbins,[BRlow BRhigh]/60,1,frame_rate,Nseg/2);
    end
else
    Hidx = ones(size(Ridx));
end

dataArray.meanHBPM = round(numel(find(diff(Hidx)<0))*60/(dt*numel(Hidx)));

if sum(abs(diff(diff(navTiming)))) >= 1e-6
    dataArray.Hidx_full = row(Hidx);
    Hidx = interp1(newTiming,Hidx,navTiming,'nearest','extrap');
    Phi_rt_small_init = Phi_rt_small_init_orig;
    Segidx = mod(navIndices-1,size(curvePhi_binning,1))+1;
end

ccL = min(10,size(curvePhi_binning,2));
Phicard = zeros(rbins,cbins,L_init,ccL);
for j = 1:rbins
    for k = 1:cbins
        Phicard(j,k,:,:) = Phi_rt_small_init(:,(dataArray.Ridx==j)&(Hidx==k))*pinv(curvePhi_binning(Segidx((dataArray.Ridx==j)&(Hidx==k)),1:ccL).');
    end
end

bestres = zeros(size(Hidx));
for j = 1:size(Phi_rt_small_init,2)
    bestres(j) = norm(Phi_rt_small_init(:,j)-squeeze(Phicard(dataArray.Ridx(j),Hidx(j),:,:))*curvePhi_binning(Segidx(j),1:ccL).');
end
    
bestresCard = accumarray(vec(Hidx(dataArray.Ridx==1)),vec(bestres(dataArray.Ridx==1)),[],@(x)norm(x)/sqrt(numel(x)));


%% wall cloack
switch ScanType
  case 'Cine'
    wallClock = 1 + cumsum(diff(Hidx)==(1-cbins));
    wallClock(end+1) = wallClock(end);
    wallClock(wallClock == wallClock(end)) = wallClock(end)-1;
    Segidx(:) = 1;
  case {'IR','SR','T2IR','T2IR_VFA','IR_VFA'}
    wallClock = 1+cumsum(diff(Hidx)==(1-cbins));
    wallClock(end+1) = wallClock(end);
    wallClock(wallClock == wallClock(end)) = wallClock(end)-1; %the last heartbeat may be incomplete, so group it with the next-to-last
    wallClock = max(wallClock-1,1); %the first heartbeat may be incomplete, so group it with the second
  case 'T2prep'
    wallClock = 1 + cumsum(diff(Hidx)==(1-cbins));
    wallClock(end+1) = wallClock(end);
    wallClock = ceil(wallClock/5);
    wallClock(wallClock == wallClock(end)) = wallClock(end)-1;
end

%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

U_init = reshape(U_init,Ny,Nx,Nz,[]);
if isCartesian
    U_init = fftshift(fftshift(U_init,1),3);
end
temp = zeros(Nydisp,Nxdisp,numel(dispSlice),cbins,rbins);
for j = 1:rbins
    for k = 1:cbins
        temp(:,:,:,k,j) = abs(reshape(reshape(dispim(U_init),[],L_init)...
                          *mean(Phi_rt_small_init(:,Hidx==k & dataArray.Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[]));
    end
end
cw = prctile(temp(:),99);
dataArray.binsCard = abs(temp)/cw;

if flagCommandLine && cbins > 1
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(dataArray.Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(Segidx(:), dataArray.Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
    subplot(2,2,1),plot(Hidx,'.-');axis([-inf inf 0 cbins+1]);title(sprintf('Hidx, mean BPM %d',dataArray.meanHBPM));
    subplot(2,2,2);
    for rphase = 1:rbins
        plot(Segidx((dataArray.Ridx==rphase)), Hidx((dataArray.Ridx==rphase)),'.');hold on;
    end
    hold off; axis([0 linesPerShot*moduleLength 0 cbins+1]);title('Hidx / Segment Index');
    
    implayZoom(squeeze(imageOrientLPS(abs(dataArray.binsCard(:,:,1,:,1)),params)),cbins);
end

%%

dataArray.Hidx = Hidx;
dataArray.bestresCard = bestresCard;
dataArray.bestres = bestres;
dataArray.wallClockCardiac = wallClock;

dataArray.diastoleIdx = 1;
dataArray.systoleIdx  = floor(max(Hidx)/2) + 1;

temporalBasis.Phi_rt = Phi_rt;
temporalBasis.Phi_rt_full  = Phi_rt_full;
temporalBasis.Phi_rt_small = Phi_rt_small;


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
