function dataArray = binningResp_edge2(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI)
%v2.0

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Invalid hash: Multitasking "license check" failed');
%     waitfor(dlg);
%     return;
% end

if nargin < 6
    selectROI = false;
end

% load variables
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

% if strcmp(params.ScanType,'SR')  % respiratory dimension for perfusion
%     rbins = 1;
% end

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

L_init = size(Phi_rt_small_init,1);
cL = size(curvePhi,2);
Segidx = mod(navIndices-1,linesPerShot*moduleLength) + 1;

vec = @(x) x(:);

if rbins > 1

    %% Prepare weighted images
    Nzshift = floor(Nz/2) - floor(Nzorig/2);
    Nyshift = floor(Ny/2) - floor(Nydisp/2);
    Nxshift = floor(Nx/2) - floor(Nxdisp/2);
   

    Phi_rt       = Phi_rt_init;
    Phi_rt_full  = Phi_rt_full_init;

    temporalBasis.Phi_rt       = Phi_rt;
    temporalBasis.Phi_rt_full  = Phi_rt_full;
%     temporalBasis.Phi_rt_small = Phi_rt_small;

    %% Realtime & masking

    Utemp = reshape(U_init,Ny,Nx,Nz,[]);
    if MBfactor == 1
        Utemp = fftshift(Utemp,3);
    end
    if params.isCartesian
        Utemp = fftshift(Utemp,1);
    end
    Phi_rt_small = Phi_rt_small_init;
    Phi_rt_small(1,:) = repmat(Phi_rt_small_init(1,round(linesPerShot*5*0.99)),1,size(Phi_rt_small_init,2));

    ntemp = ceil(totalTime/60);
    segtemp = floor(size(Phi_rt_small,2)/ntemp);
    cstep = 0;
    sig = [];
    for ctemp = 1:ntemp
        recon=permute(squeeze(abs(dispim(reshape(reshape(Utemp,Ny*Nx*Nz,[])*Phi_rt_small(:,(ctemp-1)*segtemp+1:(ctemp)*segtemp),Ny,Nx,Nz,[])))),[2 1 3]);
        cstep=cstep+1;
        sig(:,:,(cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3))=recon;
    end
    recon = permute(squeeze(abs(dispim(reshape(reshape(Utemp,Ny*Nx*Nz,[])*Phi_rt_small(:,ntemp*segtemp+1:size(Phi_rt_small,2)),Ny,Nx,Nz,[])))),[2 1 3]);
    sig = cat(3,sig,recon);
    cw = prctile(abs(sig(:)),99.5);
    sig = sig./cw;
    implay(sig);
    recon = sig;

    % Respiratiory signal extraction
    close all hidden;
    reconsl = recon(:,:,round(linesPerShot*5*0.99));
    cw = prctile(reconsl,99);
    figure;imagesc(reconsl./cw);axis equal tight;colormap('gray');title('Draw resp motion ROI')
    h = drawline('LineWidth',4,'Color','cyan');
    x = h.Position(:,1);
    y = h.Position(:,2);
    slope = -1/((y(2)-y(1))/(x(2)-x(1)));
    roiw = 5;
    xnew1 = x(1) + (roiw/2*sqrt(1/(1+slope^2)))*[-1;1];
    xnew2 = x(2) + (roiw/2*sqrt(1/(1+slope^2)))*[-1;1];
    ynew1 = y(1) + (slope*roiw/2*sqrt(1/(1+slope^2)))*[-1;1];
    ynew2 = y(2) + (slope*roiw/2*sqrt(1/(1+slope^2)))*[-1;1];
    Position = [[xnew1(1),ynew1(1)]; [xnew2(1),ynew2(1)]; [xnew2(2),ynew2(2)]; [xnew1(2),ynew1(2)]];
    h = images.roi.Freehand(gca,'Position',Position);
    respmask=createMask(h);

    roil = round(norm(Position(1,:)-Position(2,:)));
    vecl = (Position(2,:)-Position(1,:))./roil;
    vecw = (Position(4,:)-Position(1,:))./roiw;
    mcorx = zeros(roiw+1,roil+1);
    mcory = zeros(roiw+1,roil+1);
    for ii = 0: roil
        for jj = 0: roiw
            mcorx(jj+1,ii+1) = Position(1,1) + ii*vecl(1) + jj*vecw(1);
            mcory(jj+1,ii+1) = Position(1,2) + ii*vecl(2) + jj*vecw(2);
        end
    end
    [X,Y] = meshgrid(floor(min(Position(:,1),[],'all')):ceil(max(Position(:,1),[],'all')),...
                                  floor(min(Position(:,2),[],'all')):ceil(max(Position(:,2),[],'all')));
    dim = size(X);
    mind = sub2ind(size(reconsl),Y(:),X(:));
    respsig = zeros(roil+1,size(recon,3));
    for i = 1:size(recon,3)
        temp = recon(:,:,i);
        V = temp(mind);
        V = reshape(V,dim(1),[]);
        respsig(:,i) = mean(interp2(X,Y,V,mcorx,mcory),1);
    end
    figure; imagesc(respsig,'CDataMapping','scaled');
    colormap('gray');

    %respsig = histeq(respsig);
    %respsig(respsig  < 0.2) = 0;
    %figure; imagesc(respsig,'CDataMapping','scaled'); colormap('gray'); 
%{
    pause;
    prompt = {'Clamp at : '};
    dlgtitle = 'Input';
    dims = [1 35];
    id = inputdlg(prompt,dlgtitle,dims);

    if str2num(id{1,1}) <= size(respsig,1)
        respsig = respsig(1:str2num(id{1, 1}),:);
    else
        error(['Input Error']);
    end
%}
    %% Gating
    
    q=3;
    Mdl = rica(respsig',q,'NonGaussianityIndicator',ones(q,1));
    unmixed = transform(Mdl,respsig');

    figure(101);
    set(gcf,'Name','ICA results');
    for n=1:q
        figure(101); subplot(q,1,n); plot(unmixed(:,n));  hold on
    end

    temp = unmixed;
    dt  = lEchoSpacing*SGBlock;    
    fs  = 1/dt;
    fsm = 1/(lEchoSpacing*linesPerShot);
    df  = fs/size(Phi_rt_small,2);
    n   = size(Phi_rt_small,2);
    f(1:ceil(n/2)) = (0:ceil(n/2)-1)*(fs/n);            % frequency range
    f(ceil(n/2)+1:n) = (-floor(n/2):1:-1)*(fs/n);       % frequency range


    for fn=1:floor(max(f)/fsm)
        d(fn) = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',fsm*fn-df/0.5,'HalfPowerFrequency2',fsm*fn+df/0.5, ...
            'DesignMethod','butter','SampleRate',fs);
        for i=1:q
            temp(:,i)= filtfilt(d(fn),temp(:,i));
        end
    end
    figure(101); set(gcf,'Name',['ICA notch filtered results - ', num2str(ctemp)]); hold on
    for i=1:q
        subplot(q,1,i); plot(temp(:,i));
    end
    unmixed_decoupl=temp;

    for p=1:q
        x=abs(unmixed_decoupl(:,p));
        y = fft(x);
        n = length(x);          % number of samples
        f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
        f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
        Power(p,:) = abs(y).^2/n;    % power of the DFT
        frange=(f>0.001)&(f<5);
        figure(103); hold on
        subplot(q,1,p); plot(f(frange),Power(p,frange));%/sum(power(:,p)));

        respiratoryrange=(f>(BRlow/60))&(f<(BRhigh/60));
        rscore(p)=sum(Power(p,frange).*respiratoryrange(frange))/sum(Power(p,frange));
    end

    Respp=find(rscore==max(rscore));
    
    hf = designfilt('bandpassiir','FilterOrder',2, ...
        'HalfPowerFrequency1',BRlow/60,'HalfPowerFrequency2',BRhigh/60, ...
        'SampleRate',fs);
    Z= filtfilt(hf,unmixed_decoupl(:,Respp));
    %}
    % Method 1 : Min/Max location
    
%     locsp=ampd(Z);
%     locsv=ampd(-Z);
% 
%     x=[1,locsv,locsp];
%     y=[1;ones(length(locsv),1);ones(length(locsp),1)*rbins;];
%     xi=1:size(Z,1);
%     [xs,xsorder]=sort(x);
%     y=y(xsorder);
%     Ridx = interp1q(xs',y,xi');
%     Ridx=round(Ridx);
%     Ridx = round(interp1(find(~isnan(Ridx)),Ridx(find(~isnan(Ridx))),1:numel(Ridx),'linear','extrap'));
% 
%     %temp = circshift(unmixed(:,Respp),floor(length(unmixed(:,Respp))/2));
%     temp = circshift(Z,floor(length(Z)/2));
%     %zt = filtfilt(hf,temp);
%     zt = temp;
%     %range = sort(mod([min([locsp,locsv]),max([locsp,locsv])]+floor(length(unmixed(:,Respp))/2),length(unmixed(:,Respp))));
%     range = sort(mod([min([locsp,locsv]),max([locsp,locsv])]+floor(length(Z)/2),length(Z)));
% 
%     locspt=ampd(zt);
%     locsvt=ampd(-zt);
% 
%     x=[1,locsvt,locspt];
%     y=[1;ones(length(locsvt),1);ones(length(locspt),1)*rbins;];
%     xi=1:size(Z,1);
%     [xs,xsorder]=sort(x);
%     y=y(xsorder);
%     Ridxt = interp1q(xs',y,xi');
%     Ridxt=round(Ridxt);
%     %Ridx(1:min([locsp,locsv])) = Ridxt(floor(length(unmixed(:,Respp))/2)+1:range(2));
%     %Ridx(max([locsp,locsv]):end) = Ridxt(range(1):floor(length(unmixed(:,Respp))/2));
%     Ridx(1:min([locsp,locsv])) = Ridxt(floor(length(Z)/2)+1:range(2));
%     Ridx(max([locsp,locsv]):end) = Ridxt(range(1):floor(length(Z)/2));
%     figure; plot(Ridx);
    %}
    % Method 2 : Histogram
    
    Z = (Z-min(Z))/range(Z)*(rbins-1)+1;
    [Zn,~] = hist(Z*10,1:round(max(Z)*10));
    Zn = cumsum(Zn)/sum(Zn);
    Z = Zn(round(Z*10));
    Z = ceil(Z*rbins);
    Z(Z<1) = 1;
    Z(Z>rbins) = rbins;
    Z = Z(:);
    Ridx = Z;
    figure; plot(Ridx);

    % Method 3 : SVD
%     [~,curveS,curveU] = svd(double(respsig),'econ');
%     Z = curveU(:,1);
%     fs = 1/2/lEchoSpacing; % sample frequency (Hz)
%     hf = designfilt('bandpassiir','FilterOrder',2, ...
%         'HalfPowerFrequency1',BRlow/60,'HalfPowerFrequency2',BRhigh/60, ...
%         'SampleRate',fs);
%     Z= filtfilt(hf,Z);
%     Z = (Z-min(Z))/range(Z)*(rbins-1)+1;
%     [Zn,~] = hist(Z*10,1:round(max(Z)*10));
%     Zn = cumsum(Zn)/sum(Zn);
%     Z = Zn(round(Z*10));
%     Z = ceil(Z*rbins);
%     Z(Z<1) = 1;
%     Z(Z>rbins) = rbins;
%     Z = Z(:);
%     Ridx = Z;
%     figure; plot(Ridx);
    %}

    clear BW
end

dataArray.Ridx   = (Ridx)';
dataArray.Segidx = (Segidx)';

temp = fft(Ridx);
temp(1) = 0;
[~,idx] = max(abs(temp));
dataArray.meanRPeriod = (lEchoSpacing*SGBlock*numel(Ridx))/idx;

%% Calculate residual

Phiresp = zeros(rbins,L_init,cL);
for j = 1:rbins
    Phiresp(j,:,:) = Phi_rt_small(:,Ridx==j) * pinv(curvePhi(Segidx(Ridx==j),:).');
end
bestres = sqrt(sum(abs(Phi_rt_small.' - sum(Phiresp(Ridx,:,:).*permute(curvePhi(Segidx,:),[1 3 2]),3)).^2,2));
resnorm = norm(bestres);
bestresResp = accumarray(Ridx(:),bestres,[],@(x)norm(x)/sqrt(numel(x)));
dataArray.bestresResp = bestresResp;

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
end
cw = prctile(temp(:),99.9);
dataArray.binsRespMean = abs(temp)/cw;

temp2 = cell(rbins,1);
for j = 1:rbins
    temp2{j,1} = abs(reshape(reshape(dispim(Utemp),[],L_init)*Phi_rt_small_init(:,Ridx==j),Nydisp,Nxdisp,numel(dispSlice),[]));
    temp2{j,1} = temp2{j,1}./prctile(temp2{j,1}(:),99.9);
end
dataArray.binsResp = temp2;

% dataArray.binsResp = cell(rbins,1);
% for j = 1:rbins
%     temp = abs(reshape(reshape(dispim(Utemp),[],L_init)...
%                     *Phi_rt_small_init(:,Ridx==j),Nydisp,Nxdisp,numel(dispSlice),[]));
%     cw = prctile(temp(:),99);
%     dataArray.binsResp{j} = abs(temp)/cw;
% end

if flagCommandLine && rbins > 1
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
    
    implayZoom(dataArray.binsRespMean(:,:,1,:),2);
    implayZoom(dataArray.binsResp{max(Ridx(:))});
end

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
