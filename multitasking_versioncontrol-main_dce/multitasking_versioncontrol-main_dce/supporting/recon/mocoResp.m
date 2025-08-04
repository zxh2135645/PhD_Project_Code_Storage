function [params,reconOptions,dataArray,temporalBasis,spatialCoeff] = mocoResp(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI)
%v0.2
if nargin < 6
    selectROI = false;
end

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L_init = size(Phi_rt_small_init,1);
cL = size(curvePhi,2);

vec = @(x) x(:);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end

templateim = zeros(Ny*Nx*Nz,rbins);
for j = 1:rbins
    templateim(:,j) = reshape(Utemp,[],L_init)*mean(Phi_rt_small_init(:,Ridx(:)==j),2);
end
templateim = reshape(templateim,Ny,Nx,Nz,[]);
templateim = templateim(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:,:);  % crop to center
% if flagCommandLine
%     implay(abs(templateim(:,:,DC_kz,:))/max(abs(vec(templateim(:,:,DC_kz,:)))),2);
% end

% Define k-space coordinates
[ky,kx,kz] = ndgrid(-Nydisp/2:Nydisp/2-1,-Nxdisp/2:Nxdisp/2-1,-floor(Nz/2):ceil(Nz/2)-1);
ky = ifftshift(ky)*2*pi/Nydisp;
kx = ifftshift(kx)*2*pi/Nxdisp;
kz = ifftshift(kz)*2*pi/Nz;

% Fourier transform template images
if Nz > 1
    fft3  = @(x) fft(fft2(x),[],3);
    ifft3 = @(x) ifft(ifft2(x),[],3);
else
    fft3  = @(x) fft2(x);
    ifft3 = @(x) ifft2(x);
end

% Define window
if selectROI
    temp = fbpComposite(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:,:);
    h = figure;imagesc(abs(temp(:,:,floor(Nz/2)+1,1)));axis equal tight;colormap('gray');title('Draw heart ROI')
    roiResp = imellipse;
    roiPosition = roiResp.getPosition();
    roiPosition = round(roiPosition);
    close(h);
    Nxshift = roiPosition(1);
    Nyshift = roiPosition(2);
    ROINx   = roiPosition(3);
    ROINy   = roiPosition(4);

    if Nz > 1
        h = figure;imagesc(squeeze(abs(temp(:,Nxshift+floor(ROINx/2),:,1)))');axis equal tight;colormap('gray');title('Draw heart ROI')
        roiResp = imellipse;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        Nzshift = roiPosition(1);
        Nyshift = min(roiPosition(2),Nyshift);
        ROINz   = roiPosition(3);
        ROINy   = max(roiPosition(4),ROINy);

        h = figure;imagesc(squeeze(abs(temp(Nyshift+floor(ROINy/2),:,:,1))));axis equal tight;colormap('gray');title('Draw heart ROI')
        roiResp = imellipse;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        Nzshift = min(roiPosition(1),Nzshift);
        Nxshift = min(roiPosition(2),Nxshift);
        ROINz   = max(roiPosition(3),ROINz);
        ROINx   = max(roiPosition(4),ROINx);
    else
        Nzshift = 0;
        ROINz   = 1;
    end
    
    windowfun = zeros(Nydisp,Nxdisp,Nz);
    windowfun(Nyshift+(1:ROINy),Nxshift+(1:ROINx),:) = repmat(hanning(ROINy)*hanning(ROINx).',[1 1 Nz]);
    windowZ = zeros(Nz,1);
    windowZ(Nzshift+(1:ROINz)) = hann(ROINz);
    windowfun = sqrt(sqrt(bsxfun(@times,windowfun,reshape(windowZ,1,1,Nz))));
else
    Nywindow = floor(3*Nydisp/4);
    Nxwindow = floor(3*Nxdisp/4);
    windowfun = zeros(Nydisp,Nxdisp);
    windowfun(floor(Nydisp/2)-floor(Nywindow/2) + (1:Nywindow), floor(Nxdisp/2)-floor(Nxwindow/2) + (1:Nxwindow)) = hann(Nywindow)*hann(Nxwindow).';
    windowfun = sqrt(sqrt(bsxfun(@times,windowfun,reshape(hann(Nz),1,1,Nz))));
end


%% Register images to template image (max correlation)
fprintf('Registering images to the 1st bin... ');
templateim = fft3(templateim);
xyzshifts = zeros(rbins,3);
opts = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
shiftkernel = @(x)exp(1i*(ky*x(1)+kx*x(2)+kz*x(3)));
targetim = vec(bsxfun(@times,ifft3(templateim(:,:,:,1)),windowfun));

% Register images to the last bin image (max correlation)
% for j = rbins-1:-1:1
%     cost = @(x)1-abs(corr( targetim, vec(bsxfun(@times,ifft3(templateim(:,:,:,j).*shiftkernel(x)),windowfun)) ));
%     xyzshifts(j,:) = fminsearch(cost,xyzshifts(j+1,:),opts);
% end
% xyzshifts(end,:) = [];

% Register images to the first bin image (max correlation)
for j = 2:rbins
    cost = @(x)1-abs(corr( targetim, vec(bsxfun(@times,ifft3(templateim(:,:,:,j).*shiftkernel(x)),windowfun)) ));
    xyzshifts(j,:) = fminsearch(cost,xyzshifts(j-1,:),opts);
end
xyzshifts(1,:) = [];

%% Do fully joint registration?
resh4 = @(x)reshape(x,1,1,1,[]);
%shiftkernel = @(x)exp(1i*(bsxfun(@times,kx,resh4([x(:,1); 0]))+bsxfun(@times,ky,resh4([x(:,2); 0]))+bsxfun(@times,kz,resh4([x(:,3); 0]))));
shiftkernel = @(x)exp(1i*(bsxfun(@times,ky,resh4([0; x(:,1)]))+bsxfun(@times,kx,resh4([0; x(:,2)]))+bsxfun(@times,kz,resh4([0; x(:,3)]))));
cost = @(x)sum(svde(reshape(bsxfun(@times,ifft3(templateim.*shiftkernel(x)),windowfun),[],rbins)));
opts = optimset('PlotFcns',@optimplotfval);
xyzshifts = fminsearch(cost,xyzshifts,opts);
cost(xyzshifts);

%% View result
fprintf('displaying results.\n');
templateim_reg = imageOrientLPS(ifft3(templateim.*shiftkernel(xyzshifts)),params);
templateim     = imageOrientLPS(ifft3(templateim),params);

% Relative slice thickness
slthick = dThickness_mm*Norig/dReadoutFOV_mm;
newNz   = round(slthick*Nz);

if reconOptions.flagCommandLine
    if Nz > 1
        if abs(zDir) == 3 
            vecOrder = abs([zDir xDir 4 yDir]);
        else
            vecOrder = [3 2 4 1];
        end
        templateim2_reg = permute(templateim_reg,vecOrder);
        templateim2     = permute(templateim,vecOrder);
        if abs(zDir) == 3 && sign(zDir) > 0
            templateim2_reg = flip(templateim2_reg,1);
            templateim2     = flip(templateim2,1);
        end
        ky = floor(size(templateim2,4)/2) + 1;
        temp = [squeeze([templateim(:,:,DC_kz,:) templateim_reg(:,:,DC_kz,:) templateim(:,:,DC_kz,:)-templateim_reg(:,:,DC_kz,:)]);
                imresize([templateim2(:,:,:,ky) templateim2_reg(:,:,:,ky) templateim2(:,:,:,ky)-templateim2_reg(:,:,:,ky)],[round(newNz) Nxdisp*3])];
        temp = abs(temp)/max(abs(temp(:)));
    else
        temp = abs([templateim(:,:,DC_kz,:) templateim_reg(:,:,DC_kz,:) templateim(:,:,DC_kz,:)-templateim_reg(:,:,DC_kz,:)])/max(abs(vec(templateim(:,:,DC_kz,:))));
    end
    implayZoom(temp*1.2,rbins/2,2,'Resp bins. Left: original, Mid: corrected, Right: Difference');
    % implayZoom(abs([templateim(:,:,DC_kz,:) templateim_reg(:,:,DC_kz,:) templateim(:,:,DC_kz,:)-templateim_reg(:,:,DC_kz,:)])/max(abs(vec(templateim(:,:,DC_kz,:)))),rbins/2);
    % implay(abs(interpft([permute(templateim(Norig/2+1,:,:,:),[3 2 4 1]) permute(templateim_reg(Norig/2+1,:,:,:),[3 2 4 1]) permute(templateim(Norig/2+1,:,:,:)-templateim_reg(Norig/2+1,:,:,:),[3 2 4 1])],round(slthick*Nz)))/max(abs(vec(templateim(Norig/2+1,:,:,:)))))
    % implay(abs(interpft([permute(templateim(:,Norig/2+1,:,:),[3 1 4 2]) permute(templateim_reg(:,Norig/2+1,:,:),[3 1 4 2]) permute(templateim(:,Norig/2+1,:,:)-templateim_reg(:,Norig/2+1,:,:),[3 1 4 2])],round(slthick*Nz)))/max(abs(vec(templateim(:,Norig/2+1,:,:)))))
end

% disp('Check MoCo quality')
% keyboard;

%% Correct data

xyzshifts = [0 0 0;xyzshifts];

% Correct nav data
fprintf('Correcting navData... ');

xyzshifts_nav = zeros(size(navData,1),3);
for j = 2:rbins
    xyzshifts_nav(Ridx==j,:) = repmat(xyzshifts(j,:),[sum(Ridx==j), 1]);
end

navk = fft(navData,[],2);
k = -size(navk,2)/2:size(navk,2)/2-1;
k = ifftshift(k)*2*pi/size(navk,2);
navk = bsxfun(@times,navk,exp(1i*xyzshifts_nav(:,2)*k(:).'));
navdata_reg = ifft(navk,[],2);
try
    dataArray.navData = navdata_reg(:,DC+((-Norig/2):(Norig/2-1)),:,:,:);
catch
    dataArray.navData = navdata_reg(:,ceil(max(abs(xyzshifts(:,2)))):end-floor(max(abs(xyzshifts(:,2)))),:,:,:);
end

%% Correct k-space data
fprintf('correcting kspaceData...');

xyzshifts_k = interp1(navIndices,xyzshifts_nav,1:Ntpoint,'nearest','extrap');
xyzshifts_k(navIndices,:)=[];
xyzshifts_k = [ sin(thetas(:)*pi/180).*xyzshifts_k(:,1) ...
              + cos(thetas(:)*pi/180).*xyzshifts_k(:,2), ...
                xyzshifts_k(:,3) ];

k = (1-DC_kx):(Nkx-DC_kx);
%k  = -Norig:Norig-1;
k  = k*pi/Norig;
kz = (parOrder-DC_kz)*2*pi/Nz;

dataArray.kspaceData = bsxfun(@times,kspaceData,exp(1i*bsxfun(@plus,xyzshifts_k(:,1)*k(:).',xyzshifts_k(:,2).*kz(:))));

%% Recalculate sensitivites and new "initial" reconstruction

flagCommandLineTemp = flagCommandLine;
reconOptions.flagCommandLine = false;

[params, reconOptions, dataArray] = estimateSensitivities(params, reconOptions, dataArray);

close all;

[params, reconOptions, dataArray,temporalBasis] = realtimeSubspace(params, reconOptions, dataArray, temporalBasis);

reconOptions.flagUseInitialGuess = false; 
[reconOptions,dataArray,temporalBasis,spatialCoeff] = reconLeastSquares(params, reconOptions, dataArray, temporalBasis);

flagCommandLine = flagCommandLineTemp;
reconOptions.flagCommandLine = flagCommandLineTemp;

close all;

%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

Phi_rt_small_init = temporalBasis.Phi_rt_small_init;
L_init = size(Phi_rt_small_init,1);

Utemp = reshape(spatialCoeff.U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end

temp = zeros(Nydisp,Nxdisp,numel(dispSlice),rbins);
for j = 1:rbins
    temp(:,:,:,j) = abs(reshape(reshape(dispim(Utemp),[],L_init)...
                    *mean(Phi_rt_small_init(:,Ridx==j),2),Nydisp,Nxdisp,numel(dispSlice),[]));
end
cw = prctile(temp(:),99);
dataArray.binsRespMean = abs(temp)/cw;

dataArray.binsResp = cell(rbins,1);
for j = 1:rbins
    temp = abs(reshape(reshape(dispim(Utemp),[],L_init)...
           *Phi_rt_small_init(:,Ridx==j),Nydisp,Nxdisp,numel(dispSlice),[]));
    cw = prctile(temp(:),99);
    dataArray.binsResp{j} = abs(temp)/cw;
end

% h = findall(groot,'Type','figure','Name','Binning Results');
% if isempty(h)
%     figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
% else
%     figure(h);
% end
% subplot(2,2,3),plot(Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
% subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
if flagCommandLine 
    implayZoom(imageOrientLPS(dataArray.binsRespMean(:,:,1,:),params),2);
end

reconOptions.flagRespMoco = true;

