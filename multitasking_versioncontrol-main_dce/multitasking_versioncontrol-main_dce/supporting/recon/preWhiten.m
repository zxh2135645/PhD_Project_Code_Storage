function [params, reconOptions, dataArray] = preWhiten(params, reconOptions, dataArray, cropRO)

if nargin < 3
    cropRO = 1;
end
    
vec = @(x) x(:);

extractVarFromStruct(params);
extractVarFromStruct(dataArray);

Ncoils = size(kspaceData, 4);

if ~flagIsPrewhitened

    %% Noise estimation
    if Nz > 3
        if isCartesian
            t_ind = abs(linOrder-DC_ky)>=Ny/4 & abs(parOrder-DC_kz)>=Nz/4;
        else
            t_ind = abs(parOrder-DC_kz)>=Nz/4;
        end
    else
        if isCartesian
            t_ind = abs(linOrder-DC_ky)>=Ny/4 & (parOrder==1);
        else
            t_ind = (parOrder==1);
        end
    end
        
    fprintf('Estimating noise covariance matrix...')
    if exist('noiseData','var')
        Psi = reshape(permute(noiseData(2:end,:,:),[1 3 2]),[],Ncoils);
    else
        Psi = reshape(kspaceData(t_ind,end,:,:,end),[],Ncoils);
    end
    Psi   = cov(Psi);
    msdev = sqrt(trace(Psi)/Ncoils);
    Psi   = Psi/msdev^2;

    %% Pre-whiten and transform data
    NechoNav = size(navData,5);
    
    kspaceData = reshape(permute(kspaceData,[5 1 2 3 4]),[],Ncoils);
    navData    = reshape(permute(navData,[5 1 2 3 4]),[],Ncoils);
    
    fprintf('Prewhitening...')
%     temp = sqrtm(Psi);
%     kspaceData = kspaceData/temp; %/sqrtm(Psi), *inv(sqrtm(Psi)), *sqrtm(inv(Psi));
%     navData    = navData/temp;
    
    kspaceData = permute(reshape(kspaceData,Necho,[],Nkx,size(kspaceData,3),Ncoils),[2 3 4 5 1]);
    navData    = permute(reshape(navData,NechoNav,[],Nkx,size(navData,3),Ncoils),[2 3 4 5 1]);
    
    msdev = std(vec(kspaceData(t_ind,end,:,:,end)));

    fprintf('Cropping oversampled (readout) dimension...')
    newNx = floor(Norig/cropRO/2)*2;

    if ~strcmp(Trajectory,'Spiral')
        temp = zeros(size(navData,1),Nx,1,Ncoils,NechoNav);
        temp(:,(Nx-Nkx+1):Nx,:,:,:) = navData;
        navData = ift1d(temp,2);
        navData = navData(:, DC_kx + ((-newNx/2):(newNx/2-1)),:,:,:);
    end
    
    % For Cartesian, apply FFT in kx direction on kspaceData and crop
    if isCartesian      % ifft and crop out readout oversampling
        temp = zeros(size(kspaceData,1),Nx,1,Ncoils,Necho);
        temp(:,(Nx-Nkx+1):Nx,:,:,:) = kspaceData;
        kspaceData = ift1d(temp,2);
        kspaceData = kspaceData(:, DC_kx + ((-newNx/2):(newNx/2-1)), :, :, :);
        params.Nx     = newNx;
        params.Nxdisp = newNx;
        params.N     = newNx;
        params.Norig = newNx;
        params.DC    = floor(newNx/2) + 1; 
        params.DC_kx = params.DC;
        reconOptions.st.Nd(2)    = newNx;
        reconOptions.st.Ndisp(2) = newNx;
    end

    % For radial, orthogonalize nav data to GA phase variation.
    if strcmp(Trajectory,'Radial')
        thetas_temp = 180 * mod(trajInc * ((1:Nread) - 1 + (cutoff - ceil(cutoff/SGBlock))), 2*Ntrajs) / Ntrajs;
        thetas_temp = thetas_temp(cutoff_shift:SGBlock-1:end);
        thetas_temp = thetas_temp(1:size(navData,1));
        nuisance    = exp(bsxfun(@times,[1i; -1i],thetas_temp*pi/180)).';
        for echo = 1:size(navData,5)    
            temp = navData(:,:,:,:,echo);
            navData(:,:,:,:,echo) = reshape(temp(:,:) - nuisance*(pinv(nuisance)*temp(:,:)),size(temp));
        end
        clear thetas_temp
    end

    dataArray.kspaceData = kspaceData;
    dataArray.navData    = navData;
    dataArray.msdev    = msdev;
    dataArray.Psi      = Psi;
    dataArray.Psi_orig = Psi;
    dataArray.flagIsPrewhitened = true;
    
    fprintf(' done.\n')
end

% POCS
function img = pocs(tempk,dim,iterations)
asymmask = logical(abs(tempk));
N = size(tempk,dim);
newdim = [1 2 3 4 5]; newdim(1) = dim; newdim(dim) = 1;
temp = permute(tempk,newdim);
centerwinwidth = (N/2 - find(sum(sum(temp(:,:,:,1,1),2),3)==0,1,'last'))*2;
centerwin = zeros(N,1);
centerwin(N/2-centerwinwidth/2+(1:centerwinwidth)) = hamming(centerwinwidth); 
centerwin = permute(centerwin,newdim);
img = ift1d(tempk,dim);
imgphase = sign(ift1d(bsxfun(@times,tempk,centerwin),dim));
for j = 1:iterations    % POCS iterations
    temp = abs(img).*imgphase;
    temp = ft1d(temp,dim);
    temp(asymmask) = tempk(asymmask);
    img  = ift1d(temp,dim);
end
