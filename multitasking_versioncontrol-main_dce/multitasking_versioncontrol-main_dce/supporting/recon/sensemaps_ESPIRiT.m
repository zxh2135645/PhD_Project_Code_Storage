%% ESPIRiT for Multitasking w./ CAIPIRINHA sampling 
% Written by Sen Ma, BIRI, Cedars-Sinai Medical Center 05/18/2020
% Code is based on the 2D demo provided by Michale Lustig group, https://people.eecs.berkeley.edu/~mlustig/software/
% This version adds the option of 3D ESPIRiT

% Input:
    % ims: multi-channel images [Ny Nx Nz Ncoil] 
    % calib_size: k-space ACS dimension
    %       3D: calib_size = [caliby, calibz]; 
    %       2D: calib_size = [caliby];
    % kernel_size: kernel dimension
    %       3D: kernel_size = [kernely,kernelz,kernelx];
    %       2D: kernel_size = [kernely,kernelx];

% Output:
    % maps: weighted sensitivity maps, [Ny Nx Nz Ncoil] or [Ny Nx Ncoil]
    % weights: sensitivity map weighting


function [maps, weights] = sensemaps_ESPIRiT(ims, calib_size, kernel_size)

if (length(kernel_size) ~= length(size(ims))-1) || (length(kernel_size) ~= length(calib_size))
    disp('Kernel size and image size do not match!')
    return
end

maps = [];
weights = [];

dims = length(kernel_size);
if dims == 2
    ims = squeeze(ims);
    [Ny,Nx,Nc] = size(ims);
    % k_space: ky centered; kz centered; kx centered
    fbp_k = ft1d(ft1d(ims,1),2);
    
    %display magnitude coil image
    figure, imshow3(abs(ims),[],[1,Nc]); 
    title('magnitude of physical coil images');
    colormap((gray(256))); colorbar;

    %display phase coil image
    figure, imshow3(squeeze(angle(ims)),[],[1,Nc]); 
    title('phase of physical coil images');
    colormap('hsv'); colorbar;
else
    [Ny,Nx,Nz,Nc] = size(ims);
    % k_space: ky centered; kz centered; kx centered
    fbp_k = ft1d(ft1d(ft1d(ims,1),2),3);
    
    %display magnitude coil image
    figure, imshow3(squeeze(abs(ims(:,:,ceil(Nz/2),:))),[],[1,Nc]); 
    title('magnitude of physical coil images');
    colormap((gray(256))); colorbar;

    %display phase coil image
    figure, imshow3(squeeze(angle(ims(:,:,ceil(Nz/2),:))),[],[1,Nc]); 
    title('phase of physical coil images');
    colormap(hsv(256)); colorbar;
end


if dims == 2
    eigThresh_k  = 0.02; % threshold of eigenvectors in k-space
    eigThresh_im = 0.9; % threshold of eigenvectors in image space
    calib = crop(fbp_k,calib_size(1),Nx,Nc); % ACS region for calibration
else
    eigThresh_k  = 0.05; % threshold of eigenvectors in k-space
    eigThresh_im = 0.9; % threshold of eigenvectors in image space
    calib = crop(fbp_k,calib_size(1),calib_size(2),calib_size(3),Nc); % ACS region for calibration
end

%% Compute Eigen-Value Maps
% Maps are computed in two steps. 

% First, compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
disp('STEP 1: Compute calibration matrix and k-space kernels (row space)...')
tic;[k,S] = dat2Kernel(calib,kernel_size);toc;

idx = max(find(S >= S(1)*eigThresh_k));

% Display the singular vectors and values of the calibration matrix

kdisp = reshape(k,[prod(kernel_size)*Nc,prod(kernel_size)*Nc]);
figure, subplot(211), plot([1:prod(kernel_size)*Nc],S,'LineWidth',2);
hold on, 
plot([1:prod(kernel_size)*Nc],S(1)*eigThresh_k,'r-','LineWidth',2);
plot([idx,idx],[0,S(1)],'g--','LineWidth',2)
legend('signular vector value','threshold')
title('Singular Vectors')
subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
xlabel('Singular value #');
title('Singular vectors')

%%
% Second, crop kernels and compute eigen-value decomposition in image space to get
% maps
disp('STEP 2: Compute eigenvalue decompostions in image space...')
if dims == 3
    tic;    [M,W] = kernelEig(k(:,:,:,:,1:idx),[Ny Nx Nz]);  toc;
else
    tic;    [M,W] = kernelEig(k(:,:,:,1:idx),[Ny Nx]);  toc;
end

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps
if dims == 3
    Wtemp = squeeze(W(:,:,ceil(Nz/2),:));
    Mtemp = squeeze(M(:,:,ceil(Nz/2),:,:));
else
    Wtemp = W;
    Mtemp = M;
end

figure, imshow3(abs(Wtemp),[],[1,Nc]); 
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

figure, imshow3(abs(Mtemp),[],[Nc,Nc]); 
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(Mtemp),[],[Nc,Nc]); 
title('Phase of Eigen Vectors');
colormap(hsv(256)); colorbar;

clear Wtemp Mtemp

%%
% crop sensitivity maps 
if dims == 3
    %maps = M(:,:,:,:,end).*repmat(W(:,:,:,end)>eigThresh_im,[1,1,1,Nc]);
    maps = M(:,:,:,:,end);
    weights = repmat(W(:,:,:,end)>eigThresh_im,[1,1,1,Nc]);
    maps_temp = squeeze(maps(:,:,ceil(Nz/2),:));
else
    maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);
    weights = repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);
    maps_temp = maps;
end

figure, imshow3(abs(maps_temp),[],[1,Nc]); 
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps_temp),[],[1,Nc]); 
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;

clear maps_temp



