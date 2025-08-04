%% ESPIRiT for Multitasking w./ CAIPIRINHA sampling 
% Written by Sen Ma, BIRI, Cedars-Sinai Medical Center 05/18/2020
% Code is based on the 2D demo provided by Michale Lustig group, https://people.eecs.berkeley.edu/~mlustig/software/
% This version adds the option of 3D ESPIRiT


% load fbp, which is generated from the ACS region of CAIPI sampling pattern
% fbp: Ny shifted; Nz shifted, phase modulated
load fbp
[Ny,Nz,Nx,Nc] = size(fbp);

% fbp: Ny centered; Nz centered, phase corrected
fbp=fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(fft(fft(fbp,[],1),[],2),1),2),[],1),[],2),1),2);

% k_space: ky centered; kz centered; kx centered
fbp_k=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(fbp,1),2),3),[],1),[],2),[],3),1),2),3)/sqrt(Ny*Nz*Nx);

%display magnitude coil image
figure, imshow3(squeeze(abs(fbp(:,end/2+1,:,:))),[],[1,Nc]); 
title('magnitude of physical coil images');
colormap((gray(256))); colorbar;

%display phase coil image
figure, imshow3(squeeze(angle(fbp(:,end/2+1,:,:))),[],[1,Nc]); 
title('phase of physical coil images');
colormap('hsv'); colorbar;

% ACS size for ky and kz (kx fully sampled)
ncaliby = 2*ceil(Ny/16); ncalibz = 2*ceil(Nz/16);
% ESPIRiT kernel-window-size
kernel_size = [6,3,6];

eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space

calib = crop(fbp_k,ncaliby,ncalibz,Nx,Nc); % ACS region for calibration

%% Compute Eigen-Value Maps
% Maps are computed in two steps. 

% First, compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
disp('STEP 1: Compute calibration matrix and k-space kernels (row space)...')
tic;[k,S] = dat2Kernel(calib,kernel_size);toc;

idx = max(find(S >= S(1)*eigThresh_k));

% Display the singular vectors and values of the calibration matrix


kdisp = reshape(k,[kernel_size(1)*kernel_size(2)*kernel_size(3)*Nc,kernel_size(1)*kernel_size(2)*kernel_size(3)*Nc]);
figure, subplot(211), plot([1:kernel_size(1)*kernel_size(2)*kernel_size(3)*Nc],S,'LineWidth',2);
hold all, 
plot([1:kernel_size(1)*kernel_size(2)*kernel_size(3)*Nc],S(1)*eigThresh_k,'r-','LineWidth',2);
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
tic;[M,W] = kernelEig(k(:,:,:,:,1:idx),[Ny Nz Nx]);toc;

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps
Wtemp = squeeze(W(:,end/2+1,:,:));

figure, imshow3(abs(Wtemp),[],[1,Nc]); 
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

Mtemp = squeeze(M(:,end/2+1,:,:,:));

figure, imshow3(abs(Mtemp),[],[Nc,Nc]); 
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(Mtemp),[],[Nc,Nc]); 
title('Phase of Eigen Vectors');
colormap(hsv(256)); colorbar;

clear Wtemp Mtemp

%%
% crop sensitivity maps 
maps = M(:,:,:,:,end).*repmat(W(:,:,:,end)>eigThresh_im,[1,1,1,Nc]);

maps_temp=squeeze(maps(:,end/2+1,:,:));

figure, imshow3(abs(maps_temp),[],[1,Nc]); 
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps_temp),[],[1,Nc]); 
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;

clear maps_temp