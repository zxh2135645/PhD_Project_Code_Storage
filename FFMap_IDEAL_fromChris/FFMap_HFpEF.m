clear all;
close all;

%%
clear img_grappa
for echo = 1:size(kspace,5)
    img_grappa(:,:,:,:,echo) = fftshift(fftshift(fftshift(squeeze(ifft(ifftshift(ifft(ifftshift(ifft(ifftshift(squeeze(kspace(:,:,:,:,echo)),2),[],2),3),[],3),4),[],4)),4),3),2);
end

% img_grappa = zeros(size(kspace,2), size(kspace,3), size(kspace,4), size(kspace,5));
% for echo = 1:size(kspace,5)
%     img_grappa(:,:,:,echo) = fftshift(fftshift(fftshift(squeeze(sqrt(sum(ifft(ifftshift(ifft(ifftshift(ifft(ifftshift(squeeze(kspace(:,:,:,:,echo)),2),[],2),3),[],3),4),[],4),1))),3),2),1);
% end

%%
RawdataFT_5D = img_grappa;
iField = permute(RawdataFT_5D,[2 3 4 5 1])*100; % NumRO, NumPE, NumSlices, NumCh, NumEcho =>  NumRO, NumPE, NumSlices, NumEcho, NumCh

Mask_pre = ones((size(iField,1)-2), (size(iField, 2)-2), size(iField, 3));
Mask = zeros((size(iField,1)), (size(iField, 2)), size(iField, 3));
Mask(2:end-1, 2:end-1, :) = Mask_pre;

if size(iField,4)>1
    % combine multiple coils together, assuming the coil is the fifth dimension
    iField = sqrt(sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5));  %
    mag_iField = abs(iField);
    mag_iField_norm = mag_iField ./ max(mag_iField(:));

    Fieldmap_eddy  = sqrt(sqrt((iField(:,:,:,2)./iField(:,:,:,1))./(iField(:,:,:,3)./iField(:,:,:,2))));
    iField_uneddy = zeros(size(iField));
    iField_uneddy(:,:,:,1:2:end) = iField(:,:,:,1:2:end).*Fieldmap_eddy;
    iField_uneddy(:,:,:,2:2:end) = iField(:,:,:,2:2:end)./Fieldmap_eddy;
    iField = mag_iField_norm.*exp(1i*angle(iField_uneddy));
end

%%
voxel_size = [2.1 2.1 2.1];
f_central = 123.2*10e6;
SUBSAMPLE = 1;
TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
dfat = -3.47e-6*f_central;
LABEL = 1;

[water,fat,iFreq,unwph_uf,unwph,N_std,R2s] = ...
    spurs_gc(iField(:,:,:,:),TE,f_central,voxel_size, Mask>0,SUBSAMPLE,dfat,LABEL);

%iMag = sqrt(sum(abs(iField_combined(:,:,:,1:8)).^2,4));

[wwater wfat wfreq R2s iter model fitting_error] = fit_IDEAL_R2(iField(:,:,:,:), TE, dfat);

%%
figure(); imagesc(angle(iField(:,:,30,6)));
figure(); imagesc(abs(iField(:,:,30,6)))

%%
figure;
subplot(2,2,1);
imagesc(abs(fat(:,:,36)));colorbar;colormap jet;caxis([0,1]);
subplot(2,2,2);
imagesc(abs(water(:,:,36)));colorbar;colormap jet;caxis([0,1]);
subplot(2,2,3);
imagesc(abs(iFreq(:,:,36)));colorbar;colormap jet;caxis([-pi pi])
subplot(2,2,4);
imagesc(unwph_uf(:,:,36));colorbar;colormap jet;
%%
figure();
imagesc(abs(fat(:,:,36)) ./ (abs(fat(:,:,36)) + abs(water(:,:,36))));
colorbar;colormap jet;caxis([0,1]);
%%
figure();
imagesc(abs(wfat(:,:,36)) ./ (abs(wfat(:,:,36)) + abs(wwater(:,:,36))));
colorbar;colormap jet;caxis([0,1]);

%% CREAM
addpath(genpath(pwd));
rmpath(genpath('./Example/'));
addpath(genpath('../reconstructionPipeline/'));
addpath('../function/');

%%
algoParams=ReadYaml('../CREAM_PDFF-main/algoParams/Hernando_algoParams.yml');

modelParams = struct;
[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

imDataParams = struct;

clear img_grappa
for echo = 1:size(kspace,5)
    img_grappa(:,:,:,:,echo) = fftshift(fftshift(fftshift(squeeze(ifft(ifftshift(ifft(ifftshift(ifft(ifftshift(squeeze(kspace(:,:,:,:,echo)),2),[],2),3),[],3),4),[],4)),4),3),2);
end
imDataParams.images = permute(img_grappa, [2 3 4 1 5]);
imDataParams.TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
imDataParams.PrecessionIsClockwise = 1;
imDataParams.FieldStrength = 2.895;
imDataParams.voxelSize =[2.1 2.1 2.1];
%%
[params, sse] = hernando_main(imDataParams,algoParams);
% save('hernando_simu_SNR50.mat','params','sse'); 