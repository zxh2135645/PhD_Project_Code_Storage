%%
clear all;
close all;
baseDir = '/Users/jameszhang/Documents/Data';
slice_list = 1:64;

data_cell = cell(numel(slice_list),1);
S = [];
expectedNEcho = [];

for k = 1:numel(slice_list)
    s = slice_list(k);
    files = dir(fullfile(baseDir, 'Echo_*', 'xspace-532', sprintf('Slice_%d.mat', s)));
    %files = dir(fullfile(baseDir, 'Echo_*', 'xspace-797', sprintf('Slice_%d.mat', s)));
    if isempty(files)
        warning('No file found for Slice %d', s);
        continue;
    end
    nEcho = length(files);

    for eco = 1:length(files)
        fp = fullfile(files(eco).folder, files(eco).name);
        data_struct = load(fp, 'data');
        data_cell{k} = data_struct;
        d = data_struct.data;
        dsz = size(d);
        if numel(dsz) < 3
            error('Loaded data for slice %d has fewer than 3 dimensions', s);
        end
        nx = dsz(1); ny = dsz(2); nz = dsz(3);
        % if numel(dsz) >= 4
        %     nEcho = dsz(4);
        % else
        %     nEcho = 1;
        % end
        if isempty(S)
            % Create S as 5D: [Nx, Ny, Nz, NEcho, NSlices]
            S = zeros(nx, ny, nz, nEcho, numel(slice_list));
            expectedNEcho = nEcho;
        else
            if nEcho ~= expectedNEcho
                error('Inconsistent number of echoes across slices: expected %d, got %d (slice %d)', expectedNEcho, nEcho, s);
            end
        end
        if nEcho == 1
            S(:,:,:,1,k) = d;
        else
            S(:,:,:,eco,k) = d;
        end
    end
end

%%
iField = permute(S, [2 3 5 4 1]); % RO x PE x Par x Eco x Ch
% combine multiple coils together, assuming the coil is the fifth dimension
iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
iField = sqrt(abs(iField)).*exp(1i*angle(iField));

%%
iField = permute(S, [2 3 5 4 1]); % RO x PE x Par x Eco x Ch
iField = sqrt(sum(iField.*iField,5));  %

%%
figure();
for slc = 1:size(iField, 4)
    subplot(2,3,slc)
    imagesc(abs(iField(:,:,64,slc)));
end

figure();
for slc = 1:size(iField, 4)
    subplot(2,3,slc)
    imagesc(angle(iField(:,:,64,slc)));
end

%% Bi-polar correction
addpath(genpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/function/Chris/'));
TE_o = [1.17e-3, 2.34e-3, 3.51e-3, 4.68e-3, 5.85e-3, 7.02e-3];

TE_o = TE_o - TE_o(1);
TE_o = [1.17e-3, 2.34e-3, 3.51e-3, 4.68e-3, 5.85e-3, 7.02e-3];

voxel_size = [2.1 2.1 2.1];

%%%%%%  Creat B0map Normalization and correct for phase drift
iField_combined_angle_norm = angle(iField.*conj(iField(:,:,:,1)));
iField_combined_mag = abs(iField);
%%%% Bipolar eddy current correction %%%%
iField_drift_corrected_odd = phase_drift_correction(iField(:,:,:,1:2:end), TE_o(1:2:end), voxel_size);
iField_drift_corrected_even = phase_drift_correction(iField(:,:,:,2:2:end), TE_o(2:2:end), voxel_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iField_combined_angle = ones(size(iField));
iField_combined_angle(:,:,:,1:2:end) = angle(iField_drift_corrected_odd);
iField_combined_angle(:,:,:,2:2:end) = angle(iField_drift_corrected_even);
iField_combined_angle(isnan(iField_combined_angle)) = 0;
iField_combined_mag_vec_temp = iField_combined_mag(:);
iField_combined_mag_norm = iField_combined_mag./max(iField_combined_mag_vec_temp(:));
iField_combined = iField_combined_mag_norm.*exp(1i*iField_combined_angle);
iField_combined_norm = iField_combined_mag_norm.*exp(1i*iField_combined_angle_norm);
%TE = TE_o;
iField_combined = conj(iField_combined);
% iField_combined = iField_combined_norm;
iField_combined_norm = conj(iField_combined_norm);

%%
figure();
for slc = 1:size(iField_combined_norm, 4)
    subplot(2,3,slc)
    imagesc(abs(iField_combined_norm(:,:,64,slc)));
end

figure();
for slc = 1:size(iField_combined_norm, 4)
    subplot(2,3,slc)
    imagesc(angle(iField_combined_norm(:,:,64,slc)));
end

%%
TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
TE = TE - TE(1);
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField_combined_norm(:,:,:,1:3),TE(1:3));
%TE = [1.24, 3.69, 6.15]*1e-3;
%[iFreq_raw N_std] = Fit_ppm_complex_TE(iField(:,:,:,1:3),TE(1:3));

iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
figure(); imagesc(iFreq_raw(:,:,64))

delta_TE = (TE(2) - TE(1));

voxel_size = [2.1 2.1 2.1];
f_central = 123.2*10e6;
dfat = -3.47e-6*f_central;

dyna_range = 1/delta_TE;
effect_fat_Hz = dfat + floor( (0.5*dyna_range-dfat)/dyna_range)*dyna_range;
effect_fat_rad = effect_fat_Hz/dyna_range*2*pi;


p = 2;
w1 = effect_fat_rad/pi;
z = size(iFreq_raw,3);
w1
%%
w = -w1;
%iMag = sqrt(sum(abs(iField).^2,4));
iMag = sqrt(sum(abs(iField_combined_norm).^2,4));
Mask=autoMask(iMag,voxel_size);
[unwphw,iter,erglist] = phase_unwrap_3d_UNIC_Chris_1002(iFreq_raw,p,iMag,voxel_size,Mask);   
%%
figure();
imagesc(abs(unwphw(:,:,64)));
%%
%phase_3d = angle(iField(:,:,:,2)./iField(:,:,:,1));

phase_3d = angle(iField_combined_norm(:,:,:,2)./iField_combined_norm(:,:,:,1));
[xx yy zz] = size(phase_3d);
R2s = zeros([1 xx*yy*zz]);
[wwater wfat wfreq R2s iter model] = fit_IDEAL_R2(iField_combined_norm(:,:,:,:), TE, dfat, (unwphw-2*pi)./(2*pi*delta_TE), R2s, 30);
%[wwater wfat wfreq R2s iter model] = fit_IDEAL_R2(iField(:,:,:,[1 3 5]), TE, dfat, (unwphw-2*pi)./(2*pi*delta_TE), R2s, 30);


figure;
subplot(2,2,1);
imagesc(abs(wfat(:,:,64)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,2);
imagesc(abs(wwater(:,:,64)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,3);
imagesc(abs(wfreq(:,:,64)));colorbar;colormap jet;caxis([-pi pi])
subplot(2,2,4);
imagesc(unwphw(:,:,64));colorbar;colormap jet;

figure();
imagesc(abs(wfat(:,:)) ./ (abs(wfat(:,:)) + abs(wwater(:,:))));
colorbar;colormap jet;caxis([0,1]);
