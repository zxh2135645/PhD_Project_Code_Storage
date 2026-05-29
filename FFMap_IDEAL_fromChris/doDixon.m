
%%

% f_glob = glob(cat(2,'/Users/jameszhang/Documents/Data/Echo_*/xspace-532/Slice_', '64', '.mat'));
% 
% data_cell = cell(length(f_glob), 1);
% for i = 1:length(f_glob)
% 
%     data_cell{i} = load(f_glob{i}, 'data');
% 
%     if i == 1
%         S = zeros(size(data_cell{1}.data,1), size(data_cell{1}.data,2),  size(data_cell{1}.data,3), length(f_glob));
%     end
%     S(:,:,:,i) = data_cell{i}.data;
% end

baseDir = '/Users/jameszhang/Documents/Data';
slice_list = 1:104;

data_cell = cell(numel(slice_list),1);
S = [];
expectedNEcho = [];

for k = 1:numel(slice_list)
    s = slice_list(k);
    files = dir(fullfile(baseDir, 'Echo_*', 'xspace-532', sprintf('Slice_%d.mat', s)));
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

% data1 = squeeze(sqrt(sum(data_cell{1}.data .* data_cell{1}.data,1)));
% data2 = squeeze(sqrt(sum(data_cell{2}.data .* data_cell{2}.data,1)));
% figure(); imagesc(abs((data1+data2)/2));
% figure(); imagesc(abs((data1-data2)/2));

%%
%%
iField = permute(S, [2 3 5 4 1]); % RO x PE x Par x Eco x Ch
% combine multiple coils together, assuming the coil is the fifth dimension
iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
iField = sqrt(abs(iField)).*exp(1i*angle(iField));

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
%TE_o = TE_o - TE_o(1);
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
TE = TE_o;
iField_combined = conj(iField_combined);
% iField_combined = iField_combined_norm;
iField_combined_norm = conj(iField_combined_norm);



% %%
% %if size(iField,5)>1
%     % combine multiple coils together, assuming the coil is the fifth dimension
%     % iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
%     mag_iField = abs(iField);
%     mag_iField_norm = mag_iField ./ max(mag_iField(:));
% 
%     Fieldmap_eddy  = sqrt(sqrt((iField(:,:,2)./iField(:,:,1))./(iField(:,:,3)./iField(:,:,2))));
%     iField_uneddy = zeros(size(iField));
%     iField_uneddy(:,:,1:2:end) = iField(:,:,1:2:end).*Fieldmap_eddy;
%     iField_uneddy(:,:,2:2:end) = iField(:,:,2:2:end)./Fieldmap_eddy;
%     iField = mag_iField_norm.*exp(1i*angle(iField_uneddy));
% %end

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

TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField_combined_norm(:,:,:,1:3),TE(1:3));

iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
figure(); imagesc(iFreq_raw(:,:,64))
%%
addpath(genpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/function/Chris/'));

voxel_size = [2.1 2.1 2.1];
f_central = 123.2*10e6;
SUBSAMPLE = 1;
TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
TE = TE - TE(1);
dfat = -3.47e-6*f_central;
LABEL = 2;
L_Fat = 1;
[wunwph_uf ,unwphw ,N_std ,iFreq_raw,iFreq_ref,Mask1,CenterF_ref] = spurs_gc_UNIC_RY_Chris_WorkFlow(L_Fat,iField_combined_norm,TE,f_central,voxel_size,SUBSAMPLE,dfat)

%
% img_grappa = zeros(size(kspace,1), size(kspace,2), size(kspace,3), size(kspace,4), size(kspace,5));
% for echo = 1:size(kspace,5)
%     img_grappa(:,:,:,:,echo) = fftshift(fftshift(fftshift(ifft(ifftshift(ifft(ifftshift(ifft(ifftshift(squeeze(kspace(:,:,:,:,echo)),2),[],2),3),[],3),4),[],4),4),3),2);
% end
% 
% iField = permute(img_grappa, [2 3 4 5 1])*100; % NumRO, NumPE, NumSlices, NumCh, NumEcho =>  NumRO, NumPE, NumSlices, NumEcho, NumCh


% %%
% iField = sqrt(sum(iField.* iField,5));
%%
if size(iField,5)>1
        % combine multiple coils together, assuming the coil is the fifth dimension
        iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
        iField = sqrt(abs(iField)).*exp(1i*angle(iField));
end

% %%
% figure();
% for slc = 1:size(iField, 4)
%     subplot(2,3,slc)
%     imagesc(abs(iField(:,:,44,slc)));
% end
% 
% figure();
% for slc = 1:size(iField, 4)
%     subplot(2,3,slc)
%     imagesc(angle(iField(:,:,44,slc))); 
% end
%%
% Complex data at 3 echoes (same geometry), e.g. [Ny Nx 3]:
% S(:,:,k) must include both magnitude and phase (no abs()!)
% TEs in seconds:
TEs = [1.17e-3, 2.34e-3, 3.51e-3, 4.68e-3, 5.85e-3, 7.02e-3];
TEs = TEs - TEs(1);
% Run (use default ~428 Hz fat-water shift for 3T):
S0 = squeeze(iField_combined_norm(:,:,64,:));
[W,F,fB0,R2s] = dixon3point(S0, TEs, 'DeltaFHz', 428, 'RefineB0', true, 'DoR2s', false);

% Display (magnitude images):
figure; subplot(1,3,1); imshow(abs(W),[]); title('Water');
subplot(1,3,2); imshow(abs(F),[]); title('Fat');
subplot(1,3,3); imshow(fB0,[]);    title('B0 [Hz]');

figure(); imagesc(abs(F)./(abs(F)+abs(W))); caxis([0 1]);colormap gray;
figure(); imagesc(abs(W)./(abs(F)+abs(W))); caxis([0 1]);colormap gray;


%%
rmpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/FFMap_IDEAL_fromChris/Example/MEDIfunctions/_spurs_gc/')
addpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/FFMap_IDEAL_fromChris/MEDIfunctions/_spurs_gc/');
voxel_size = [2.1 2.1 2.1];
f_central = 123.2*10e6;
SUBSAMPLE = 1;
TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
dfat = -3.47e-6*f_central;
LABEL = 2;

Mask_pre = ones((size(iField,1)-2), (size(iField, 2)-2), 1);
Mask = zeros((size(iField,1)), (size(iField, 2)), 1);
Mask(2:end-1, 2:end-1, :) = Mask_pre;

sizes = size(iField);

[water,fat,iFreq,unwph_uf,unwph,N_std,R2s] = ...
    spurs_gc(reshape(iField(:,:,:), sizes(1), sizes(2), 1, sizes(3)),TE,f_central,voxel_size, Mask>0,SUBSAMPLE,dfat,LABEL);

figure;
subplot(2,2,1);
imagesc(abs(fat(:,:)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,2);
imagesc(abs(water(:,:)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,3);
imagesc(abs(iFreq(:,:)));colorbar;colormap jet;caxis([-pi pi])
subplot(2,2,4);
imagesc(unwph_uf(:,:));colorbar;colormap jet;

figure();
imagesc(abs(fat(:,:)) ./ (abs(fat(:,:)) + abs(water(:,:))));
colorbar;colormap jet;caxis([0,1]);

%%

TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField(:,:,:,1:3),TE(1:3));

iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
figure(); imagesc(iFreq_raw(:,:,44))

%% Phase unwrapping

unwrap_map1 = squeeze(tv_pc_silent(iField(:,:,:,:),0.1,5,1));

figure();
for slc = 1:size(iField, 4)
    subplot(2,3,slc)
    imagesc(angle(unwrap_map1(:,:,44,slc)));
end
%%
[iFreq_raw1 N_std] = Fit_ppm_complex_TE(unwrap_map1(:,:,:,1:3),TE(1:3));
iFreq_raw1(isnan(iFreq_raw1))=0;
iFreq_raw1(isinf(iFreq_raw1))=0;
figure(); imagesc(iFreq_raw1(:,:,44))
%%
addpath(genpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/FFMap_IDEAL_fromChris/CREAM_PDFF-main'));
algoName = 'Berglund';
B0_wrapped = iField;
unwrap_map2 = unwrap_B0_algorithms(TE,iFreq_raw,algoName)

%%
[xx yy zz] = size(wunwph_uf);
R2s = zeros([1 xx*yy*zz]);
[wwater wfat wfreq R2s] = fit_IDEAL_R2((iField(:,:,:,:)), TE, dfat, (wunwph_uf-2*pi)/(2*pi*delta_TE),R2s,30);
