%%

baseDir = 'F:\PhD files backups\Postdoc project\Dixon for James\Data_24P094C2';
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
addpath(genpath('F:\PhD files backups\BOX\PhD projects\QSM\HDR-QSM\Patient study code\supporting'));
TE_o = [1.17e-3, 2.34e-3, 3.51e-3, 4.68e-3, 5.85e-3, 7.02e-3];
TE_o = TE_o - TE_o(1);
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

figure;montage(angle(iField_combined(:,:,end/2,:)./iField_combined(:,:,end/2,1)));colormap jet;clim([-pi,pi]);
figure;montage(angle(iField_combined_norm(:,:,end/2,:)./iField_combined_norm(:,:,end/2,1)));colormap jet;clim([-pi,pi]);
%%
figure(101);
for slc = 1:size(iField_combined_norm, 4)
    subplot(2,3,slc)
    imagesc(angle(iField_combined_norm(:,:,30,slc))); 
end

figure(102);
for slc = 1:size(iField_combined, 4)
    subplot(2,3,slc)
    imagesc(angle(iField_combined(:,:,30,slc)));
end

%%
TE = [1.24, 2.46, 3.69, 4.32, 6.15, 7.38]*1e-3;
TE = TE - TE(1);
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField_combined(:,:,:,1:3),TE(1:3));

iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
figure(); imagesc(iFreq_raw(:,:,64))
%%
%% Draw Mask
iMag = sqrt(sum(abs(iField_combined(:,:,:,1:end)).^2,4));
iMag = iMag./max(iMag(:));
% Mask_o= logical(ROI3D_LRT(iMag./.5));
% Mask = Mask_o(:,:,1:end-1);

% iMag_mask = sqrt(sum(abs(iField_combined(:,:,:,1)).^2,4));
% iMag_mask = iMag_mask./max(iMag_mask(:));
Mask_auto = autoMask(iMag, voxel_size);
Mask = Mask_auto;

Mask([1,end],[1,end],[1,end]) = 0;
% %% SPURS
% f_central = 42.577*1e6*3;
% tic  
%     [unwph_uf,unwph,N_std,iFreq_raw,Mask,wwater,wfat,wfreq,R2s, mask_fat,fat_rad] = ...
%         QG_spursfat_IDEAL_UNIC_RY_Chris(conj(iField_combined(:,:,:,1:end)),TE(1:end)*1e6,f_central,voxel_size,iMag, Mask); 
% t1 = toc;
% N_std_short = N_std;
% fat_map_spurs = abs(wfat)./(abs(wfat)+ abs(wwater));
% water_map_spurs = abs(wwater)./(abs(wfat)+ abs(wwater));
% % %% R2* standard fit
% % [R2s_fitted, rsquare_map, adjrsquare_map] = r2star_fitting(iField_combined_norm, TE);
% %% IDEAL fitting
% delta_TE = (TE(2) - TE(1));
% iField = iField_combined;
% [xx, yy, zz, necho] = size(iField);
% dfat = -3.4e-6*f_central;%sum(fat_shift_specturm.*f_central.*fat_weight);
% R2s_init = 25*ones([1 xx*yy*zz]);
% 
% [wwater_ideal, wfat_ideal, wfreq_ideal, R2s_ideal,~,model_ideal, fitting_error_matrix_ideal, fitting_error_ideal] = ...
%     fit_IDEAL_R2_Chris(conj(iField(:,:,:,1:end)), (TE(1:end)).*1e-3, dfat, (unwph_uf)/(2*pi*delta_TE), R2s_init, 1e1, Mask, mask_fat);
% wfreq_ideal = (wfreq_ideal)*2*pi*delta_TE;
% 
% fat_map_ideal = abs(wfat_ideal)./(abs(wfat_ideal)+ abs(wwater_ideal));
% water_map_ideal = abs(wwater_ideal)./(abs(wfat_ideal)+ abs(wwater_ideal));
% 
% params_ideal.B0 = wfreq_ideal./2./pi./delta_TE;
% params_ideal.R2 = R2s_ideal;
% params_ideal.F  = wfat_ideal;
% params_ideal.W  = wwater_ideal;
% params_ideal.FF  = fat_map_ideal*1e2;
%% Hernando-GC
addpath('C:\Users\yhhua\Box\HDR-QSM\Patient study code\supporting\fat_water\CREAM_PDFF-main');
algoParams_GC=ReadYaml('C:\Users\yhhua\Box\PhD projects\Fat and water seperation\CREAM_PDFF-main\algoParams\Hernando_algoParams.yml');
modelParams_GC =ReadYaml('C:\Users\yhhua\Box\PhD projects\Fat and water seperation\CREAM_PDFF-main\modelParams\CustommodelParams.yml');

[species_GC ,FWspectrum_GC]= setupModelParams(modelParams_GC);
algoParams_GC.gyro=FWspectrum_GC.gyro;
algoParams_GC.species =species_GC;

imDataParams_GC.images = reshape(iField_combined_norm(:,:,:,1:end),[xx, yy, zz,1,necho]);
[nx,ny, nz, ncoils, nTE] = size(imDataParams_GC.images);

% imDataParams.TE = TE(1:end)*1e-3;
imDataParams_GC.TE = (TE(1:end));
imDataParams_GC.FieldStrength = f_central./algoParams_GC.gyro/1e6;
imDataParams_GC.PrecessionIsClockwise = 1;
imDataParams_GC.mask = Mask;
algoParams_GC.range_r2star = [0,300];

[params_GC, sse_GC] = hernando_main(imDataParams_GC,algoParams_GC);
%% IDEAL-CE
algoParams_idealce=ReadYaml('C:\Users\yhhua\Box\HDR-QSM\Patient study code\supporting\fat_water\CREAM_PDFF-main\algoParams\Bydder_algoParams.yml');
modelParams_idealce =ReadYaml('C:\Users\yhhua\Box\HDR-QSM\Patient study code\supporting\fat_water\CREAM_PDFF-main\modelParams\Hodson2008modelParams.yml');

[species_idealce ,FWspectrum_idealce]= setupModelParams(modelParams_idealce);
algoParams_idealce.gyro=FWspectrum_idealce.gyro;
algoParams_idealce.species = species_idealce;
% algoParams_idealce.species(1).name = 'water';
% algoParams_idealce.species(1).frequency = 0;
% algoParams_idealce.species(1).relAmps = 1;
% algoParams_idealce.species(2).name = 'fat';
% algoParams_idealce.species(2).frequency = -3.4;%-3.5;%-3.27;%-dfat_min_shift;%-2.9276;% -3.0138;%[-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
% algoParams_idealce.species(2).relAmps = 1;%[0.087 0.693 0.128 0.004 0.039 0.048];
% algoParams_idealce.nonnegR2 = 201;
% algoParams_idealce.psi = (((-phase_hdr_fit_uf_cen_shift)/(delta_TE*2*pi))) + 1i*R2s_fitted;%(((wfreq)/(2*pi*delta_TE))) + 1i*R2s_o./(2*pi);
algoParams_idealce.muB = 1e-5;%1e-1;
algoParams_idealce.muR = 1e-5;% 5e-2;
algoParams_idealce.filter = ones(3);
algoParams_idealce.smooth_phase = 0;
algoParams_idealce.smooth_field = 1;
algoParams_idealce.maxit = [10, 10, 5];
imDataParams.images = reshape(iField_combined_norm(:,:,:,1:end).*Mask,[xx, yy, zz,1, necho]);
[nx,ny, nz, ncoils, nTE] = size(imDataParams.images);
imDataParams.TE = (TE(1:end));
% imDataParams.TE = (TE(1:end) - TE(1));
imDataParams.FieldStrength = f_central./algoParams_idealce.gyro/1e6;
imDataParams.PrecessionIsClockwise = 1;
imDataParams.mask = Mask;

[params_idealce, sse_idealce, sse1_idealce, YM_idealce] = Bydder_main(imDataParams,algoParams_idealce);