clear all;
close all;
clc;
%% Need output from exvivo_invivo_register
addpath('../function/');
base_dir = uigetdir;
name = '18D16';
name = 'Sahara';
name = 'Mojave';
time_point = 'EXVIVO';

reg_info_f = cat(2, base_dir, '/data/', name, '/RegInfo.mat');
load(reg_info_f);

dicom_dir = uigetdir; % Exvivo DICOM of T2* weighted image
folder_glob = glob(cat(2, dicom_dir, '\*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);

strings = strsplit(list_to_read{order_to_read}, '/');
f_name = strings{end-1};

strings = strsplit(dicom_dir, '/');
name_folder = strings{end};
series_desc = f_name(1:end-5);
% ff_dir = cat(2, base_dir, '/FF_Data/', name, '/', name, '_', time_point, '_', f_name, '.mat');
ff_dir = GetFullPath(cat(2, dicom_dir, '/../FF_Data_06232022/', name_folder, '/', series_desc, '.mat'));
load(ff_dir, 'fwmc_ff', 'fwmc_r2star');

ff_map = fwmc_ff;
r2star_map = fwmc_r2star;

exvivo_match_flip = flip(exvivo_match);
ff_map_match = ff_map(:,:,exvivo_match_flip);
r2star_map_match = r2star_map(:,:,exvivo_match_flip);

ff_hb_dir = GetFullPath(cat(2, dicom_dir, '/../Result_FromWorkstation_07052022/AllPhasemap_', upper(name), '.mat'));
load(ff_hb_dir, 'fat_r2s', 'water_r2s', 'R2s');
ff_map_hb = abs(fat_r2s) ./ (abs(fat_r2s) + abs(water_r2s));
ff_map_match_hb = ff_map_hb(:,:,exvivo_match_flip);
R2s_match_hb = R2s(:,:,exvivo_match_flip);

%% Load Exvivo Masks
addpath('../AHA16Segment/');
strings = strsplit(list_to_read{1}, '/');
strings = strsplit(strings{end-1}, '_');
num_label = strings{end};
num_label = num2str(str2num(num_label) + 2, '%.4d');

mask_f_lv = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_LV.mat');
load(mask_f_lv);

mask_f_blood = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_Blood.mat');
load(mask_f_x = lsqnonlin(fun,x0,lb,ub)blood);

mask_f_insert = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_InsertionPt.mat');
load(mask_f_insert);

se2 = strel('disk',2);
mask_lv_erode = imerode(mask_lv, se2);
mask_myocardium = double(mask_lv_erode) - mask_blood;
mask_myocardium_match = mask_myocardium(:,:,exvivo_match_flip);

mask_myocardium_match_nan = mask_myocardium_match;
mask_myocardium_match_nan(mask_myocardium_match_nan == 0) = nan;

figure();
for i = 1:size(ff_map_match, 3)
    subplot(2,5,i);
   imagesc(ff_map_match(:,:,i) .* mask_myocardium_match_nan(:,:,i)); axis image; caxis([-2 20]);
   title(cat(2, 'Slice ', num2str(i)));
   subplot(2,5,i+5);
   imagesc(transpose(ff_map_match_hb(:,:,i)*100) .* mask_myocardium_match_nan(:,:,i)); axis image; caxis([-2 20]);
   title(cat(2, 'Slice ', num2str(i)));
end

figure();
for i = 1:size(r2star_map_match, 3)
    subplot(2,5,i);
   imagesc(r2star_map_match(:,:,i) .* mask_myocardium_match_nan(:,:,i)); axis image; caxis([-10 100]);
   title(cat(2, 'Slice ', num2str(i)));
   subplot(2,5,i+5);
   imagesc(transpose(R2s_match_hb(:,:,i)) .* mask_myocardium_match_nan(:,:,i)); axis image; caxis([-10 100]);
   title(cat(2, 'Slice ', num2str(i)));
end

ff_map_match(ff_map_match > 100) = 100;
ff_map_match(ff_map_match < 0) = 0;
r2star_map_match(r2star_map_match < 0) = 0;
r2star_map_match(r2star_map_match > 200) = 200;

centroids_exvivo_match = centroids_exvivo(exvivo_match);
thera_i_array_exvivo_match = -thera_i_array_exvivo(exvivo_match) * 180 / pi;
thera_s_array_exvivo_match = -thera_s_array_exvivo(exvivo_match) * 180 / pi;

se = strel('disk', 1);
Segn = 50;

ff_mask_3d = zeros(size(mask_myocardium_match));
r2star_mask_3d = zeros(size(mask_myocardium_match));
ff_mask_3d_hb = zeros(size(mask_myocardium_match));
r2star_mask_3d_hb = zeros(size(mask_myocardium_match));

chords = zeros(Segn, 2, size(mask_myocardium_match, 3));
chords_r2star = zeros(Segn, 2, size(mask_myocardium_match, 3));
chords_hb = zeros(Segn, 2, size(mask_myocardium_match, 3));
chords_r2star_hb = zeros(Segn, 2, size(mask_myocardium_match, 3));

for slc = 1:size(mask_myocardium_match, 3)
    BW_skel = bwmorph(mask_myocardium_match(:,:,slc), 'skel', Inf);
    center_fixed = imfill(BW_skel, 'hole');
    center_fixed = imopen(center_fixed, se); % Removing spikes
    myo_epi = mask_myocardium_match(:,:,slc) - center_fixed > 0;
    myo_endo = center_fixed + mask_myocardium_match(:,:,slc) > 1;
    % Groove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
    Groove = thera_s_array_exvivo_match(slc); % Only true for the current 18D16
    
    [Segmentpix_endo, stats, Mask_Segn_endo] = AHASegmentation(ff_map_match(:,:,slc), myo_endo, Segn, Groove);
    [Segmentpix_epi, stats, Mask_Segn_epi] = AHASegmentation(ff_map_match(:,:,slc), myo_epi, Segn, Groove);
    
    [Segmentpix_endo_r2star, stats, Mask_Segn_endo_r2star] = AHASegmentation(r2star_map_match(:,:,slc), myo_endo, Segn, Groove);
    [Segmentpix_epi_r2star, stats, Mask_Segn_epi_r2star] = AHASegmentation(r2star_map_match(:,:,slc), myo_epi, Segn, Groove);
    
    [Segmentpix_endo_hb, stats, Mask_Segn_endo_hb] = AHASegmentation(transpose(ff_map_match_hb(:,:,slc)), myo_endo, Segn, Groove);
    [Segmentpix_epi_hb, stats, Mask_Segn_epi_hb] = AHASegmentation(transpose(ff_map_match_hb(:,:,slc)), myo_epi, Segn, Groove);
    
    [Segmentpix_endo_r2star_hb, stats, Mask_Segn_endo_r2star_hb] = AHASegmentation(transpose(R2s_match_hb(:,:,slc)), myo_endo, Segn, Groove);
    [Segmentpix_epi_r2star_hb, stats, Mask_Segn_epi_r2star_hb] = AHASegmentation(transpose(R2s_match_hb(:,:,slc)), myo_epi, Segn, Groove);

    ff_mask = zeros(size(Mask_Segn_endo));
    r2star_mask = zeros(size(Mask_Segn_endo));
    ff_mask_hb = zeros(size(Mask_Segn_endo));
    r2star_mask_hb = zeros(size(Mask_Segn_endo));
    for i = 1:Segn
        seg_mean_endo = mean(Segmentpix_endo{i});
        seg_mean_epi = mean(Segmentpix_epi{i});
        chords(i, 1, slc) = seg_mean_endo;
        chords(i, 2, slc) = seg_mean_epi;
        
        seg_mean_endo_r2star = mean(Segmentpix_endo_r2star{i});
        seg_mean_epi_r2star = mean(Segmentpix_epi_r2star{i});
        chords_r2star(i, 1, slc) = seg_mean_endo_r2star;
        chords_r2star(i, 2, slc) = seg_mean_epi_r2star;
        
        ff_mask(Mask_Segn_endo == i) = seg_mean_endo;
        ff_mask(Mask_Segn_epi == i) = seg_mean_epi;
        r2star_mask(Mask_Segn_endo == i) = seg_mean_endo_r2star;
        r2star_mask(Mask_Segn_epi == i) = seg_mean_epi_r2star;

        seg_mean_endo_hb = mean(Segmentpix_endo_hb{i});
        seg_mean_epi_hb = mean(Segmentpix_epi_hb{i});
        chords_hb(i, 1, slc) = seg_mean_endo_hb;
        chords_hb(i, 2, slc) = seg_mean_epi_hb;
        
        seg_mean_endo_r2star_hb = mean(Segmentpix_endo_r2star_hb{i});
        seg_mean_epi_r2star_hb = mean(Segmentpix_epi_r2star_hb{i});
        chords_r2star_hb(i, 1, slc) = seg_mean_endo_r2star_hb;
        chords_r2star_hb(i, 2, slc) = seg_mean_epi_r2star_hb;
        
        ff_mask_hb(Mask_Segn_endo_hb == i) = seg_mean_endo_hb;
        ff_mask_hb(Mask_Segn_epi_hb == i) = seg_mean_epi_hb;
        r2star_mask_hb(Mask_Segn_endo_hb == i) = seg_mean_endo_r2star_hb;
        r2star_mask_hb(Mask_Segn_epi_hb == i) = seg_mean_epi_r2star_hb;
    end
    
    ff_mask_3d(:,:,slc) = ff_mask;
    r2star_mask_3d(:,:,slc) = r2star_mask;
    ff_mask_3d_hb(:,:,slc) = ff_mask_hb;
    r2star_mask_3d_hb(:,:,slc) = r2star_mask_hb;
end

LABEL = 1;
if LABEL == 1
    chords = flip(chords, 1);
    chords_r2star = flip(chords_r2star, 1);
    chords_hb = flip(chords_hb, 1);
    chords_r2star_hb = flip(chords_r2star_hb, 1);
end
%chords = flip(chords, 1);
%chords_r2star = flip(chords_r2star, 1);

figure();
subplot(1,2,1);
imagesc(ff_mask_3d(:,:,3)); caxis([0 10]);
subplot(1,2,2);
imagesc(ff_mask_3d_hb(:,:,3)*100); caxis([0 10]);

figure();
subplot(1,2,1);
imagesc(r2star_mask_3d(:,:,3)); caxis([0 100]);
subplot(1,2,2);
imagesc(r2star_mask_3d_hb(:,:,3)); caxis([0 100]);
%% Compare to invivo FF map
%time_point = '9MO';
name = 'Mojave';
time_point = '1YR';
sequence_label = {'T1', 'T2star', 'LGE', 'T2'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
label_t2 = sequence_label{4};

tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
if ~exist(tp_dir, 'dir')
    disp(cat(2, 'No folder at: ', name, ' ', time_point));
else
    % T1
    % tp_count = tp_count+1;
    myo_glob = glob(cat(2, tp_dir, label_t1, '/', anatomy_label{5}, '/*'));
    roi_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{3}, '/*'));
    remote_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{6}, '/*'));
    noreflow_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{7}, '/*'));
    
    load(cat(2, tp_dir, label_t1, '/', label_t1, '_vol_img_3D.mat'));
    load(myo_glob{1});
    load(roi_glob{1});
    load(remote_glob{1});
    load(noreflow_glob{1});
    load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
    [slc_array_t1, idx_reordered_t1] = sort(slc_array);
    
    roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
    remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
    hemo_core_t1 = noReflowMask_3D .* mask_myocardium_3D;
    roi_t1 = roi_in_myo_t1 .* vol_img_3D;
    remote_t1 = remote_in_myo_t1 .* vol_img_3D;
    t1 = vol_img_3D;
    myo_t1 = imerode(mask_myocardium_3D,se);
    
    roi_in_myo_t1 = roi_in_myo_t1(:,:,idx_reordered_t1);
    remote_in_myo_t1 = remote_in_myo_t1(:,:,idx_reordered_t1);
    roi_t1 = roi_t1(:,:,idx_reordered_t1);
    remote_t1 = remote_t1(:,:,idx_reordered_t1);
    t1 = t1(:,:,idx_reordered_t1);
    myo_t1 = myo_t1(:,:,idx_reordered_t1);
    hemo_core_t1 = hemo_core_t1(:,:,idx_reordered_t1);
    
    myo_t1_nan = myo_t1;
    myo_t1_nan(myo_t1 == 0) = nan;
    
    % FF
    myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
    roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
    remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
    
    load(myo_glob{1});
    load(roi_glob{1});
    load(remote_glob{1});
    load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat'));
    load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
    slc_array_ff = slc_array;
    
    ff_map = cell(1, length(glob_names));
    for f = 1:length(ff_map)
        ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
    end
    % convert ff_map to matrix
    ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2),size(ff_map{1}.fwmc_ff,3),  length(ff_map));
    
    for f = 1:length(ff_map)
        ff(:,:,:,f) = ff_map{f}.fwmc_ff;
    end
    ff = squeeze(ff);
    ff(ff > 100) = 100;
    ff(ff < 0) = 0;
    
    roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
    remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
    %hemo_core_t2star = mask_myocardium_3D .* freeROIMask_3D_te8;
    myo_ff = imerode(mask_myocardium_3D,se);
    %myo_ff = mask_myocardium_3D;
    
    roi_ff = roi_in_myo_ff .* ff;
    remote_ff = remote_in_myo_ff .* ff;

    idx_reordered_ff = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
    ff = ff(:,:,idx_reordered_ff);
    myo_ff = myo_ff(:,:,idx_reordered_ff);
    remote_ff = remote_ff(:,:,idx_reordered_ff);
    roi_ff = roi_ff(:,:,idx_reordered_ff);
    remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered_ff);
    roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered_ff);
    %hemo_core_t2star = hemo_core_t2star(:,:,idx_reordered);
    myo_ff_nan = myo_ff;
    myo_ff_nan(myo_ff == 0) = nan;

    myo_ff_hb = imtranslate(myo_ff,[-1, -1]);
    myo_ff_hb_nan = myo_ff_hb;
    myo_ff_hb_nan(myo_ff_hb == 0) = nan;

    idx_reordered_inex = Func_AlignSliceLoc(slc_array_ff(idx_reordered_ff), cine_SliceLoc);
    thera_s_array_cine_match = -thera_s_array_cine(idx_reordered_inex) * 180 / pi;
    
    % R2star Map
    r2star_map = cell(1, length(glob_names));
    for f = 1:length(r2star_map)
        r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
    end
    
    % convert ff_map to matrix
    r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{2}.fwmc_r2star, 2), length(r2star_map));
    for f = 1:length(r2star_map)
        r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
    end
    
    r2star = r2star(:,:,idx_reordered_ff);
    r2star(r2star < 0) = 0;
    r2star(r2star > 200) = 200;
    
    % T2
    load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
    myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
    roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
    remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));
    % noreflow_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{7}, '/*'));
    load(myo_glob{1});
    load(roi_glob{1});
    load(remote_glob{1});
    % load(noreflow_glob{1});
    load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat'));
    slc_array_t2 = slc_array;
    
    roi_in_myo_t2 = mask_myocardium_3D .* freeROIMask_3D;
    remote_in_myo_t2 = mask_myocardium_3D .* myoRefMask_3D;
    roi_t2 = roi_in_myo_t2 .* vol_img_3D;
    remote_t2 = remote_in_myo_t2 .* vol_img_3D;
    t2 = vol_img_3D;
    myo_t2 = imerode(mask_myocardium_3D,se);
    %hemo_core_t2 = noReflowMask_3D .* mask_myocardium_3D;
    
    idx_reordered_t2 = Func_AlignSliceLoc(slc_array_t1, slc_array_t2);
    roi_in_myo_t2 = roi_in_myo_t2(:,:,idx_reordered_t2);
    remote_in_myo_t2 = remote_in_myo_t2(:,:,idx_reordered_t2);
    roi_t2 = roi_t2(:,:,idx_reordered_t2);
    remote_t2 = remote_t2(:,:,idx_reordered_t2);
    t2 = t2(:,:,idx_reordered_t2);
    myo_t2 = myo_t2(:,:,idx_reordered_t2);
    %hemo_core_t2 = hemo_core_t2(:,:,idx_reordered);
    
    myo_t2_nan = myo_t2;
    myo_t2_nan(myo_t2 == 0) = nan;
    
    % load homebrew FF and r2star
    ideal_hb = load(GetFullPath(cat(2, dicom_dir, '/../../T1FP_Data/Dharmakumar_London_', name, '_12month_Invivo/Result/AllPhasemap.mat')));
    ff_hb = abs(ideal_hb.fat) ./ (abs(ideal_hb.fat) + abs(ideal_hb.water));
    r2star_hb = ideal_hb.R2s;
    ff_hb = ff_hb(:,:,idx_reordered_ff);
    r2star_hb = r2star_hb(:,:,idx_reordered_ff);
    
    for slc = 1:size(ff, 3)
        figure();
        subplot(1,2,1);
        imagesc(ff(:,:,slc)); caxis([0 20]); axis image;
        subplot(1,2,2);
        imagesc(ff_hb(:,:,slc)*100); caxis([0 20]); axis image;
    end

    for slc = 1:size(ff, 3)
        figure();
        subplot(1,2,1);
        imagesc(r2star(:,:,slc)); caxis([0 100]); axis image;
        subplot(1,2,2);
        imagesc(r2star_hb(:,:,slc)); caxis([0 100]); axis image;
    end

    ff_mask_invivo_3d = zeros(size(ff));
    ff_mask_exvivo_3d = zeros(size(ff));
    ff_mask_invivo_3d_hb = zeros(size(ff));
    ff_mask_exvivo_3d_hb = zeros(size(ff));

    t1_mask_invivo_3d = zeros(size(t1));
    r2star_mask_invivo_3d = zeros(size(ff));
    r2star_mask_exvivo_3d = zeros(size(ff));
    r2star_mask_invivo_3d_hb = zeros(size(ff));
    r2star_mask_exvivo_3d_hb = zeros(size(ff));

    t2_mask_invivo_3d = zeros(size(t2));
    chords_invivo = zeros(Segn, 2, size(ff, 3));
    chords_invivo_t1 = zeros(Segn, 2, size(t1, 3));
    chords_invivo_r2star = zeros(Segn, 2, size(r2star, 3));
    chords_invivo_t2 = zeros(Segn, 2, size(t2, 3));
    mi_label_invivo = zeros(Segn, 2, size(ff, 3));
    
    chords_invivo_hb = zeros(Segn, 2, size(ff, 3));
    chords_invivo_r2star_hb = zeros(Segn, 2, size(r2star, 3));
    
    for slc = 1:size(ff, 3)
        BW_skel = bwmorph(myo_ff(:,:,slc), 'skel', Inf);
        center_fixed = imfill(BW_skel, 'hole');
        center_fixed = imopen(center_fixed, se); % Removing spikes
        myo_epi_ff = myo_ff(:,:,slc) - center_fixed > 0;
        myo_endo_ff = center_fixed + myo_ff(:,:,slc) > 1;
        
        BW_skel = bwmorph(myo_t1(:,:,slc), 'skel', Inf);
        center_fixed = imfill(BW_skel, 'hole');
        center_fixed = imopen(center_fixed, se); % Removing spikes
        myo_epi_t1 = myo_t1(:,:,slc) - center_fixed > 0;
        myo_endo_t1 = center_fixed + myo_t1(:,:,slc) > 1;
        
        BW_skel = bwmorph(myo_t2(:,:,slc), 'skel', Inf);
        center_fixed = imfill(BW_skel, 'hole');
        center_fixed = imopen(center_fixed, se); % Removing spikes
        myo_epi_t2 = myo_t2(:,:,slc) - center_fixed > 0;
        myo_endo_t2 = center_fixed + myo_t2(:,:,slc) > 1;
        
        BW_skel = bwmorph(myo_ff_hb(:,:,slc), 'skel', Inf);
        center_fixed = imfill(BW_skel, 'hole');
        center_fixed = imopen(center_fixed, se); % Removing spikes
        myo_epi_ff_hb = myo_ff_hb(:,:,slc) - center_fixed > 0;
        myo_endo_ff_hb = center_fixed + myo_ff_hb(:,:,slc) > 1;

        % Groove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
        Groove = thera_s_array_cine_match(slc);
        
        [Segmentpix_endo, stats, Mask_Segn_endo] = AHASegmentation(ff(:,:,slc), myo_endo_ff, Segn, Groove);
        [Segmentpix_epi, stats, Mask_Segn_epi] = AHASegmentation(ff(:,:,slc), myo_epi_ff, Segn, Groove);
        
        [Segmentpix_endo_t1, stats, Mask_Segn_endo_t1] = AHASegmentation(t1(:,:,slc), myo_endo_t1, Segn, Groove);
        [Segmentpix_epi_t1, stats, Mask_Segn_epi_t1] = AHASegmentation(t1(:,:,slc), myo_epi_t1, Segn, Groove);
        
        [Segmentpix_endo_r2star, stats, Mask_Segn_endo_r2star] = AHASegmentation(r2star(:,:,slc), myo_endo_ff, Segn, Groove);
        [Segmentpix_epi_r2star, stats, Mask_Segn_epi_r2star] = AHASegmentation(r2star(:,:,slc), myo_epi_ff, Segn, Groove);
        
        [Segmentpix_endo_t2, stats, Mask_Segn_endo_t2] = AHASegmentation(t2(:,:,slc), myo_endo_t2, Segn, Groove);
        [Segmentpix_epi_t2, stats, Mask_Segn_epi_t2] = AHASegmentation(t2(:,:,slc), myo_epi_t2, Segn, Groove);
        
        [Segmentpix_endo_hb, stats, Mask_Segn_endo_hb] = AHASegmentation(ff_hb(:,:,slc), myo_endo_ff_hb, Segn, Groove);
        [Segmentpix_epi_hb, stats, Mask_Segn_epi_hb] = AHASegmentation(ff_hb(:,:,slc), myo_epi_ff_hb, Segn, Groove);

        [Segmentpix_endo_r2star_hb, stats, Mask_Segn_endo_r2star_hb] = AHASegmentation(r2star_hb(:,:,slc), myo_endo_ff_hb, Segn, Groove);
        [Segmentpix_epi_r2star_hb, stats, Mask_Segn_epi_r2star_hb] = AHASegmentation(r2star_hb(:,:,slc), myo_epi_ff_hb, Segn, Groove);
        
        ff_mask_invivo = zeros(size(Mask_Segn_endo));
        ff_mask_exvivo = zeros(size(Mask_Segn_endo));
        t1_mask_invivo = zeros(size(Mask_Segn_endo_t1));
        t2_mask_invivo = zeros(size(Mask_Segn_endo_t2));
        r2star_mask_invivo = zeros(size(Mask_Segn_endo_r2star));
        r2star_mask_exvivo = zeros(size(Mask_Segn_endo_r2star));
        
        ff_mask_invivo_hb = zeros(size(Mask_Segn_endo));
        ff_mask_exvivo_hb = zeros(size(Mask_Segn_endo));
        r2star_mask_invivo_hb = zeros(size(Mask_Segn_endo_r2star));
        r2star_mask_exvivo_hb = zeros(size(Mask_Segn_endo_r2star));

        for i = 1:Segn
            seg_mean_endo = mean(Segmentpix_endo{i});
            seg_mean_epi = mean(Segmentpix_epi{i});
            chords_invivo(i,1,slc) = seg_mean_endo;
            chords_invivo(i,2,slc) = seg_mean_epi;
            
            seg_mean_endo_t1 = mean(Segmentpix_endo_t1{i});
            seg_mean_epi_t1 = mean(Segmentpix_epi_t1{i});
            chords_invivo_t1(i,1,slc) = seg_mean_endo_t1;
            chords_invivo_t1(i,2,slc) = seg_mean_epi_t1;
            
            seg_mean_endo_r2star = mean(Segmentpix_endo_r2star{i});
            seg_mean_epi_r2star = mean(Segmentpix_epi_r2star{i});
            chords_invivo_r2star(i,1,slc) = seg_mean_endo_r2star;
            chords_invivo_r2star(i,2,slc) = seg_mean_epi_r2star;
            
            seg_mean_endo_t2 = mean(Segmentpix_endo_t2{i});
            seg_mean_epi_t2 = mean(Segmentpix_epi_t2{i});
            chords_invivo_t2(i,1,slc) = seg_mean_endo_t2;
            chords_invivo_t2(i,2,slc) = seg_mean_epi_t2;
            
            seg_mean_endo_hb = mean(Segmentpix_endo_hb{i});
            seg_mean_epi_hb = mean(Segmentpix_epi_hb{i});
            chords_invivo_hb(i,1,slc) = seg_mean_endo_hb;
            chords_invivo_hb(i,2,slc) = seg_mean_epi_hb;

            seg_mean_endo_r2star_hb = mean(Segmentpix_endo_r2star_hb{i});
            seg_mean_epi_r2star_hb = mean(Segmentpix_epi_r2star_hb{i});
            chords_invivo_r2star_hb(i,1,slc) = seg_mean_endo_r2star_hb;
            chords_invivo_r2star_hb(i,2,slc) = seg_mean_epi_r2star_hb;

            if any(any(Mask_Segn_endo.*roi_in_myo_ff(:,:,slc) == i))
                mi_label_invivo(i,1,slc) = 1;
            end
            
            if any(any(Mask_Segn_epi.*roi_in_myo_ff(:,:,slc) == i))
                mi_label_invivo(i,2,slc) = 1;
            end
            
            ff_mask_invivo(Mask_Segn_endo == i) = seg_mean_endo;
            ff_mask_invivo(Mask_Segn_epi == i) = seg_mean_epi;
            t1_mask_invivo(Mask_Segn_endo_t1 == i) = seg_mean_endo_t1;
            t1_mask_invivo(Mask_Segn_epi_t1 == i) = seg_mean_epi_t1;
            r2star_mask_invivo(Mask_Segn_endo_r2star == i) = seg_mean_endo_r2star;
            r2star_mask_invivo(Mask_Segn_epi_r2star == i) = seg_mean_epi_r2star;
            t2_mask_invivo(Mask_Segn_endo_t2 == i) = seg_mean_endo_t2;
            t2_mask_invivo(Mask_Segn_epi_t2 == i) = seg_mean_epi_t2;
            
            ff_mask_exvivo(Mask_Segn_endo == i) = chords(i,1,idx_reordered_inex(slc));
            ff_mask_exvivo(Mask_Segn_epi == i) = chords(i,2,idx_reordered_inex(slc));
            r2star_mask_exvivo(Mask_Segn_endo_r2star == i) = chords_r2star(i,1,idx_reordered_inex(slc));
            r2star_mask_exvivo(Mask_Segn_epi_r2star == i) = chords_r2star(i,2,idx_reordered_inex(slc));

            ff_mask_invivo_hb(Mask_Segn_endo_hb == i) = seg_mean_endo_hb;
            ff_mask_invivo_hb(Mask_Segn_epi_hb == i) = seg_mean_epi_hb;
            r2star_mask_invivo_hb(Mask_Segn_endo_r2star_hb == i) = seg_mean_endo_r2star_hb;
            r2star_mask_invivo_hb(Mask_Segn_epi_r2star_hb == i) = seg_mean_epi_r2star_hb;
            ff_mask_exvivo_hb(Mask_Segn_endo_hb == i) = chords_hb(i,1,idx_reordered_inex(slc));
            ff_mask_exvivo_hb(Mask_Segn_epi_hb == i) = chords_hb(i,2,idx_reordered_inex(slc));
            r2star_mask_exvivo_hb(Mask_Segn_endo_r2star_hb == i) = chords_r2star_hb(i,1,idx_reordered_inex(slc));
            r2star_mask_exvivo_hb(Mask_Segn_epi_r2star_hb == i) = chords_r2star_hb(i,2,idx_reordered_inex(slc));
        end
        
        ff_mask_invivo_3d(:,:,slc) = ff_mask_invivo;
        ff_mask_exvivo_3d(:,:,slc) = ff_mask_exvivo;
        t1_mask_invivo_3d(:,:,slc) = t1_mask_invivo;
        t2_mask_invivo_3d(:,:,slc) = t2_mask_invivo;
        r2star_mask_invivo_3d(:,:,slc) = r2star_mask_invivo;
        r2star_mask_exvivo_3d(:,:,slc) = r2star_mask_exvivo;

        ff_mask_invivo_3d_hb(:,:,slc) = ff_mask_invivo_hb;
        ff_mask_exvivo_3d_hb(:,:,slc) = ff_mask_exvivo_hb;
        r2star_mask_invivo_3d_hb(:,:,slc) = r2star_mask_invivo_hb;
        r2star_mask_exvivo_3d_hb(:,:,slc) = r2star_mask_exvivo_hb;
    end
    
    slc_for_plot = 3;
    figure();
    subplot(2,3,1);
    imagesc(ff_mask_invivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([-2 15]); title('Invivo FF'); colorbar;
    subplot(2,3,2);
    imagesc(ff_mask_exvivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([-2 15]); title('Exvivo FF'); colorbar;
    subplot(2,3,4);
    imagesc(r2star_mask_invivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([10 80]); title('Invivo R2star'); colorbar;
    subplot(2,3,5);
    imagesc(r2star_mask_exvivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([10 80]); title('Exvivo R2star'); colorbar;
    subplot(2,3,3);
    imagesc(t1_mask_invivo_3d(:,:,slc_for_plot).* myo_t1_nan(:,:,slc_for_plot)); title('Invivo T1'); colorbar;
    subplot(2,3,6);
    imagesc(t2_mask_invivo_3d(:,:,slc_for_plot)./10.* myo_t2_nan(:,:,slc_for_plot)); title('Invivo T2'); caxis([20 50]); colorbar;
    colormap(brewermap([],'*RdYlBu'));

%     figure();
%     subplot(2,3,1);
%     imagesc(ff(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([-2 10]); title('Invivo FF'); colorbar;
%     subplot(2,3,2);
%     imagesc(ff_map_match(:,:,slc_for_plot).* mask_myocardium_match(:,:,slc_for_plot)); caxis([-2 10]); title('Exvivo FF'); colorbar;
%     subplot(2,3,4);
%     imagesc(r2star(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([10 50]); title('Invivo R2star'); colorbar;
%     subplot(2,3,5);
%     imagesc(r2star_map_match(:,:,slc_for_plot).* mask_myocardium_match(:,:,slc_for_plot)); caxis([10 50]); title('Exvivo R2star'); colorbar;
%     colormap(brewermap([],'*RdYlBu'));
    
    figure();
    subplot(2,4,1);
    imagesc(ff(:,:,slc_for_plot)); caxis([-2 20]); title('Invivo FF'); colorbar;
    subplot(2,4,3);
    imagesc(ff_map_match(:,:,slc_for_plot)); caxis([-2 20]); title('Exvivo FF'); colorbar;
    subplot(2,4,2);
    imagesc(ff_hb(:,:,slc_for_plot)*100); caxis([-2 20]); title('Invivo FF HB'); colorbar;
    subplot(2,4,4);
    imagesc(transpose(ff_map_match_hb(:,:,slc_for_plot))*100); caxis([-2 20]); title('Exvivo FF HB'); colorbar;

    subplot(2,4,5);
    imagesc(r2star(:,:,slc_for_plot)); caxis([10 80]); title('Invivo R2star'); colorbar;
    subplot(2,4,7);
    imagesc(r2star_map_match(:,:,slc_for_plot)); caxis([10 80]); title('Exvivo R2star'); colorbar;
    subplot(2,4,6);
    imagesc(r2star_hb(:,:,slc_for_plot)); caxis([10 80]); title('Invivo R2star HB'); colorbar;
    subplot(2,4,8);
    imagesc(transpose(R2s_match_hb(:,:,slc_for_plot))); caxis([10 80]); title('Exvivo R2star HB'); colorbar;
    colormap(brewermap([],'*RdYlBu'));

    figure();
    subplot(2,4,1);
    imagesc(ff(:,:,slc_for_plot).*myo_ff_nan(:,:,slc_for_plot)); caxis([-2 20]); title('Invivo FF'); colorbar;
    subplot(2,4,3);
    imagesc(ff_map_match(:,:,slc_for_plot).*mask_myocardium_match(:,:,slc_for_plot)); caxis([-2 20]); title('Exvivo FF'); colorbar;
    subplot(2,4,2);
    imagesc(ff_hb(:,:,slc_for_plot)*100.*myo_ff_hb_nan(:,:,slc_for_plot)); caxis([-2 20]); title('Invivo FF HB'); colorbar;
    subplot(2,4,4);
    imagesc(transpose(ff_map_match_hb(:,:,slc_for_plot))*100.*mask_myocardium_match(:,:,slc_for_plot)); caxis([-2 20]); title('Exvivo FF HB'); colorbar;


    figure();
    subplot(2,2,1);
    imagesc(ff_mask_invivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([-2 10]); title('Invivo FF'); colorbar;
    subplot(2,2,2);
    imagesc(ff_mask_exvivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([-2 10]); title('Exvivo FF'); colorbar;
    subplot(2,2,3);
    imagesc(ff_mask_invivo_3d_hb(:,:,slc_for_plot).* myo_ff_hb_nan(:,:,slc_for_plot)*100); caxis([-2 10]); title('Invivo FF HB'); colorbar;
    subplot(2,2,4);
    imagesc(ff_mask_exvivo_3d_hb(:,:,slc_for_plot).* myo_ff_hb_nan(:,:,slc_for_plot)*100); caxis([-2 10]); title('Exvivo FF HB'); colorbar;

    figure();
    subplot(2,2,1);
    imagesc(r2star_mask_invivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([10 50]); title('Invivo R2star'); colorbar;
    subplot(2,2,2);
    imagesc(r2star_mask_exvivo_3d(:,:,slc_for_plot).* myo_ff_nan(:,:,slc_for_plot)); caxis([10 50]); title('Exvivo R2star'); colorbar;
    subplot(2,2,3);
    imagesc(r2star_mask_invivo_3d_hb(:,:,slc_for_plot).* myo_ff_hb_nan(:,:,slc_for_plot)); caxis([10 50]); title('Invivo R2star HB'); colorbar;
    subplot(2,2,4);
    imagesc(r2star_mask_exvivo_3d_hb(:,:,slc_for_plot).* myo_ff_hb_nan(:,:,slc_for_plot)); caxis([10 50]); title('Exvivo R2star HB'); colorbar;
    %chords_invivo = flip(chords_invivo,1);
    %chords_invivo_r2star = flip(chords_invivo_r2star, 1);
    
    % Save as
    reg_info_f = cat(2, base_dir, '/data/', name, '/Maps_and_Chords.mat');
    all_in_here = struct;
    all_in_here.ff_mask_invivo_3d = ff_mask_invivo_3d;
    all_in_here.ff_mask_exvivo_3d = ff_mask_exvivo_3d;
    all_in_here.r2star_mask_invivo_3d = r2star_mask_invivo_3d;
    all_in_here.r2star_mask_exvivo_3d = r2star_mask_exvivo_3d;
    all_in_here.t1_mask_invivo_3d = t1_mask_invivo_3d;
    all_in_here.t2_mask_invivo_3d = t2_mask_invivo_3d;
    
    all_in_here.chords_invivo_t2 = chords_invivo_t2;
    all_in_here.chords_invivo_t1 = chords_invivo_t1;
    all_in_here.chords_r2star = chords_r2star;
    all_in_here.chords = chords;
     
    all_in_here.mi_label_invivo = mi_label_invivo;
    all_in_here.idx_reordered_inex = idx_reordered_inex;
    save(reg_info_f, '-struct', 'all_in_here');
end

%% Save
%% Bland-Altman
addpath('../function/BlandAltman/');
% Heart muscle territories per patient
% territories = {'9MO'};
territories = {time_point};
nterritories = length(territories);

% Patient states during measurement
states = {'Remote', 'MI'};
nstates = length(states);

mi_label_invivo_reshape = mi_label_invivo(:);
chords_invivo_reshape = chords_invivo(:);
remote_label_invivo_reshape = double(~mi_label_invivo_reshape);

mi_label_invivo_reshape(mi_label_invivo_reshape == 0) = nan;
remote_label_invivo_reshape(remote_label_invivo_reshape == 0) = nan;

chords_invivo_reshape_mi = mi_label_invivo_reshape .* chords_invivo_reshape;
chords_invivo_reshape_remote = remote_label_invivo_reshape .* chords_invivo_reshape;

% Exvivo
chords_exvivo = chords(:,:,idx_reordered_inex);
chords_exvivo_reshape = chords_exvivo(:);
chords_exvivo_reshape_mi = mi_label_invivo_reshape .* chords_exvivo_reshape;
chords_exvivo_reshape_remote = remote_label_invivo_reshape .* chords_exvivo_reshape;

% Baseline data with noise
data1 = cat(3, chords_invivo_reshape_remote, chords_invivo_reshape_mi);
data2 = cat(3, chords_exvivo_reshape_remote, chords_exvivo_reshape_mi);

% BA plot paramters
tit = 'Invivo-Exvivo Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'Invivo Fat Fraction','Exvivo Fat Fraction','%'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 1 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [0 0 1;... % or RGB triplets
		      1 0 0];
end

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on');


data1_mi = cat(3, chords_invivo_reshape_mi);
data2_mi = cat(3, chords_exvivo_reshape_mi);


% BA plot paramters
states_mi = {'MI'};
tit = 'Invivo-Exvivo Comparison'; % figure title
gnames = {territories, states_mi}; % names of groups in data {dimension 1 and 2}
label = {'Invivo Fat Fraction','Exvivo Fat Fraction','%'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
colors = 'r';      % character codes

[cr_mi, fig_mi, statsStruct_mi] = BlandAltman(data1_mi, data2_mi,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on');

% % Generate figure with numbers of the data points (patients) and fixed
% % Bland-Altman difference data axes limits
% BlandAltman(data1, data2,label,[tit ' (numbers, forced 0 intercept, and fixed BA y-axis limits)'],gnames,'corrInfo',corrinfo,'axesLimits',limits,'colors',colors,'symbols','Num','baYLimMode','square','forceZeroIntercept','on')
% 
% 
% % Generate figure with differences presented as percentages and confidence
% % intervals on the correlation plot
% BAinfo = {'RPC'};
% BlandAltman(data1, data2,label,[tit ' (show fit confidence intervals and differences as percentages)'],gnames,'diffValueMode','percent', 'showFitCI','on','baInfo',BAinfo)
% 
% 
% % Display statistical results that were returned from analyses
% disp('Statistical results:');
% disp(statsStruct);