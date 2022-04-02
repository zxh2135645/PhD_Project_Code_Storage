clear all;
close all;

%% Pixel Image recon: (the main body)
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
%time_points = {'6MO', '9MO', '1YR', '15YR'};

dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    };

sequence_label = {'T1', 'T2star', 'LGE'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
label_t2 = 'T2';

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Before analysis, parse pre_QualControl
load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));

status_check = struct;
for n = 1:length(Names)
    name = Names{n};
    status_check(n).Name = name;
    status_check(n).status = [];
    status_check(n).status_final = [];
    tp_count = 1;
    for tp = 1:length(time_points)
        
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
                slc_loc = cat(2, 'Slice', num2str(i));
                if pre_QualControl(n).status(end-tp+1).(slc_loc) == 1
                    status_check(n).status(tp_count, i) = 1;
                elseif pre_QualControl(n).status(end-tp+1).(slc_loc) == 0
                    status_check(n).status(tp_count, i) = 0;
                end
            end
            
            tp_count = tp_count + 1;
        end
    end
    for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
        if i <= size(status_check(n).status,2)
            status_check(n).status_final(1,i) = all(status_check(n).status(:,i));
            % The timepoints of status_check goes from end to beginning
        end
    end
end

%%
fibrosis_group = {};
mean_ff_array_f = [];
mean_t1_array_f = [];
mean_r2star_array_f = [];

Segn = 50;
Groove = 0;     
se = strel('disk', 1);
%figure();
for n = 1:length(Names)
%for n = 1:1
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    name_data_save_dir = cat(2, data_save_dir, name);
    if ~exist(name_data_save_dir, 'dir')
        mkdir(name_data_save_dir);
    end
    tp_count = 0;
    for tp = 1:length(time_points)
    %for tp = 4:4
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            tp_count = tp_count+1;
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
            [slc_array_t1, idx_reordered] = sort(slc_array);
            
            roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
            roi_t1 = roi_in_myo_t1 .* vol_img_3D;
            remote_t1 = remote_in_myo_t1 .* vol_img_3D;
            t1 = vol_img_3D;
            myo_t1 = mask_myocardium_3D;
            hemo_core_t1 = noReflowMask_3D .* mask_myocardium_3D;
            
            roi_in_myo_t1 = roi_in_myo_t1(:,:,idx_reordered);
            remote_in_myo_t1 = remote_in_myo_t1(:,:,idx_reordered);
            roi_t1 = roi_t1(:,:,idx_reordered);
            remote_t1 = remote_t1(:,:,idx_reordered);
            t1 = t1(:,:,idx_reordered);
            myo_t1 = myo_t1(:,:,idx_reordered);
            hemo_core_t1 = hemo_core_t1(:,:,idx_reordered);
            
            % FF
            myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
            
            load(myo_glob{1});
            load(roi_glob{1});
            load(roi_glob{2});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat'));
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
            slc_array_ff = slc_array;
            
            ff_map = cell(1, length(glob_names));
            for f = 1:length(ff_map)
                ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
            end
            
            % convert ff_map to matrix
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            for f = 1:length(ff_map)
                ff(:,:,f) = ff_map{f}.fwmc_ff;
            end
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            hemo_core_t2star = mask_myocardium_3D .* freeROIMask_3D_te8;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            ff = ff(:,:,idx_reordered);
            myo_ff = myo_ff(:,:,idx_reordered);
            remote_ff = remote_ff(:,:,idx_reordered);
            roi_ff = roi_ff(:,:,idx_reordered);
            remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
            roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            hemo_core_t2star = hemo_core_t2star(:,:,idx_reordered);
            
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
            
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;
            
            r2star = r2star(:,:,idx_reordered);
            myo_r2star = myo_r2star(:,:,idx_reordered);
            remote_r2star = remote_r2star(:,:,idx_reordered);
            roi_r2star = roi_r2star(:,:,idx_reordered);
            remote_in_myo_r2star = remote_in_myo_r2star(:,:,idx_reordered);
            roi_in_myo_r2star = roi_in_myo_r2star(:,:,idx_reordered);
            
            % T2
            load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
            myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));
            noreflow_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{7}, '/*'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(noreflow_glob{1});
            load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat'));
            [slc_array_t2, idx_reordered] = sort(slc_array);
            
            roi_in_myo_t2 = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_t2 = mask_myocardium_3D .* myoRefMask_3D;
            roi_t2 = roi_in_myo_t2 .* vol_img_3D;
            remote_t2 = remote_in_myo_t2 .* vol_img_3D;
            t2 = vol_img_3D;
            myo_t2 = mask_myocardium_3D;
            hemo_core_t2 = noReflowMask_3D .* mask_myocardium_3D;
            
            roi_in_myo_t2 = roi_in_myo_t2(:,:,idx_reordered);
            remote_in_myo_t2 = remote_in_myo_t2(:,:,idx_reordered);
            roi_t2 = roi_t2(:,:,idx_reordered);
            remote_t2 = remote_t2(:,:,idx_reordered);
            t2 = t2(:,:,idx_reordered);
            myo_t2 = myo_t2(:,:,idx_reordered);
            hemo_core_t2 = hemo_core_t2(:,:,idx_reordered);
            
            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            
            LR_mdl_fname = cat(2, name_data_save_dir, '/LinearRegression_', name, '_', time_point, '.mat');

            status = status_check(n).status(tp_count,:);
            
            for slc = 1:size(hemo_core_t2star,3)
                if status(slc) == 1
                    temp = hemo_core_t2star(:,:,slc);
                    if ~any(temp(:))
                        fibrosis_group{end+1} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
%                         figure();
%                         subplot(2,2,1);
%                         imagesc(t1(:,:,slc)); axis image; caxis([0, 1500]);
%                         title(cat(2, name, '_', time_point, '_Slice', num2str(slc)));
%                         subplot(2,2,2);
%                         imagesc(r2star(:,:,slc)); axis image; caxis([0 100]);
%                         subplot(2,2,3);
%                         imagesc(ff(:,:,slc)); axis image; caxis([0 100]);
                        
                        %Func_T1FP_Chord_ReAnalysis_EndoEpi(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, tp_dir2, name, time_point, LR_mdl_fname, chord_values_fname, chord_values_fname2,status);
                        chord_values_fname = cat(2, name_data_save_dir, '/Chord_values_pixelwise_', name, '_', time_point, '.mat');
                        chord_values_fname2 = cat(2, name_data_save_dir, '/Chord_values2_', name, '_', time_point, '.mat');
                        
                        load(chord_values_fname2);
                        mean_t1_array_endo_slc = mean_t1_array_endo(:,slc);
                        mean_t1_array_epi_slc = mean_t1_array_epi(:,slc);
                        mean_ff_array_endo_slc = mean_ff_array_endo(:,slc);
                        mean_ff_array_epi_slc = mean_ff_array_epi(:,slc);
                        mean_r2star_array_endo_slc = mean_r2star_array_endo(:,slc);
                        mean_r2star_array_epi_slc = mean_r2star_array_epi(:,slc);
                        
                        label_endo = (mean_t1_array_endo_slc ~= 0 & ~isnan(mean_t1_array_endo_slc));
                        label_epi = (mean_t1_array_epi_slc ~= 0 & ~isnan(mean_t1_array_endo_slc));
                        
                        mean_ff_array_slc = [mean_ff_array_endo_slc(label_endo);mean_ff_array_epi_slc(label_epi)];
                        mean_t1_array_slc = [mean_t1_array_endo_slc(label_endo);mean_t1_array_epi_slc(label_epi)];
                        mean_r2star_array_slc = [mean_r2star_array_endo_slc(label_endo);mean_r2star_array_epi_slc(label_epi)];
                        
                        mean_ff_array_f = [mean_ff_array_f; mean_ff_array_slc];
                        mean_t1_array_f = [mean_t1_array_f; mean_t1_array_slc];
                        mean_r2star_array_f = [mean_r2star_array_f; mean_r2star_array_slc];
                        
                        figure();
                        scatter(mean_ff_array_slc, mean_t1_array_slc);
                        hold on;
                        
                        mean_remote = mean(nonzeros(remote_t1(:,:,slc)));
                        % myo_t1_eroded = Mask_Segn_epi+Mask_Segn_endo>0;
                        figure(); 
                        subplot(2,2,1); imagesc(t1(:,:,slc).*myo_t1(:,:,slc)); axis image; caxis([1000 1800]);
                        subplot(2,2,2); imagesc(t2(:,:,slc).*myo_t2(:,:,slc)); axis image; caxis([-200 1000]);
                        subplot(2,2,3); imagesc(ff(:,:,slc).*myo_ff(:,:,slc)); axis image; caxis([-10 50]);
                        subplot(2,2,4); imagesc(r2star(:,:,slc).*myo_r2star(:,:,slc)); axis image; caxis([-20 100]);
                        
                        univ_myo = myo_t1;
                        fixed_eroded = imerode(myo_t1(:,:,slc), se);
                        BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
                        center_fixed = imfill(BW_skel, 'hole');
                        center_fixed = imopen(center_fixed, se); % Removing spikes
                        fixed_epi = fixed_eroded - center_fixed > 0;
                        fixed_endo = center_fixed + fixed_eroded > 1;
                        
                        [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(t1(:,:,slc), fixed_endo, Segn, Groove);
                        [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(t1(:,:,slc), fixed_epi, Segn, Groove);
                        
                        idx_endo = find(label_endo);
                        idx_epi = find(label_epi);
                    end
                end
            end
            
            % AHA Segment
            %Segn = 50;
            %Groove = 0;


%             if (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
%                 disp(cat(2, 'Skipped: ', name, ' ', time_point))
%             else
% 
%             end
        end
    end
    % close all;
end

%%
figure(); scatter(mean_ff_array_f, mean_t1_array_f);
mdl = fitlm(mean_ff_array_f, mean_t1_array_f);

figure(); scatter(mean_r2star_array_f, mean_t1_array_f);
mdl2 = fitlm(mean_r2star_array_f, mean_t1_array_f);

figure(); scatter(mean_ff_array_f, mean_r2star_array_f);
mdl3 = fitlm(mean_ff_array_f, mean_r2star_array_f);

%% 
t1myo = mean_remote; % ms
t1fib = 1150:10:1500; % ms
%mean_ff_array_f
%mean_t1_array_f
w1 = zeros(length(mean_ff_array_f), length(t1fib));
w2 = zeros(length(mean_ff_array_f), length(t1fib));

for i = 1:length(t1fib)
    %w1(:,i) = ((1 - mean_ff_array_f/100) .* t1fib(i) - mean_t1_array_f ./ (1 + mean_ff_array_f/100 /100)) ./ (t1fib(i) - t1myo);
    w2(:,i) = 1 - w1(:,i) - mean_ff_array_f/100;
end
%%
t1myo = mean_remote; % ms
%t1fib = 1150:10:1500; % ms
%t1_fib
%mean_ff_array_f
%mean_t1_array_f
w1_2 = 0:0.01:1;
%w1 = zeros(length(mean_ff_array_f), length(t1fib));
%w2 = zeros(length(mean_ff_array_f), length(t1fib));
t1_fib = zeros(length(mean_ff_array_f), length(w1_2));
for i = 1:length(w1_2)
    %w1(:,i) = ((1 - mean_ff_array_f/100) .* t1fib(i) - mean_t1_array_f ./ (1 + mean_ff_array_f/100 /100)) ./ (t1fib(i) - t1myo);
    %w2(:,i) = 1 - w1(:,i) - mean_ff_array_f/100;
    t1_fib(:,i) = (w1_2(i) * t1myo - mean_t1_array_f ./ (1 + mean_ff_array_f/100/100)) ./ (w1_2(i) - (1 - mean_ff_array_f/100));
end

t1_fib(t1_fib>2000) = nan;
t1_fib(t1_fib<500) = nan;
t1_fib_max = max(t1_fib, [], 2);
t1_fib_min = min(t1_fib, [], 2);

y = (1:1:18)';
figure(); 
line([t1_fib_min,t1_fib_max]', [y,y]'); xlim([400 2000]);
hold on;
%scatter([t1_fib_min,t1_fib_max]', [y,y]');
%%
mean_ff_array_endo_slc2 = mean_ff_array_endo_slc;
mean_ff_array_endo_slc2(~label_endo) = nan;
mean_ff_array_epi_slc2 = mean_ff_array_epi_slc;
mean_ff_array_epi_slc2(~label_epi) = nan;
mean_t1_array_endo_slc2 = mean_t1_array_endo_slc;
mean_t1_array_endo_slc2(~label_endo) = nan;
mean_t1_array_epi_slc2 = mean_t1_array_epi_slc;
mean_t1_array_epi_slc2(~label_epi) = nan;

MaskSegn_w1 = double((Mask_Segn_endo + Mask_Segn_epi) > 0);
MaskSegn_w1(MaskSegn_w1==0) = nan;
MaskSegn_w1(MaskSegn_w1==1) = -0.1;
MaskSegn_w2 = MaskSegn_w1;
MaskSegn_t1 = MaskSegn_w1;
MaskSegn_t1(MaskSegn_t1==-0.1) = 1000;
MaskSegn_ff = MaskSegn_w1;
MaskSegn_ff(MaskSegn_ff==-0.1) = -3

for i = 1:length(idx_endo)
    MaskSegn_w1(Mask_Segn_endo == idx_endo(i)) = w1_array(idx_endo(i),1);
    MaskSegn_w2(Mask_Segn_endo == idx_endo(i)) = w2_array(idx_endo(i),1);
    MaskSegn_ff(Mask_Segn_endo == idx_endo(i)) = mean_ff_array_endo_slc(idx_endo(i),1);
    MaskSegn_t1(Mask_Segn_endo == idx_endo(i)) = mean_t1_array_endo_slc(idx_endo(i),1);
end
for i = 1:length(idx_epi)
    MaskSegn_w1(Mask_Segn_epi == idx_epi(i)) = w1_array(idx_epi(i),2);
    MaskSegn_w2(Mask_Segn_epi == idx_epi(i)) = w2_array(idx_epi(i),2);
    MaskSegn_ff(Mask_Segn_epi == idx_epi(i)) = mean_ff_array_epi_slc(idx_epi(i),1);
    MaskSegn_t1(Mask_Segn_epi == idx_epi(i)) = mean_t1_array_epi_slc(idx_epi(i),1);
end
figure(); 
subplot(2,2,1); imagesc(MaskSegn_w1); caxis([-0.2 1]);colormap(brewermap([],'*RdYlBu'));axis off; colorbar; title('Viable Myo');
subplot(2,2,2); imagesc(MaskSegn_w2); caxis([-0.2 1]);axis off;colorbar;title('Fibrosis');
subplot(2,2,3); imagesc(MaskSegn_ff); caxis([-5 25]);axis off;colorbar;title('Fat Fraction (%)');
subplot(2,2,4); imagesc(MaskSegn_t1); caxis([800 1500]);axis off;colorbar;title('T1 (ms)');

figure(); 
subplot(2,2,1); imagesc(w1_array); caxis([-0.2 1]);colormap(brewermap([],'*RdYlBu'));axis off; colorbar; title('Viable Myo');
subplot(2,2,2); imagesc(w2_array); caxis([-0.2 1]);axis off;colorbar;title('Fibrosis');
subplot(2,2,3); imagesc([mean_ff_array_endo_slc2, mean_ff_array_epi_slc2]); caxis([-5 25]);axis off;colorbar;title('Fat Fraction (%)');
subplot(2,2,4); imagesc([mean_t1_array_endo_slc2, mean_t1_array_epi_slc2]); caxis([800 1500]);axis off;colorbar;title('T1 (ms)');

figure();
subplot(2,2,1); imagesc(t1(:,:,slc).*myo_t1(:,:,slc)); axis image; caxis([1000 1800]);
subplot(2,2,2); imagesc(t2(:,:,slc).*myo_t2(:,:,slc)/10); axis image; caxis([-10 50]);
subplot(2,2,3); imagesc(ff(:,:,slc).*myo_ff(:,:,slc)); axis image; caxis([-10 50]);
subplot(2,2,4); imagesc(r2star(:,:,slc).*myo_r2star(:,:,slc)); axis image; caxis([-20 100]);
%% 
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = 640;
imageSizeY = 640;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 320;
centerY = 320;
radius = 100;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
radius = 50;
circlePixels_2 = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
figure();
image(circlePixels&~circlePixels_2) ;
colormap([0 0 0; 1 1 1]);
title('Binary image of a circle');