close all;
clear all;

addpath('function/');
base_dir = uigetdir; % Example/
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = {'11D05'};
time_points = {'7D', '8WK', '12WK', '14WK' '4MO', '6MO', '9MO', '1YR', '15YR'};

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

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

%% Pull LGE, T1 Map
label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis = struct;


for n = 1:length(Names)
% for n = 1:1
% for n = starting_point:starting_point
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    
    for tp = 1:length(time_points)
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            myo_glob = glob(cat(2, tp_dir, label_t1, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{6}, '/*'));
            
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
            [slc_array_t1, idx_reordered] = sort(slc_array);
            
            roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
            roi_t1 = roi_in_myo_t1 .* vol_img_3D;
            remote_t1 = remote_in_myo_t1 .* vol_img_3D;
            t1 = vol_img_3D;
            myo_t1 = mask_myocardium_3D;
            
            roi_in_myo_t1 = roi_in_myo_t1(:,:,idx_reordered);
            remote_in_myo_t1 = remote_in_myo_t1(:,:,idx_reordered);
            roi_t1 = roi_t1(:,:,idx_reordered);
            remote_t1 = remote_t1(:,:,idx_reordered);
            t1 = t1(:,:,idx_reordered);
            myo_t1 = myo_t1(:,:,idx_reordered);
            
            % LGE
            myo_glob = glob(cat(2, tp_dir, label_lge, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{6}, '/*'));
            
            load(cat(2, tp_dir, label_lge, '/', label_lge, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_lge, '/', label_lge, '_SliceLoc.mat'));
            slc_array_lge = slc_array;
            
            roi_in_myo_lge = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_lge = mask_myocardium_3D .* myoRefMask_3D;
            roi_lge = roi_in_myo_lge .* vol_img_3D;
            remote_lge = remote_in_myo_lge .* vol_img_3D;
            lge = vol_img_3D;
            myo_lge = mask_myocardium_3D;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_lge);
            lge = lge(:,:,idx_reordered);
            myo_lge = myo_lge(:,:,idx_reordered);
            remote_lge = remote_lge(:,:,idx_reordered);
            roi_lge = roi_lge(:,:,idx_reordered);
            remote_in_myo_lge = remote_in_myo_lge(:,:,idx_reordered);
            roi_in_myo_lge = roi_in_myo_lge(:,:,idx_reordered);
            
         
            for i = 1:size(roi_t1, 3)
                figure();
                subplot(1,2,1);
                imagesc(t1(:,:,i)); axis image; colorbar;
                title(['T1 Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                colormap(brewermap([],'*RdYlBu'));
                subplot(1,2,2);
                imagesc(lge(:,:,i)); axis image; colorbar;
                title(['LGE: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                
                overview_dir = cat(2, name_save_dir, '/overview/');
                if ~exist(overview_dir, 'dir')
                    mkdir(overview_dir);
                end
                saveas(gcf, cat(2, overview_dir, name, '_', time_point, '_Slice', num2str(i), '.png'));
            end

            for i = 1:size(roi_t1, 3)
                
                
                figure();
                subplot(1,2,1);
                roi_edg_t1  = edge(squeeze(roi_in_myo_t1(:,:,i)),'Canny');
                remote_edg_t1  = edge(squeeze(remote_in_myo_t1(:,:,i)),'Canny');
                myo_edg_t1 = edge(squeeze(myo_t1(:,:,i)), 'Canny');
                RGB_t1 = Func_Display_As_RGB(t1(:,:,i), roi_edg_t1, remote_edg_t1, myo_edg_t1);
                
                imagesc(RGB_t1); axis image;
                title(['T1 Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                %colormap(brewermap([],'*RdYlBu'));
                
                subplot(1,2,2);
                roi_edg_lge  = edge(squeeze(roi_in_myo_lge(:,:,i)),'Canny');
                remote_edg_lge  = edge(squeeze(remote_in_myo_lge(:,:,i)),'Canny');
                myo_edg_lge = edge(squeeze(myo_lge(:,:,i)), 'Canny');
                RGB_lge = Func_Display_As_RGB(lge(:,:,i), roi_edg_lge, remote_edg_lge, myo_edg_lge);
                
                imagesc(RGB_lge); axis image; 
                title(['LGE: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                
                overview_dir = cat(2, name_save_dir, '/overview/');
                if ~exist(overview_dir, 'dir')
                    mkdir(overview_dir);
                end
                saveas(gcf, cat(2, overview_dir, name, '_', time_point, '_Contour_Slice', num2str(i), '.png'));
            end
            % close all;
            
            name_tp = cat(2, name, '_', time_point);
            for_analysis(n).Name = name;
            for_analysis(n).metrics = struct;
            for_analysis(n).metrics.time_point = time_point;
            for_analysis(n).metrics.mean_roi_t1 = mean(nonzeros(roi_t1));
            for_analysis(n).metrics.sd_roi_t1 = std(nonzeros(roi_t1));
            for_analysis(n).metrics.mean_remote_t1 = mean(nonzeros(remote_t1));
            for_analysis(n).metrics.sd_remote_t1 = std(nonzeros(remote_t1));
            
        end
        
        
%         figure('Position', [1000 500 900 600]);
%         errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
%         hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
%         xticklabels(time_points);
%         set(gca, 'FontSize', 16);
%         grid on;
%         legend({'MI', 'Remote'}, 'Location', 'SouthEast');
        
    end
end

%% Example for AHA segmentation
% Simple showcase
Imgin = t1;
Maskin = myo_t1;
Segn = 50;
Groove = 0;

[Segmentpix, stats, Mask_Segn] =AHASegmentation(Imgin,Maskin,Segn,Groove);
figure(); 
for i = 1:size(Mask_Segn, 3)
   subplot(2,2,i);
   imagesc(Mask_Segn(:,:,i)); axis image;
end
%% AHA Segmentation and Plot BullsEye
% ========================== To simplify, name the reference point at (100, 100)
% ========================== To simplify, centroid is centroid of the first slice 
x = 100;
y = 100;
C = regionprops(myo_t1(:,:,1));
x_centroid = C.Centroid(2);
y_centroid = C.Centroid(1);
BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;

n = size(t1, 3);
mode = mod(n,3);
integ = fix(n/3);
if n >= 3
    switch mode
        case {0}
            aha_slice = cat(2, repmat([1], [1, integ]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {1}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {2}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ+1]), repmat([3], [1, integ]));
    end
else
    error("Available slice numbers are smaller than 3.");
end

% Basal
LocPixCount1 = zeros(6, 1);
SegTotalPixCount1 = zeros(6, 1);
Groove = BaseGroove + 60;
[~, idx_array] = sort(slc_array_t1); % assume minimal value is basal slice location
basal_idx = idx_array(aha_slice == 1);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(t1(:,:,basal_idx), myo_t1(:,:,basal_idx), 6, Groove);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount1(i) = LocPixCount1(i) + sum(Segmentpix{i,j});
        SegTotalPixCount1(i) = SegTotalPixCount1(i) + length(Segmentpix{i,j});
    end
end

% Mid-ventricular
LocPixCount2 = zeros(6, 1);
SegTotalPixCount2 = zeros(6, 1);
Groove = BaseGroove + 60;
mid_idx = idx_array(aha_slice == 2);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(t1(:,:,mid_idx), myo_t1(:,:,mid_idx), 6, Groove);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount2(i) = LocPixCount2(i) + sum(Segmentpix{i,j});
        SegTotalPixCount2(i) = SegTotalPixCount2(i) + length(Segmentpix{i,j});
    end
end

% Apical
LocPixCount3 = zeros(4, 1);
SegTotalPixCount3 = zeros(4, 1);
Groove = BaseGroove + 75;
apical_idx = idx_array(aha_slice == 2);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(t1(:,:,apical_idx), myo_t1(:,:,apical_idx), 4, Groove);

for i = 1:4
    for j = 1:size(Segmentpix, 2)
        LocPixCount3(i) = LocPixCount3(i) + sum(Segmentpix{i,j});
        SegTotalPixCount3(i) = SegTotalPixCount3(i) + length(Segmentpix{i,j});
    end
end

SegPixCount = [LocPixCount1; LocPixCount2; LocPixCount3];
SegTotalPixCount = [SegTotalPixCount1; SegTotalPixCount2; SegTotalPixCount3];
SegPixPerc = SegPixCount ./ SegTotalPixCount;


% =================== Create the AHA 17-segment bullseye =================
% ========================================================================
figure('Position', [400 400 800 800]);
PlotBullsEye(SegPixPerc);
