clear all;
close all;
%% 
% Merry 1 YR
%roi_ff_array = data_storage_rim(1).data(2).slices.roi_ff_array;
%roi_t1_array = data_storage_rim(1).data(2).slices.roi_t1_array;
% roi_ff_array = [MI_Chord_Analysis2(4).Mipix_mean2_hemo_n; MI_Chord_Analysis2(4).Mipix_mean2_hemo_p];
% roi_t1_array = [MI_Chord_Analysis2(4).Mipix_mean_hemo_n;MI_Chord_Analysis2(4).Mipix_mean_hemo_p];
% %%
% figure();
% % scatter(roi_ff_array(:), roi_t1_array(:));
% for i = 1:4
% subplot(2,2,i);
% scatter(roi_ff_array(:,i), roi_t1_array(:,i))
% end
%% Re-do AHA analysis: (the main body)
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};

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


%for n = 1:length(Names)
for n = 6:6
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            for f = 1:length(ff_map)
                ff(:,:,f) = ff_map{f}.fwmc_ff;
            end
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            ff = ff(:,:,idx_reordered);
            myo_ff = myo_ff(:,:,idx_reordered);
            remote_ff = remote_ff(:,:,idx_reordered);
            roi_ff = roi_ff(:,:,idx_reordered);
            remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
            roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            
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
            
            
            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            LR_mdl_fname = cat(2, name_data_save_dir, '/LinearRegression_', name, '_', time_point, '.mat');
            %chord_values_fname = cat(2, name_data_save_dir, '/Chord_values_', name, '_', time_point, '.mat');
            chord_values_fname = cat(2, name_data_save_dir, '/Chord_values_pixelwise_', name, '_', time_point, '.mat');
            chord_values_fname2 = cat(2, name_data_save_dir, '/Chord_values2_', name, '_', time_point, '.mat');
            % if condition
            center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');
            
            strat_fname = cat(2, name_data_save_dir, '/FF_Stratify_', name, '_', time_point, '.mat');
            
            %if ~exist(center_mask_fname, 'file')
                center_mask_ff = zeros(size(roi_in_myo_ff));
                BW_skel = zeros(size(roi_in_myo_ff));
                %caxis_rg = [0 1];
                %for i = 1:size(roi_in_myo_ff, 3)
                    % Draw center line on myocardium ff
                    %disp(cat(2, 'Center line roi: ', name, '  ', time_point, 'Slice', num2str(i)));

                    %moving = myo_ff(:,:,i);
                    %center_mask_ff(:,:,i) = Func_DrawCenterLine(moving, caxis_rg);
                    %BW_skel(:,:,i) = bwmorph(moving, 'skel', Inf);
                    %center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
                %end
                
                figure();
                sz = ceil(sqrt(size(roi_in_myo_ff, 3)));
                se = strel('disk', 1);
                for i = 1:size(roi_in_myo_ff, 3)
                    subplot(sz,sz,i);
                    myo_ff_eroded = imerode(myo_ff(:,:,i), se);
                    BW_skel(:,:,i) = bwmorph(myo_ff_eroded, 'skel', Inf);
                    center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
                    
                    % Found 18D16 8WK is not a closed shape (Hard-Coded)
                    % Felicity 6MO
                    if (n == 14 && tp == 8 && i == 3) || (n == 11 && tp == 4 && i == 3)
                        center_mask_ff(:,:,3) = bwconvhull(center_mask_ff(:,:,3));
                    end
                    epi = myo_ff_eroded - center_mask_ff(:,:,i) > 0;
                    endo = center_mask_ff(:,:,i) + myo_ff_eroded > 1;
                    imagesc(endo*2 + epi);
                    colormap(brewermap([],'*RdYlBu'));
                end
                save(center_mask_fname, 'center_mask_ff');
                saveas(gcf, cat(2, tp_dir2, 'CenterLineMask.png'))
            %else
            %    load(center_mask_fname);
            %end
            status = status_check(n).status(tp_count,:);
            % AHA Segment
            Segn = 50;
            Groove = 0;
            % MI_Chord_Analysis_fname = cat(2, name_data_save_dir, '/MIChordAnalysis_', name, '_', time_point, '.mat');
            % Remote_Chord_Analysis_fname = cat(2, name_data_save_dir, '/RemoteChordAnalysis_', name, '_', time_point, '.mat');
            if (strcmp(name, 'Queenie') && strcmp(time_point, '7D')) % || (strcmp(name, 'Evelyn') && strcmp(time_point, '6MO'))
                %Func_T1FP_Chord_ReAnalysis2_EndoEpi(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, tp_dir2, name, time_point, LR_mdl_fname, chord_values_fname, chord_values_fname2,status);
                %Func_T1FP_Chord_ReAnalysis2_Pixelwise(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,LR_mdl_fname,chord_values_fname,status);
                Func_T1FP_Stratification_Analysis2(t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,strat_fname,status);
            elseif (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
                disp(cat(2, 'Skipped: ', name, ' ', time_point))
            else
                Func_T1FP_Chord_ReAnalysis(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, tp_dir2, name, time_point, LR_mdl_fname, chord_values_fname);
                Func_T1FP_Chord_ReAnalysis_EndoEpi(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, tp_dir2, name, time_point, LR_mdl_fname, chord_values_fname, chord_values_fname2,status);
                Func_T1FP_Chord_ReAnalysis_Pixelwise(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,LR_mdl_fname,chord_values_fname,status);
                Func_T1FP_Stratification_Analysis(t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,strat_fname,status);
            end
            %[MI_Chord_Analysis, center_mask_ff] = Func_MI_Chord_Analysis(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff,...
            %    roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, MI_Chord_Analysis_fname);
            %[MI_Chord_Analysis2, ~] = Func_MI_Chord_Analysis2(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff,...
            %    roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, MI_Chord_Analysis_fname);
            %[Remote_Chord_Analysis2, ~] = Func_Remote_Chord_Analysis2(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff,...
            %    remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, Remote_Chord_Analysis_fname);
            % Plot figures and save
            %Func_plot_chord_analysis_general(MI_Chord_Analysis, tp_dir2);
            %Func_plot_chord_analysis_EpiEndo(MI_Chord_Analysis, tp_dir2);
            %Func_plot_chord_analysis_general2(MI_Chord_Analysis2, tp_dir2);
            %Func_plot_chord_analysis_EpiEndo2(MI_Chord_Analysis2, tp_dir2);
            
            
        end
    end
    close all;
end

%% From Debug_MI_Chord_Analysis.m

t1_cell = cell(50,size(roi_in_myo_t1, 3));
ff_cell = cell(50,size(roi_in_myo_t1, 3));
mean_t1_array = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array = zeros(50,size(roi_in_myo_t1, 3));
sd_ff_array = zeros(50,size(roi_in_myo_t1, 3));
sd_r2star_array = zeros(50,size(roi_in_myo_t1, 3));

mean_t1_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_endo = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array_endo = zeros(50,size(roi_in_myo_t1, 3));
sd_ff_array_endo = zeros(50,size(roi_in_myo_t1, 3));
sd_r2star_array_endo = zeros(50,size(roi_in_myo_t1, 3));  

mean_t1_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_epi = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array_epi = zeros(50,size(roi_in_myo_t1, 3));
sd_ff_array_epi = zeros(50,size(roi_in_myo_t1, 3));
sd_r2star_array_epi = zeros(50,size(roi_in_myo_t1, 3)); 

mean_t1_array_remote = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_remote = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_remote = zeros(50,size(roi_in_myo_t1, 3));
sd_t1_array_remote = zeros(50,size(roi_in_myo_t1, 3));
sd_ff_array_remote = zeros(50,size(roi_in_myo_t1, 3));
sd_r2star_array_remote = zeros(50,size(roi_in_myo_t1, 3));

for i = 1:size(roi_in_myo_t1, 3)
%for i = 1:1
img = t1(:,:,i);
img2 = ff(:,:,i);
img3 = r2star(:,:,i);
fixed = myo_t1(:,:,i);
moving = myo_ff(:,:,i);
%fixed = roi_in_myo_t1(:,:,i);
%moving = roi_in_myo_ff(:,:,i);

img2(img2 > 100) = 100;
img2(img2 < 0) = 0;

figure();
subplot(3,2,1);
imagesc(fixed); axis image; title('T1 (fixed)'); colormap(brewermap([],'*RdYlBu'));axis off;
subplot(3,2,2);
imagesc(moving); axis image; title('FF (moving)');axis off;

se = strel('disk', 1);
%[optimizer, metric] = imregconfig('multimodal');
%tform = imregtform(moving, fixed, 'rigid', optimizer, metric);

I1 = moving; I2 = fixed;
% Set static and moving image
S=I2; M=I1;

[movingRegistered,Bx,By,Fx,Fy] = register_images(M,S);
movingRegistered = movingRegistered > 0.5;
img2 = movepixels(img2,Bx,By);
img3 = movepixels(img3,Bx,By);

subplot(3,2,3); imagesc(img); axis image; title('T1 map'); axis off;
subplot(3,2,4); imagesc(img2); axis image; caxis([0 50]); title('FF map'); axis off;

movingRegistered_myo_ff = movepixels(myo_ff(:,:,i),Bx,By)>0.5;
movingRegistered_roi_ff = movepixels(roi_in_myo_ff(:,:,i),Bx,By)>0.5;
movingRegistered_roi_r2star = movepixels(roi_in_myo_ff(:,:,i),Bx,By)>0.5;
movingRegistered_remote_ff = movepixels(remote_in_myo_ff(:,:,i),Bx,By)>0.5;
%img2 = imwarp(img2,tform,'OutputView',imref2d(size(fixed)));
%img3 = imwarp(img3,tform,'OutputView',imref2d(size(fixed)));
%movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
% movingRegistered_myo_ff = imwarp(myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
% movingRegistered_roi_ff = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
% movingRegistered_roi_r2star = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
univ_roi = roi_in_myo_t1(:,:,i) & movingRegistered_roi_ff;

subplot(3,2,5);
%imagesc(movingRegistered+fixed);
imshowpair(fixed,movingRegistered,'Scaling','joint'); title('Registered');
subplot(3,2,6); imagesc(double(movingRegistered_myo_ff&fixed) + double(univ_roi) + 2*double(movingRegistered_remote_ff)); axis image;
title(cat(2, 'Slice = ', num2str(i)));
saveas(gcf, cat(2, tp_dir2, 'MyocardiumRegistration_demon_Slice', num2str(i), '.png'));

univ_myo = movingRegistered_myo_ff&fixed;
%movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;

fixed_eroded = imerode(myo_t1(:,:,i), se);
BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
center_fixed = imfill(BW_skel, 'hole');
center_fixed = imopen(center_fixed, se); % Removing spikes
fixedRegistered_epi = fixed_eroded - center_fixed > 0;
fixedRegistered_endo = center_fixed + fixed_eroded > 1;


movingRegistered_eroded = imerode(movingRegistered_myo_ff, se);
BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
center_moving = imfill(BW_skel, 'hole');
center_moving = imopen(center_moving, se); % Removing spikes
movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
movingRegistered_endo = center_moving + movingRegistered_eroded > 1;

%movingRegistered_roi_ff = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
%movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
%movingRegistered_epi = imwarp(myo_epi_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
%movingRegistered_endo = imwarp(myo_endo_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;

univ_myo_eroded = imerode(univ_myo, se);
% univ_myo_eroded = univ_myo;

% [Segmentpix, stats, Mask_Segn] = AHASegmentation(img, fixed_eroded, Segn, Groove);
% [Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, movingRegistered_eroded, Segn, Groove);
% [Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, movingRegistered_eroded, Segn, Groove);
[Segmentpix, stats, Mask_Segn] = AHASegmentation(img, univ_myo_eroded, Segn, Groove);
[Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, univ_myo_eroded, Segn, Groove);
[Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, univ_myo_eroded, Segn, Groove);


[Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(img, fixedRegistered_epi, Segn, Groove);
[Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(img, fixedRegistered_endo, Segn, Groove);
[Segmentpix, stats, Mask_Segn2_epi] = AHASegmentation(img2, movingRegistered_epi, Segn, Groove);
[Segmentpix, stats, Mask_Segn2_endo] = AHASegmentation(img2, movingRegistered_endo, Segn, Groove);
[Segmentpix, stats, Mask_Segn3_epi] = AHASegmentation(img3, movingRegistered_epi, Segn, Groove);
[Segmentpix, stats, Mask_Segn3_endo] = AHASegmentation(img3, movingRegistered_endo, Segn, Groove);


%figure('Position', [100 0 1600 1600]);
seg_mask_overlap = zeros(size(img));
%vid_save = cat(2, tp_dir2, 'myVideoFile_Slice', num2str(i),'.mov');
%myVideo = VideoWriter(vid_save); %open video file
%myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
%open(myVideo)

for j = 1:Segn
    
%     Mipix{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
%     Mipix_epi{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
%     Mipix_endo{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
%     
%     
%     Mipix2{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
%     Mipix2_epi{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
%     Mipix2_endo{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
%     
%     Mipix3{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
%     Mipix3_epi{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
%     Mipix3_endo{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
    
    Mipix{j,i} = img(Mask_Segn .* univ_roi .* fixed_eroded == j);
    Mipix_epi{j,i} = img(Mask_Segn_epi .* univ_roi .* fixed_eroded == j);
    Mipix_endo{j,i} = img(Mask_Segn_endo .* univ_roi .* fixed_eroded == j);
    
    
    Mipix2{j,i} = img2(Mask_Segn2 .* univ_roi .* movingRegistered_eroded == j);
    Mipix2_epi{j,i} = img2(Mask_Segn2_epi .* univ_roi .* movingRegistered_eroded == j);
    Mipix2_endo{j,i} = img2(Mask_Segn2_endo .* univ_roi .* movingRegistered_eroded == j);
    
    Mipix3{j,i} = img3(Mask_Segn3 .* univ_roi .* movingRegistered_eroded == j);
    Mipix3_epi{j,i} = img3(Mask_Segn3_epi .* univ_roi .* movingRegistered_eroded == j);
    Mipix3_endo{j,i} = img3(Mask_Segn3_endo .* univ_roi .* movingRegistered_eroded == j);

    remote_mean_t1 = mean(nonzeros(remote_in_myo_t1(:,:,i) .* img));
    remote_sd_t1 = std(nonzeros(remote_in_myo_t1(:,:,i) .* img));
    thresh = remote_mean_t1 - 2*remote_sd_t1;
    hemo_mask = (img<thresh).*univ_roi;
    
    % For debugging    
    seg_mask_t1 = Mask_Segn .* (univ_roi-hemo_mask) .* fixed_eroded == j;
    seg_mask_ff = Mask_Segn2 .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    seg_mask_r2star = Mask_Segn3 .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    
    seg_mask_t1_endo = Mask_Segn_endo .* (univ_roi-hemo_mask) .* fixed_eroded == j;
    seg_mask_t1_epi = Mask_Segn_epi .* (univ_roi-hemo_mask) .* fixed_eroded == j;
    
    seg_mask_ff_endo = Mask_Segn2_endo .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    seg_mask_ff_epi = Mask_Segn2_epi .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    
    seg_mask_r2star_endo = Mask_Segn3_endo .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    seg_mask_r2star_epi = Mask_Segn3_epi .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
    
    seg_mask_t1_remote = Mask_Segn .* movingRegistered_remote_ff .* fixed_eroded == j;
    seg_mask_ff_remote = Mask_Segn2 .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
    seg_mask_r2star_remote = Mask_Segn3 .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
    
    %{
    if ~isempty(nonzeros(seg_mask_t1)) && ~isempty(nonzeros(seg_mask_ff))
        seg_mask_overlap = seg_mask_overlap + seg_mask_t1 + 2*seg_mask_ff;
        ax1 = subplot(1,2,1);
        imagesc(ax1, seg_mask_overlap); axis image; caxis([0 3]);
        title(['Segment = ', num2str(j)]);
        colormap(brewermap([],'*RdYlBu'));
        ax2 = subplot(1,2,2);
        imagesc(ax2, seg_mask_t1 + 2*seg_mask_ff); axis image; caxis([0 3]);
        title(['Segment = ', num2str(j)]); %colorbar;
        colormap(brewermap([],'*RdYlBu'));
        pause(.5);
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
    end
    %}
    
    % 
    if any(seg_mask_t1(:)) 
        t1_cell{j,i} = nonzeros(img.*seg_mask_t1);
        ff_cell{j,i} = nonzeros(img2.*seg_mask_ff);
        mean_t1_array(j,i) = mean(nonzeros(img.*seg_mask_t1));
        mean_ff_array(j,i) = mean(nonzeros(img2.*seg_mask_ff));
        mean_r2star_array(j,i) = mean(nonzeros(img3.*seg_mask_r2star));
        sd_t1_array(j,i) = std(nonzeros(img.*seg_mask_t1));
        sd_ff_array(j,i) = std(nonzeros(img2.*seg_mask_ff));
        sd_r2star_array(j,i) = std(nonzeros(img3.*seg_mask_r2star));
        
        mean_t1_array_endo(j,i) = mean(nonzeros(img.*seg_mask_t1_endo));
        mean_t1_array_epi(j,i) = mean(nonzeros(img.*seg_mask_t1_epi));
        mean_ff_array_endo(j,i) = mean(nonzeros(img2.*seg_mask_ff_endo));
        mean_ff_array_epi(j,i) = mean(nonzeros(img2.*seg_mask_ff_epi));
        mean_r2star_array_endo(j,i) = mean(nonzeros(img3.*seg_mask_r2star_endo));
        mean_r2star_array_epi(j,i) = mean(nonzeros(img3.*seg_mask_r2star_epi));
        
        sd_t1_array_endo(j,i) = std(nonzeros(img.*seg_mask_t1_endo));
        sd_t1_array_epi(j,i) = std(nonzeros(img.*seg_mask_t1_epi));
        sd_ff_array_endo(j,i) = std(nonzeros(img2.*seg_mask_ff_endo));
        sd_ff_array_epi(j,i) = std(nonzeros(img2.*seg_mask_ff_epi));
        sd_r2star_array_endo(j,i) = std(nonzeros(img3.*seg_mask_r2star_endo));
        sd_r2star_array_epi(j,i) = std(nonzeros(img3.*seg_mask_r2star_epi));
        
    end
    if any(seg_mask_t1_remote(:))
        % remote
        mean_t1_array_remote(j,i) = mean(nonzeros(img.*seg_mask_t1_remote));
        mean_ff_array_remote(j,i) = mean(nonzeros(img2.*seg_mask_ff_remote));
        mean_r2star_array_remote(j,i) = mean(nonzeros(img3.*seg_mask_r2star_remote));
        sd_t1_array_remote(j,i) = std(nonzeros(img.*seg_mask_t1_remote));
        sd_ff_array_remote(j,i) = std(nonzeros(img2.*seg_mask_ff_remote));
        sd_r2star_array_remote(j,i) = std(nonzeros(img3.*seg_mask_r2star_remote));
    end
    
end
end

%% Register image

%[Ireg,Bx,By,Fx,Fy] = register_images(moving, fixed);
%% Plot one linear regression
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_ff_array_nz = nonzeros(mean_ff_array);
mean_t1_array_nz = nonzeros(mean_t1_array);
mean_ff_array_nz_new = mean_ff_array_nz;
mean_t1_array_nz_new = mean_t1_array_nz;
mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];
mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote);
mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote);
mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];

mdl = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);

figure('Position', [100 100 600 800]);
title(cat(2, name, ' ', time_point))
rows = ceil(size(roi_in_myo_t1, 3) / 2) + 1;
subplot(rows,1,rows);
scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_nz_new .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(2)-20, yl(1)+200, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

% figure();
% for each slice
mdl_slc = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz = nonzeros(mean_ff_array(:,i));
    mean_t1_array_nz = nonzeros(mean_t1_array(:,i));
    mean_ff_array_nz_new = mean_ff_array_nz;
    mean_t1_array_nz_new = mean_t1_array_nz;
    mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote(:,i));
    mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote(:,i));
    mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
    mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
    mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
    mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];

    mdl_slc = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);
    
    subplot(rows,2,i);
    scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_ff_array_nz_new .* mdl_slc.Coefficients.Estimate(2) + mdl_slc.Coefficients.Estimate(1);
    hold on;
    plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    text(xl(2)-30, yl(1)+100, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
end

%% Plot two linear regression
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_ff_array_nz = nonzeros(mean_ff_array);
mean_t1_array_nz = nonzeros(mean_t1_array);
mean_ff_array_nz_new = mean_ff_array_nz;
mean_t1_array_nz_new = mean_t1_array_nz;
mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];
mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote);
mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote);
mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];

[t1_max, idx_max] = max(mean_t1_array_nz_new);
cutoff = mean_ff_array_nz_new(idx_max);
mean_t1_array_nz_new1 = mean_t1_array_nz_new(mean_ff_array_nz_new <= cutoff);
mean_t1_array_nz_new2 = mean_t1_array_nz_new(mean_ff_array_nz_new >= cutoff);
mean_ff_array_nz_new1 = mean_ff_array_nz_new(mean_ff_array_nz_new <= cutoff);
mean_ff_array_nz_new2 = mean_ff_array_nz_new(mean_ff_array_nz_new >= cutoff);

mdl1 = fitlm(mean_ff_array_nz_new1, mean_t1_array_nz_new1);
mean_ff_array_nz_new2_1 = mean_ff_array_nz_new2(mean_ff_array_nz_new2<=35);
mean_t1_array_nz_new2_1 = mean_t1_array_nz_new2(mean_ff_array_nz_new2<=35);
mdl2 = fitlm(mean_ff_array_nz_new2_1, mean_t1_array_nz_new2_1);

figure('Position', [100 100 600 600]);
title(cat(2, name, ' ', time_point))
subplot(3,1,3);
scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y1 = mean_ff_array_nz_new1 .* mdl1.Coefficients.Estimate(2) + mdl1.Coefficients.Estimate(1);
Y2 = mean_ff_array_nz_new2 .* mdl2.Coefficients.Estimate(2) + mdl2.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_nz_new1, Y1, 'k', 'LineWidth', 1);
plot(mean_ff_array_nz_new2, Y2, 'k', 'LineWidth', 1);
xline(cutoff);
scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(1)+1, yl(2)-100, cat(2,'Y = ', num2str(mdl1.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl1.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl1.Rsquared.Ordinary,3)), 'FontSize', 12)
text(xl(2)-20, yl(1)+200, cat(2,'Y = ', num2str(mdl2.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl2.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl2.Rsquared.Ordinary,3)), 'FontSize', 12)

% figure();
% for each slice
mdl_slc = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz = nonzeros(mean_ff_array(:,i));
    mean_t1_array_nz = nonzeros(mean_t1_array(:,i));
    mean_ff_array_nz_new = mean_ff_array_nz;
    mean_t1_array_nz_new = mean_t1_array_nz;
    mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote(:,i));
    mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote(:,i));
    mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
    mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
    mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
    mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];

    mdl_slc = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);
    
    subplot(3,2,i);
    scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_ff_array_nz_new .* mdl_slc.Coefficients.Estimate(2) + mdl_slc.Coefficients.Estimate(1);
    hold on;
    plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    text(xl(2)-30, yl(1)+100, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
end

%% Plot r2star vs t1
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_r2star_array_nz = nonzeros(mean_r2star_array);
mean_t1_array_nz = nonzeros(mean_t1_array);
mean_r2star_array_nz_new = mean_r2star_array_nz;
mean_t1_array_nz_new = mean_t1_array_nz;
mean_r2star_array_nz_new(mean_r2star_array_nz<=0) = [];
mean_t1_array_nz_new(mean_r2star_array_nz<=0) = [];


mdl = fitlm(mean_r2star_array_nz_new, mean_t1_array_nz_new);

figure('Position', [100 100 600 600]);
title(cat(2, name, ' ', time_point))
subplot(3,1,3);
scatter(mean_r2star_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
xlabel('R2star (s^{-1})');
ylabel('T1 (ms)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(2)-20, yl(1)+200, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

% figure();
% for each slice
mdl_slc = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_r2star_array_nz = nonzeros(mean_r2star_array(:,i));
    mean_t1_array_nz = nonzeros(mean_t1_array(:,i));
    mean_r2star_array_nz_new = mean_r2star_array_nz;
    mean_t1_array_nz_new = mean_t1_array_nz;
    mean_r2star_array_nz_new(mean_r2star_array_nz<=0) = [];
    mean_t1_array_nz_new(mean_r2star_array_nz<=0) = [];
    
    mdl_slc = fitlm(mean_r2star_array_nz_new, mean_t1_array_nz_new);
    
    subplot(3,2,i);
    scatter(mean_r2star_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    text(xl(2)-30, yl(1)+100, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('R2star (s^{-1})');
    ylabel('T1 (ms)');
end

%% Plot r2star vs FF (For fat analysis)
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_r2star_array_nz = nonzeros(mean_r2star_array);
mean_ff_array_nz = nonzeros(mean_ff_array);
mean_r2star_array_nz_new = mean_r2star_array_nz;
mean_ff_array_nz_new = mean_ff_array_nz;
mean_r2star_array_nz_new(mean_ff_array_nz<=0) = [];
mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];


mdl = fitlm(mean_r2star_array_nz_new, mean_ff_array_nz_new);

figure('Position', [100 100 600 600]);
title(cat(2, name, ' ', time_point))
subplot(3,1,3);
scatter(mean_r2star_array_nz_new, mean_ff_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
xlabel('R2star (s^{-1})');
ylabel('FF (%)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(2)-60, yl(1)+20, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

% figure();
% for each slice
mdl_slc = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_r2star_array_nz = nonzeros(mean_r2star_array(:,i));
    mean_ff_array_nz = nonzeros(mean_ff_array(:,i));
    mean_r2star_array_nz_new = mean_r2star_array_nz;
    mean_ff_array_nz_new = mean_ff_array_nz;
    mean_r2star_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
    
    mdl_slc = fitlm(mean_r2star_array_nz_new, mean_ff_array_nz_new);
    
    subplot(3,2,i);
    scatter(mean_r2star_array_nz_new, mean_ff_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    text(xl(2)-100, yl(1)+10, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('R2star (s^{-1})');
    ylabel('FF (%)');
end

%% Plot endo

mean_ff_array_nz_endo = nonzeros(mean_ff_array_endo);
mean_t1_array_nz_endo = nonzeros(mean_t1_array_endo);
mean_ff_array_nz_new_endo = mean_ff_array_nz_endo;
mean_t1_array_nz_new_endo = mean_t1_array_nz_endo;
mean_ff_array_nz_new_endo(mean_ff_array_nz_endo<=0) = [];
mean_t1_array_nz_new_endo(mean_ff_array_nz_endo<=0) = [];

mdl_endo = fitlm(mean_ff_array_nz_new_endo, mean_t1_array_nz_new_endo);

figure('Position', [100 100 600 600]);
title(cat(2, name, ' ', time_point))
subplot(3,1,3);
scatter(mean_ff_array_nz_new_endo, mean_t1_array_nz_new_endo, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(2)-20, yl(1)+200, cat(2,'Y = ', num2str(mdl_endo.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_endo.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_endo.Rsquared.Ordinary,3)), 'FontSize', 12)

mdl_slc_endo = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz_endo = nonzeros(mean_ff_array_endo(:,i));
    mean_t1_array_nz_endo = nonzeros(mean_t1_array_endo(:,i));
    mean_ff_array_nz_new_endo = mean_ff_array_nz_endo;
    mean_t1_array_nz_new_endo = mean_t1_array_nz_endo;
    mean_ff_array_nz_new_endo(mean_ff_array_nz_endo<=0) = [];
    mean_t1_array_nz_new_endo(mean_ff_array_nz_endo<=0) = [];
    
    mdl_slc_endo = fitlm(mean_ff_array_nz_new_endo, mean_t1_array_nz_new_endo);
    
    subplot(3,2,i);
    scatter(mean_ff_array_nz_new_endo, mean_t1_array_nz_new_endo, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    text(xl(2)-30, yl(1)+100, cat(2,'Y = ', num2str(mdl_slc_endo.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc_endo.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc_endo.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
end
%% Plot epi

mean_ff_array_nz_epi = nonzeros(mean_ff_array_epi);
mean_t1_array_nz_epi = nonzeros(mean_t1_array_epi);
mean_ff_array_nz_new_epi = mean_ff_array_nz_epi;
mean_t1_array_nz_new_epi = mean_t1_array_nz_epi;
mean_ff_array_nz_new_epi(mean_ff_array_nz_epi<=0) = [];
mean_t1_array_nz_new_epi(mean_ff_array_nz_epi<=0) = [];

mdl_epi = fitlm(mean_ff_array_nz_new_epi, mean_t1_array_nz_new_epi);

figure('Position', [100 100 600 600]);
title(cat(2, name, ' ', time_point))
subplot(3,1,3);
scatter(mean_ff_array_nz_new_epi, mean_t1_array_nz_new_epi, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
% plot(mdl);
yl = ylim;
xl = xlim;
text(xl(2)-20, yl(1)+200, cat(2,'Y = ', num2str(mdl_epi.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_epi.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_epi.Rsquared.Ordinary,3)), 'FontSize', 12)
T1FP_Stats_Analysis
mdl_slc_epi = struct;
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz_epi = nonzeros(mean_ff_array_epi(:,i));
    mean_t1_array_nz_epi = nonzeros(mean_t1_array_epi(:,i));
    mean_ff_array_nz_new_epi = mean_ff_array_nz_epi;
    mean_t1_array_nz_new_epi = mean_t1_array_nz_epi;
    mean_ff_array_nz_new_epi(mean_ff_array_nz_epi<=0) = [];
    mean_t1_array_nz_new_epi(mean_ff_array_nz_epi<=0) = [];
    
    mdl_slc_epi = fitlm(mean_ff_array_nz_new_epi, mean_t1_array_nz_new_epi);
    
    subplot(3,2,i);
    scatter(mean_ff_array_nz_new_epi, mean_t1_array_nz_new_epi, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    text(xl(2)-30, yl(1)+100, cat(2,'Y = ', num2str(mdl_slc_epi.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc_epi.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc_epi.Rsquared.Ordinary,3)), 'FontSize', 10)

    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
end
