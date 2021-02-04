clear all;
close all;
addpath('../function/');
addpath('../AHA16Segment/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
time_points = {'8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};

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

%%
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
    
    %for tp = 1:length(time_points)
    for tp = 2:2
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
            % if condition
            center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');
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
            
            % AHA Segment
            Segn = 50;
            Groove = 0;
            MI_Chord_Analysis_fname = cat(2, name_data_save_dir, '/MIChordAnalysis_', name, '_', time_point, '.mat');
            
            %[MI_Chord_Analysis, center_mask_ff] = Func_MI_Chord_Analysis(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff,...
            %    roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, MI_Chord_Analysis_fname);
            [MI_Chord_Analysis2, ~] = Func_MI_Chord_Analysis2(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff,...
                roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, MI_Chord_Analysis_fname);
            % Plot figures and save
            %Func_plot_chord_analysis_general(MI_Chord_Analysis, tp_dir2);
            %Func_plot_chord_analysis_EpiEndo(MI_Chord_Analysis, tp_dir2);
            Func_plot_chord_analysis_general2(MI_Chord_Analysis2, tp_dir2);
            Func_plot_chord_analysis_EpiEndo2(MI_Chord_Analysis2, tp_dir2);
        end
    end
    close all;
end

%% Overlap (register) ROI 
roi_intercpt = roi_in_myo_ff & roi_in_myo_t1;
% figure();
% for i = 1:size(roi_intercpt, 3)
%    subplot(2,2,i);
%    imagesc(roi_in_myo_ff(:,:,i)*(-1) + roi_in_myo_t1(:,:,i)*2);
% end

figure();
for i = 1:size(roi_intercpt, 3)
    subplot(2,2,i)
    fixed = roi_in_myo_t1(:,:,i);
    moving = roi_in_myo_ff(:,:,i);
    imshowpair(fixed, moving,'Scaling','joint')
end
[optimizer, metric] = imregconfig('multimodal');

figure();
for i = 1:size(roi_intercpt, 3)
    subplot(2,2,i)
    fixed = roi_in_myo_t1(:,:,i);
    moving = roi_in_myo_ff(:,:,i);
    %movingRegistered = imregister(moving, fixed, 'rigid', optimizer, metric);
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
    imshowpair(fixed, movingRegistered,'Scaling','joint');
end

%% Overlap (register) Myocardium
Segn = 50;
Groove = 0;
roi_intercpt = roi_in_myo_ff & roi_in_myo_t1;
caxis_rg = [0 1];

center_mask_t1 = zeros(size(roi_in_myo_t1));

for i = 1:size(roi_intercpt, 3)
    
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    roi_fixed = roi_in_myo_t1(:,:,i);
    roi_moving = roi_in_myo_ff(:,:,i);

    
    subplot(2,4,2*(i-1)+1);
    imshowpair(fixed, moving,'Scaling','joint');
    title(cat(2, 'slice ', num2str(i)));
    
    subplot(2,4,2*(i-1)+2);
    imshowpair(roi_fixed, roi_moving,'Scaling','joint');
    title(cat(2, 'slice ', num2str(i)));
end
%% 
center_mask_ff = zeros(size(roi_in_myo_ff));

for i = 1:size(roi_intercpt, 3)
    moving = myo_ff(:,:,i);
    center_mask_ff(:,:,i) = Func_DrawCenterLine(moving, caxis_rg);
end
%%
[optimizer, metric] = imregconfig('multimodal');
movingRegistered_center = zeros(size(roi_in_myo_t1));
figure();
for i = 1:size(roi_intercpt, 3)
    
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    roi_fixed = roi_in_myo_t1(:,:,i);
    roi_moving = roi_in_myo_ff(:,:,i);
    
    center_moving = center_mask_ff(:,:,i);
    
    % [movingRegistered, R_reg] = imregister(moving, fixed, 'rigid', optimizer, metric);
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
    movingRegistered_roi = imwarp(roi_moving,tform,'OutputView',imref2d(size(fixed)));
    movingRegistered_center(:,:,i) = imwarp(center_moving,tform,'OutputView',imref2d(size(fixed)));
    
    subplot(2,4,2*(i-1)+1);
    imshowpair(fixed, movingRegistered,'Scaling','joint');
    title(cat(2, 'slice ', num2str(i)));
    
    subplot(2,4,2*(i-1)+2);
    imshowpair(roi_fixed, movingRegistered_roi,'Scaling','joint');
    title(cat(2, 'slice ', num2str(i)));
end 



%%
Mipix = cell(Segn, size(t1, 3));
Mipix_epi = cell(Segn, size(t1, 3));
Mipix_endo = cell(Segn, size(t1, 3));
Mipix2 = cell(Segn, size(t1, 3));
Mipix2_epi = cell(Segn, size(t1, 3));
Mipix2_endo = cell(Segn, size(t1, 3));
Mipix3 = cell(Segn, size(t1, 3));
Mipix3_epi = cell(Segn, size(t1, 3));
Mipix3_endo = cell(Segn, size(t1, 3));

Mipix_mean = [];
Mipix_mean2 = [];
Mipix_mean3 = [];
Mipix_mean_epi = [];
Mipix_mean2_epi = [];
Mipix_mean3_epi = [];
Mipix_mean_endo = [];
Mipix_mean2_endo = [];
Mipix_mean3_endo = [];

se = strel('disk', 1);

for i = 1:size(t1, 3)
%for i = 4:4
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    
    center_moving = center_mask_ff(:,:,i);
    myo_epi_moving = moving - center_moving > 0;
    myo_endo_moving = center_moving + moving > 1;
    
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
    img2 = imwarp(img2,tform,'OutputView',imref2d(size(fixed)));
    
    movingRegistered_epi = imwarp(myo_epi_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    movingRegistered_endo = imwarp(myo_endo_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    
    movingRegistered_roi_ff = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;

    fixed_eroded = imerode(fixed, se);
    movingRegistered_eroded = imerode(movingRegistered, se);
    
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(img, fixed, Segn, Groove);
    [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(img, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(img, movingRegistered_endo, Segn, Groove);


    [Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, movingRegistered, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_epi] = AHASegmentation(img2, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_endo] = AHASegmentation(img2, movingRegistered_endo, Segn, Groove);
    
    img3 = r2star(:,:,i);
    img3 = imwarp(img3,tform,'OutputView',imref2d(size(fixed)));
    [Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, movingRegistered, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_epi] = AHASegmentation(img3, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_endo] = AHASegmentation(img3, movingRegistered_endo, Segn, Groove);
    
    
    
    for j = 1:Segn
        Mipix{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        Mipix_epi{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        Mipix_endo{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        
        Mipix2{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        Mipix2_epi{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        Mipix2_endo{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        
        Mipix3{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
        Mipix3_epi{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
        Mipix3_endo{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
        
    end

    Mi_idx1 = [];
    Mi_idx1_epi = [];
    Mi_idx1_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix{j,i})
            Mi_idx1 = [Mi_idx1, j];
        end
        if ~isempty(Mipix_epi{j,i})
            Mi_idx1_epi = [Mi_idx1_epi, j];
        end
        if ~isempty(Mipix_endo{j,i})
            Mi_idx1_endo = [Mi_idx1_endo, j];
        end
    end
    
    Mi_idx2 = [];
    Mi_idx2_epi = [];
    Mi_idx2_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix2{j,i})
            Mi_idx2 = [Mi_idx2, j];
        end
        if ~isempty(Mipix2_epi{j,i})
            Mi_idx2_epi = [Mi_idx2_epi, j];
        end
        if ~isempty(Mipix_endo{j,i})
            Mi_idx2_endo = [Mi_idx2_endo, j];
        end
    end
    
    Mi_idx3 = [];
    Mi_idx3_epi = [];
    Mi_idx3_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix3{j,i})
            Mi_idx3 = [Mi_idx3, j];
        end
        if ~isempty(Mipix3_epi{j,i})
            Mi_idx3_epi = [Mi_idx3_epi, j];
        end
        if ~isempty(Mipix3_endo{j,i})
            Mi_idx3_endo = [Mi_idx3_endo, j];
        end
    end

    % 
    Mi_idx_intersect = intersect(Mi_idx1, Mi_idx2);
    Mi_idx_intersect_epi = intersect(Mi_idx1_epi, Mi_idx2_epi);
    Mi_idx_intersect_endo = intersect(Mi_idx1_endo, Mi_idx2_endo);
    Mipix_flat = {};
    Mipix_flat2 = {};
    Mipix_flat3 = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect)
        ind = Mi_idx_intersect(j);
        Mipix_flat{idx} = Mipix{ind,i};
        Mipix_flat2{idx} = Mipix2{ind,i};
        Mipix_flat2{idx}(Mipix_flat2{idx} > 100) = 100;
        Mipix_flat2{idx}(Mipix_flat2{idx} < 0) = 0;
        
        Mipix_flat3{idx} = Mipix3{ind,i};
        Mipix_flat3{idx}(Mipix_flat3{idx} < 0) = 0;
        
        Mipix_mean(idx, i) = mean(Mipix_flat{idx});
        Mipix_mean2(idx, i) = mean(Mipix_flat2{idx});
        Mipix_mean3(idx, i) = mean(Mipix_flat3{idx});
        idx = idx + 1;
    end
    
    % epi
    Mipix_flat_epi = {};
    Mipix_flat2_epi = {};
    Mipix_flat3_epi = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect_epi)
        ind = Mi_idx_intersect_epi(j);
        Mipix_flat_epi{idx} = Mipix_epi{ind,i};
        Mipix_flat2_epi{idx} = Mipix2_epi{ind,i};
        Mipix_flat2_epi{idx}(Mipix_flat2_epi{idx} > 100) = 100;
        Mipix_flat2_epi{idx}(Mipix_flat2_epi{idx} < 0) = 0;
        
        Mipix_flat3_epi{idx} = Mipix3_epi{ind,i};
        Mipix_flat3_epi{idx}(Mipix_flat3_epi{idx} < 0) = 0;
        
        Mipix_mean_epi(idx, i) = mean(Mipix_flat_epi{idx});
        Mipix_mean2_epi(idx, i) = mean(Mipix_flat2_epi{idx});
        Mipix_mean3_epi(idx, i) = mean(Mipix_flat3_epi{idx});
        idx = idx + 1;
    end
    
    % endo
    Mipix_flat_endo = {};
    Mipix_flat2_endo = {};
    Mipix_flat3_endo = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect_endo)
        ind = Mi_idx_intersect_endo(j);
        Mipix_flat_endo{idx} = Mipix_endo{ind,i};
        Mipix_flat2_endo{idx} = Mipix2_endo{ind,i};
        Mipix_flat2_endo{idx}(Mipix_flat2_endo{idx} > 100) = 100;
        Mipix_flat2_endo{idx}(Mipix_flat2_endo{idx} < 0) = 0;
        
        Mipix_flat3_endo{idx} = Mipix3_endo{ind,i};
        Mipix_flat3_endo{idx}(Mipix_flat3_endo{idx} < 0) = 0;
        
        Mipix_mean_endo(idx, i) = mean(Mipix_flat_endo{idx});
        Mipix_mean2_endo(idx, i) = mean(Mipix_flat2_endo{idx});
        Mipix_mean3_endo(idx, i) = mean(Mipix_flat3_endo{idx});
        idx = idx + 1;
    end
end

%% plot
[I,J] = find(Mipix_mean);
Mipix_mean_nnz = [];
Mipix_mean2_nnz = [];
Mipix_mean3_nnz = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz(idx) = Mipix_mean(I(i), J(i));
    Mipix_mean2_nnz(idx) = Mipix_mean2(I(i), J(i));
    Mipix_mean3_nnz(idx) = Mipix_mean3(I(i), J(i));
    idx = idx + 1;
end

figure();
subplot(3,1,1);
imagesc(Mipix_mean_nnz);
colorbar; title('T1 (ms)')
caxis([1000 1600]);
subplot(3,1,2);
imagesc(Mipix_mean2_nnz);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,1,3);
imagesc(Mipix_mean3_nnz);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));
%% plot II
[I,J] = find(Mipix_mean_epi);
Mipix_mean_nnz_epi = [];
Mipix_mean2_nnz_epi = [];
Mipix_mean3_nnz_epi = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz_epi(idx) = Mipix_mean_epi(I(i), J(i));
    Mipix_mean2_nnz_epi(idx) = Mipix_mean2_epi(I(i), J(i));
    Mipix_mean3_nnz_epi(idx) = Mipix_mean3_epi(I(i), J(i));
    idx = idx + 1;
end

[I,J] = find(Mipix_mean_endo);
Mipix_mean_nnz_endo = [];
Mipix_mean2_nnz_endo = [];
Mipix_mean3_nnz_endo = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz_endo(idx) = Mipix_mean_endo(I(i), J(i));
    Mipix_mean2_nnz_endo(idx) = Mipix_mean2_endo(I(i), J(i));
    Mipix_mean3_nnz_endo(idx) = Mipix_mean3_endo(I(i), J(i));
    idx = idx + 1;
end

figure();
subplot(3,2,1);
imagesc(Mipix_mean_nnz_epi);
caxis([1000 1600]);
colorbar; title('T1 (ms)')
subplot(3,2,3);
imagesc(Mipix_mean2_nnz_epi);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,2,5);
imagesc(Mipix_mean3_nnz_epi);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));

subplot(3,2,2);
imagesc(Mipix_mean_nnz_endo);
caxis([1000 1600]);
colorbar; title('T1 (ms)')
subplot(3,2,4);
imagesc(Mipix_mean2_nnz_endo);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,2,6);
imagesc(Mipix_mean3_nnz_endo);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));
%% scatter
figure();
subplot(2,1,1);
scatter(Mipix_mean2_nnz, Mipix_mean_nnz, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20)
%ylim([800 2000])
grid on;

subplot(2,1,2);
scatter(Mipix_mean3_nnz, Mipix_mean_nnz, 48, 'filled');
xlabel('R2star (Hz)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20)
grid on;

% figure();
% scatter3(Mipix_mean2_nnz, Mipix_mean3_nnz, Mipix_mean_nnz, 'filled');
% xlabel('FF (%)'); ylabel('R2star (Hz)'); zlabel('T1 (ms)');

%% scatter II
figure();
subplot(2,2,1);
scatter(Mipix_mean2_nnz_epi, Mipix_mean_nnz_epi, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([800 2000]);
title('Epi');
grid on;

subplot(2,2,3);
scatter(Mipix_mean3_nnz_epi, Mipix_mean_nnz_epi, 48, 'filled');
xlabel('R2star (Hz)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20)
grid on;

subplot(2,2,2);
scatter(Mipix_mean2_nnz_endo, Mipix_mean_nnz_endo, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([800 2000]);
title('Endo');
grid on;

subplot(2,2,4);
scatter(Mipix_mean3_nnz_endo, Mipix_mean_nnz_endo, 48, 'filled');
xlabel('R2star (Hz)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20)
grid on;
%% Deprecated
roi_intercpt_t1 = roi_intercpt .* t1;
roi_intercpt_ff = roi_intercpt .* ff;


figure();
for i = 1:size(roi_intercpt_t1, 3)
   subplot(2,2,i)
   imagesc(roi_intercpt_t1(:,:,i))
end

figure();
roi_intercpt_ff(roi_intercpt > 100) = 100;
roi_intercpt_ff(roi_intercpt < 0) = 0;

for i = 1:size(roi_intercpt_ff, 3)
   subplot(2,2,i)
   imagesc(roi_intercpt_ff(:,:,i))
end

figure();
roi_nofat_ff = roi_intercpt_ff < 10;
for i = 1:size(roi_intercpt_ff, 3)
   subplot(2,2,i)
   imagesc(roi_nofat_ff(:,:,i))
end

roi_nofat_t1 = roi_intercpt_t1 .* roi_nofat_ff;
figure();
for i = 1:size(roi_intercpt_ff, 3)
   subplot(2,2,i)
   imagesc(roi_nofat_t1(:,:,i))
end