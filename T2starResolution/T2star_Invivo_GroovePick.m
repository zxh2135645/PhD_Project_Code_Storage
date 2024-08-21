clear all;
close all;
% This needs to be run after T2star_Invivo_ReadCVIContours.m
% The output is save in '/ContourData_Invivo/<name>/<name>_<time_label>/Results/'

addpath('../function/');
addpath('../T1NFF/');
addpath('../FatFractionMap/');

base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/ContourData_Invivo/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

time_label = {'8WK'};

names_to_rule_out = {'20P40_1Month'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
overwrite_label = 1;

%%
sel_array = [501, 153, 242001, 170, 257001, 172, 270, 267, 264, 202, 201];

% for n = 1:length(Names)
    for n = 8:8
    name = Names{n}
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);

        dst_dir = cat(2, base_dir, '/ContourData_Invivo/', name,'/', name, '_', time_label{t}, '/');

        if exist(dst_dir)
            time = time_label{t}
            save_dir = cat(2, dst_dir, '/Results/');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            for i = 1:length(sequence_label)
                label = sequence_label{i};
                load(cat(2, dst_dir, label, '_vol_img_3D.mat'));
                load(cat(2, dst_dir, 'freeROI/freeROI.mat'));
                load(cat(2, dst_dir, 'MyoReference/myoRef.mat'));
                load(cat(2, dst_dir, 'Myocardium/mask_myocardium.mat'));
                load(cat(2, dst_dir, 'noReflowArea/noReflow.mat'));
                load(cat(2, dst_dir, 'Heart/mask_heart.mat'));
                load(cat(2, dst_dir, label, '_SliceLoc.mat'));

                [slc_array_sorted, idx_sorted] = sort(slc_array);

                mask_myocardium_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                vol_img_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                noReflowMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                freeROIMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                mask_heart_tight = zeros(size(mask_heart_3D{1},1), size(mask_heart_3D{1},2), length(idx_sorted));

                for ii = 1:length(mask_myocardium_3D)
                    mask_myocardium_tight(:,:,ii) = mask_myocardium_3D{idx_sorted(ii)};
                    vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                    noReflowMask_tight(:,:,ii) = noReflowMask_3D{idx_sorted(ii)};
                    freeROIMask_tight(:,:,ii) = freeROIMask_3D{idx_sorted(ii)};
                    mask_heart_tight(:,:,ii) = mask_heart_3D{idx_sorted(ii)};
                end

                noReflowMask_tight = mask_myocardium_tight .* noReflowMask_tight;
                freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;

                mvo_perc = zeros(length(slc_array), 1);
                infarct_perc = zeros(length(slc_array), 1);

                for slc = 1:size(mask_myocardium_tight, 3)
                    mvo_perc(slc) = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                    infarct_perc(slc) = numel(nonzeros(freeROIMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                end
                % infarct_perc*100
                
                outputFileName = [save_dir, label, '_coords.mat'];
                [x, y, x_centroid, y_centroid] = GroovePickNCheckFunc3(vol_img_tight, mask_heart_tight, outputFileName, overwrite_label);

                coords.x = x; % Rows
                coords.y = y; % Column
                coords.x_centroid = x_centroid;
                coords.y_centroid = y_centroid;

                save(outputFileName, 'coords');
            end
        end
    end
end

%% For VaryingRes
name_glob = glob(cat(2, base_dir, '/ContourData_Invivo/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

time_label = {'8WK_VaryingRes'};

names_to_rule_out = {'17P73', '18P90', '18P92', '18P93', '18P94', '18P95', '20P10', '20P11', '20P40'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
overwrite_label = 1;

count = 0;
sel_cell = cell(size(Names, 1), 1); 
sel_cell{1} = [204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240];
sel_cell{2} = [216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248];


for n = 1:length(Names)
    % for n = 9:9
    name = Names{n}
    sel_array = sel_cell{n};
    coords = struct;
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);

        dst_dir = cat(2, base_dir, '/ContourData_Invivo/', name,'/', name, '_', time_label{t}, '/');

        if exist(dst_dir)
            time = time_label{t}
            save_dir = cat(2, dst_dir, '/Results/');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            for i = 1:length(sequence_label)
                label = sequence_label{i};
                for j = 1:length(sel_array)
                    sel_char = cat(2, 'series',num2str(sel_array(j), '%04.f'), '-Body/');
                    load(cat(2, dst_dir, sel_char, label, '_vol_img_3D.mat'));
                    load(cat(2, dst_dir, sel_char, 'freeROI/freeROI.mat'));
                    load(cat(2, dst_dir, sel_char, 'MyoReference/myoRef.mat'));
                    load(cat(2, dst_dir, sel_char, 'Myocardium/mask_myocardium.mat'));
                    load(cat(2, dst_dir, sel_char, 'noReflowArea/noReflow.mat'));
                    load(cat(2, dst_dir, sel_char, 'Heart/mask_heart.mat'));
                    load(cat(2, dst_dir, sel_char, label, '_SliceLoc.mat'));

                    [slc_array_sorted, idx_sorted] = sort(slc_array);

                    mask_myocardium_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    vol_img_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    noReflowMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    freeROIMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    mask_heart_tight = zeros(size(mask_heart_3D{1},1), size(mask_heart_3D{1},2), length(idx_sorted));

                    for ii = 1:length(mask_myocardium_3D)
                        mask_myocardium_tight(:,:,ii) = mask_myocardium_3D{idx_sorted(ii)};
                        vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                        noReflowMask_tight(:,:,ii) = noReflowMask_3D{idx_sorted(ii)};
                        freeROIMask_tight(:,:,ii) = freeROIMask_3D{idx_sorted(ii)};
                        mask_heart_tight(:,:,ii) = mask_heart_3D{idx_sorted(ii)};
                    end

                    noReflowMask_tight = mask_myocardium_tight .* noReflowMask_tight;
                    freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;

                    mvo_perc = zeros(length(slc_array), 1);
                    infarct_perc = zeros(length(slc_array), 1);

                    for slc = 1:size(mask_myocardium_tight, 3)
                        mvo_perc(slc) = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                        infarct_perc(slc) = numel(nonzeros(freeROIMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                    end
                    % infarct_perc*100

                    [x, y, x_centroid, y_centroid] = GroovePickNCheckFunc3(vol_img_tight, mask_heart_tight);

                    coords(j).x = x; % Rows
                    coords(j).y = y; % Column
                    coords(j).x_centroid = x_centroid;
                    coords(j).y_centroid = y_centroid;
                    coords(j).seriesNum = sel_array(j);
                end
            end
            
            outputFileName = [save_dir, label, '_coords.mat'];
            save(outputFileName, 'coords');
        end
    end
end

%% 0.3x0.3x2 ground truth (Read DICOMs)
subject_name_cell = {'DHARMAKUMAR_18P90_EXVIVO2_FORMALIN', 'DHARMAKUMAR_18P93_EXVIVO2_FORMALIN', ...
    '20P40_EXVIVO2_JAMES', 'DHARMAKUMAR_20P10_EXVIVO7_FORMALIN', 'DHARMAKUMAR_20P11_EXVIVO6_FORMALIN',...
    'DHARMAKUMAR_18P92_EXVIVO2_FORMALIN', 'DHARMAKUMAR_18P94_EXVIVO3_FORMALIN', 'DHARMAKUMAR_18P95_EXVIVO2_FORMALIN',...
    'DHARMAKUMAR_17P73_EXVIVO2_FORMALIN', 'JAMES_PHANTOM_20P48_EXVIVO'};
sel_array = [70, 73, 69, 79, 71, 65, 94, 65, 79, 85]
dicom_dir = uigetdir;


whatsinit = cell(length(sel_array), 1);
slice_data = cell(length(sel_array), 1);


for i = 1:length(sel_array)
    seriesNum = num2str(sel_array(i), '%04.f');
    dicom_glob = glob(cat(2, dicom_dir, '/', subject_name_cell{i}, '/*/T2STAR_IMAGES_', seriesNum));
    list_to_read = dicom_glob{1};
    [whatsinit{i}, slice_data{i}] = dicom23D(list_to_read);
end

%% Draw insertion points
subject_name_cell2 = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
coords = struct;
save_dir = cat(2, base_dir, '/Data/Groove_GT/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%for n = 1:length(subject_name_cell2)
    for n = 4:4
    name = subject_name_cell2{n}
    %mvo_perc = zeros(length(sequence_label), 1);
    %infarct_perc = zeros(length(sequence_label), 1);

    dst_dir = cat(2, base_dir, '/Data/', name, '/');

    if exist(dst_dir)
        load(cat(2, dst_dir, 'mask.mat'));

        fluid_mask = mask_struct(1).fluid_mask;
        img_gt = whatsinit{n};

        [x, y, x_centroid, y_centroid] = GroovePickNCheckFunc3(img_gt, fluid_mask);

        coords(n).x = x; % Rows
        coords(n).y = y; % Column
        coords(n).x_centroid = x_centroid;
        coords(n).y_centroid = y_centroid;
        coords(n).name = name;
    end
end

outputFileName = cat(2, save_dir, 'GroundTruth_');
save(outputFileName, 'coords');

% Needs to compare invivo and exvivo