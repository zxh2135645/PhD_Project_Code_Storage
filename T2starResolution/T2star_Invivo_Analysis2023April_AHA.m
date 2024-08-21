clear all;
close all;

addpath('../function/');
addpath('../T1NFF/');
addpath('../FatFractionMap/');
addpath('../AHA16Segment/');

base_dir = uigetdir;

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

name_glob = glob(cat(2, base_dir, '/ContourData_Invivo/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

time_label = {'8WK'};

names_to_rule_out = {'20P10', '20P40_1Month'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
label = labels{5};
%% 0.3x0.3x2 ground truth (Read DICOMs)
subject_name_cell = {'DHARMAKUMAR_18P90_EXVIVO2_FORMALIN', 'DHARMAKUMAR_18P93_EXVIVO2_FORMALIN', ...
    '20P40_EXVIVO2_JAMES', 'DHARMAKUMAR_20P10_EXVIVO7_FORMALIN', 'DHARMAKUMAR_20P11_EXVIVO6_FORMALIN',...
    'DHARMAKUMAR_18P92_EXVIVO2_FORMALIN', 'DHARMAKUMAR_18P94_EXVIVO3_FORMALIN', 'DHARMAKUMAR_18P95_EXVIVO2_FORMALIN',...
    'DHARMAKUMAR_17P73_EXVIVO2_FORMALIN', 'JAMES_PHANTOM_20P48_EXVIVO'};
sel_array = [70, 73, 69, 79, 71, 65, 94, 65, 79, 85]
dicom_dir = uigetdir;

whatsinit_gt = cell(length(sel_array), 1);
slice_data_gt = cell(length(sel_array), 1);


for i = 1:length(sel_array)
    seriesNum = num2str(sel_array(i), '%04.f');
    dicom_glob = glob(cat(2, dicom_dir, '/', subject_name_cell{i}, '/*/T2STAR_IMAGES_', seriesNum));
    list_to_read = dicom_glob{1};
    [whatsinit_gt{i}, slice_data_gt{i}] = dicom23D(list_to_read);
end

%% Read dicom INVIVO
dicom_dir = uigetdir;
dicom_glob = glob(cat(2, dicom_dir, '/*'));
sel_array = [501, 153, 242001, 170, 257001, 172, 270, 264, 202, 201];
sel_array = [501, 153, 242001, 170, 257001, 172, 267, 264, 202, 201];

whatsinit = cell(length(sel_array), 1);
slice_data = cell(length(sel_array), 1);
count = 0;

for i = 1:length(dicom_glob)
    if i ~= 11 && i ~= 7
        count = count + 1;
        idx_array = contains(dicom_glob, label);

        ya_glob = glob(cat(2, dicom_glob{i}, 'DICOM/*'));
        dst_names = ExtractNames(ya_glob);
        sub_name = num2str(sel_array(count), '%04.f');
        
        ind_array = zeros(size(dst_names, 1), 1);
        ind_array = ind_array + contains(dst_names, sub_name);
        ind_array2 = find(ind_array == 1);
        list_to_read = ya_glob(ind_array2);

        for j = 1:length(list_to_read)
            % [whatsinit{count}, slice_data{count}] = dicom23D(list_to_read{j});
            yya_glob = glob(cat(2, list_to_read{j}, '*'));
            whatsinit{count} = dicomread(yya_glob{1});
            slice_data{count} = dicominfo(yya_glob{1});
        end
    end
end

rescaled_t2star = cell(length(sel_array), 1);
for i = 1:length(sel_array)
    rescaled_t2star{i} = double(whatsinit{i} .* slice_data{i}.RescaleSlope);
end
%%
subject_name_cell = {'17P73', '18P90', '18P92', '18P93', '18P94_Exvivo3', '18P95', '20P10_Exvivo7', '20P11_Exvivo6', '20P40', '20P48'};
Segn = 50;
aha_mi = struct;
gt_groove_dir = cat(2, base_dir, '/data/Groove_GT/');
coords_gt = load(cat(2, gt_groove_dir, 'GroundTruth_.mat'));
idx_array = zeros(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    for j = 1:length(subject_name_cell)
        flg = strcmp(coords_gt.coords(j).name, subject_name_cell{i});
        if flg
            idx_array(i) = j;
        end
    end
end

perc_array_gt_cell = cell(length(Names),1);
perc_array_mi2_gt_cell = cell(length(Names),1);
perc_array_cell = cell(length(Names),1);
perc_array_mi2_cell = cell(length(Names),1);
res_mi_cell = cell(length(Names),1);
res_mi_gt_cell = cell(length(Names),1);

for n = 1:length(Names)
    %for n = 7:7
    name = Names{n}
    name_gt = subject_name_cell{n};
    for t = 1:length(time_label)

        dst_dir = cat(2, base_dir, '/ContourData_Invivo/', name,'/', name, '_', time_label{t}, '/');
        gt_dir = cat(2, base_dir, '/data/', name_gt, '/');
        save_dir = cat(2, dst_dir, '/Results/');

        for i = 1:length(sequence_label)
            label = sequence_label{i};
            load(cat(2, dst_dir, label, '_vol_img_3D.mat'));
            load(cat(2, dst_dir, 'Myocardium/mask_myocardium.mat'));
            load(cat(2, dst_dir, 'excludeArea/excludeArea.mat'));
            load(cat(2, dst_dir, 'MyoReference/myoRef.mat'));
            load(cat(2, dst_dir, label, '_SliceLoc.mat'));
            load(cat(2, save_dir, label, '_coords.mat'));
            load(cat(2, gt_dir, 'mask.mat'));
            
            mi_mask = mask_struct(1).mi_mask;
            remote_mask = mask_struct(1).remote_mask;
            myo_mask = mask_struct(1).myo_mask;
            fluid_mask = mask_struct(1).fluid_mask;

            [slc_array_sorted, idx_sorted] = sort(slc_array);

            mask_myocardium_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
            vol_img_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
            excludeArea_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
            MyoRef_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));

            for ii = 1:length(mask_myocardium_3D)
                mask_myocardium_tight(:,:,ii) = mask_myocardium_3D{idx_sorted(ii)};
                vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                excludeArea_tight(:,:,ii) = excludeMask_3D{idx_sorted(ii)};
                MyoRef_tight(:,:,ii) = myoRefMask_3D{idx_sorted(ii)};
            end
            
            MyoRef_tight = mask_myocardium_tight .* MyoRef_tight;

            x_gt = coords_gt.coords(idx_array(n)).x;
            y_gt = coords_gt.coords(idx_array(n)).y;
            x_centroid_gt = coords_gt.coords(idx_array(n)).x_centroid;
            y_centroid_gt = coords_gt.coords(idx_array(n)).y_centroid;
            Groove_gt = -(atan2(x_gt - x_centroid_gt, y_gt - y_centroid_gt) * 180 / pi);

            x = coords.x;
            y = coords.y;
            x_centroid = coords.x_centroid;
            y_centroid = coords.y_centroid;
            Groove = atan2(x - x_centroid, y - y_centroid) * 180 / pi; 
            % TODO Something compensated in AHASegmentation_CounterClock
            
            hemo_thresh_gt = mean(nonzeros(whatsinit_gt{idx_array(n)} .* myo_mask .* remote_mask)) - 2*std(nonzeros(whatsinit_gt{idx_array(n)} .* myo_mask .* remote_mask));
            hemo_mask_gt = (whatsinit_gt{idx_array(n)} < hemo_thresh_gt) .* mi_mask .* myo_mask;

            hemo_thresh = mean(nonzeros(rescaled_t2star{n} .* MyoRef_tight)) - 2*std(nonzeros(rescaled_t2star{n} .* MyoRef_tight));
            hemo_mask = (rescaled_t2star{n} < hemo_thresh) .* (~excludeArea_tight) .* mask_myocardium_tight;
            
            [Segmentpix_gt, stats_gt, Mask_Segn_gt] = AHASegmentation_CounterClock(whatsinit_gt{idx_array(n)}, myo_mask, Segn, Groove_gt, fluid_mask);
            [Segmentpix, stats, Mask_Segn] = AHASegmentation(rescaled_t2star{n}, mask_myocardium_tight, Segn, Groove);
            
            perc_array_gt = [];
            for j = 1:Segn
                perc_array_gt(j) = sum(Segmentpix_gt{j}<hemo_thresh_gt) / length(Segmentpix_gt{j});
            end

            Mipix_gt = cell(Segn, 1);
            for j = 1:Segn
                Mipix_gt{j,1} = whatsinit_gt{idx_array(n)}(Mask_Segn_gt .* mi_mask == j);
            end
            aha_mi_gt(n).Mipix = Mipix_gt;

            Mipix_flat_gt = {};
            idx = 1;
            for j = 1:Segn
                if ~isempty(aha_mi_gt(n).Mipix{j,1})
                    Mipix_flat_gt{idx} = aha_mi_gt(n).Mipix{j,1};
                    idx = idx + 1;
                end
            end
            aha_mi_gt(n).Mipix_flat = Mipix_flat_gt;
            
            perc_array_mi2_gt = [];
            for j = 1:length(Mipix_flat_gt)
                perc_array_mi2_gt(j) = sum(aha_mi_gt(n).Mipix_flat{j}<hemo_thresh_gt) / length(aha_mi_gt(n).Mipix_flat{j});
            end

            res_mi_gt = perc_array_mi2_gt > 0.1;


            perc_array = [];
            for j = 1:Segn
                perc_array(j) = sum(Segmentpix{j}<hemo_thresh) / length(Segmentpix{j});
            end

            Mipix = cell(Segn, 1);
            for j = 1:Segn
                Mipix{j,1} = rescaled_t2star{n}(Mask_Segn .* (~excludeArea_tight) .* mask_myocardium_tight == j);
            end
            aha_mi(n).Mipix = Mipix;

            Mipix_flat = {};
            idx = 1;
            for j = 1:Segn
                if ~isempty(aha_mi(n).Mipix{j,1})
                    Mipix_flat{idx} = aha_mi(n).Mipix{j,1};
                    idx = idx + 1;
                end
            end
            aha_mi(n).Mipix_flat = Mipix_flat;
            
            perc_array_mi2 = [];
            for j = 1:length(Mipix_flat)
                perc_array_mi2(j) = sum(aha_mi(n).Mipix_flat{j}<hemo_thresh) / length(aha_mi(n).Mipix_flat{j});
            end

            res_mi = perc_array_mi2 > 0.1;
            
            perc_array_gt_cell{n} = perc_array_gt;
            perc_array_mi2_gt_cell{n} = perc_array_mi2_gt;
            perc_array_cell{n} = perc_array;
            perc_array_mi2_cell{n} = perc_array_mi2;
            res_mi_cell{n} = res_mi;
            res_mi_gt_cell{n} = res_mi_gt;

            % [cm, order] = confusionmat(res_mi(1,:),res_mi2(i,:));
            % subplot(4,7,i);confusionchart(cm, order);
            % set(gca, 'FontSize', 18);
            % sens_mi2(i) = cm(2,2) / (cm(2,2) + cm(2,1));
            % spec_mi2(i) = cm(1,1) / (cm(1,1) + cm(1,2));
            % accu_mi2(i) = (cm(1,1)+cm(2,2)) ./ (cm(1,1)+cm(1,2)+cm(2,1)+cm(2,2));
        end
    end
end

%%
% {'17P73', '18P90', '18P92', '18P93', '18P94_Exvivo3', '18P95', '20P10_Exvivo7', '20P11_Exvivo6', '20P40', '20P48'};
opt_sens = zeros(length(res_mi_gt_cell), 1);
opt_accu = zeros(length(res_mi_gt_cell), 1);

for i = 1:length(res_mi_gt_cell)
    %for i = 4:4
    sens_mi = [];
    spec_mi = [];
    accu_mi = [];
    res_mi = res_mi_cell{i};
    res_mi_gt = res_mi_gt_cell{i};
    len_mi = length(res_mi);
    len_mi_gt = length(res_mi_gt);

    [min_num] = min(len_mi, len_mi_gt);
    [max_num] = max(len_mi, len_mi_gt);
    diff_num = max_num - min_num + 1;

    if len_mi == min_num
        for j = 1:diff_num
            res_mi_gt_trunc = res_mi_gt(j:j+min_num-1);
            [cm, order] = confusionmat(res_mi_gt_trunc, res_mi);
            sens_mi(j) = cm(2,2) / (cm(2,2) + cm(2,1));
            spec_mi(j) = cm(1,1) / (cm(1,1) + cm(1,2));
            accu_mi(j) = (cm(1,1)+cm(2,2)) ./ (cm(1,1)+cm(1,2)+cm(2,1)+cm(2,2));
        end
    else
        for j = 1:diff_num
            res_mi_trunc = res_mi(j:j+min_num-1);
            [cm, order] = confusionmat(res_mi_gt, res_mi_trunc);
            sens_mi(j) = cm(2,2) / (cm(2,2) + cm(2,1));
            spec_mi(j) = cm(1,1) / (cm(1,1) + cm(1,2));
            accu_mi(j) = (cm(1,1)+cm(2,2)) ./ (cm(1,1)+cm(1,2)+cm(2,1)+cm(2,2));
        end
    end
    [opt_sens(i), I1] = max(sens_mi);
    % [m_accu, I2] = max(accu_mi)
    opt_accu(i) = accu_mi(I1);
end

% Is there a better strategy? Shall I abort?

% Sensitivity
[0.500000000000000
0.333333333333333
0
0
1
0
0.545454545454545
0
0
0.222222222222222];

% Accuracy
[0.750000000000000
0.764705882352941
0.600000000000000
0.833333333333333
1
0.714285714285714
0.583333333333333
0.600000000000000
0.812500000000000
0.608695652173913];