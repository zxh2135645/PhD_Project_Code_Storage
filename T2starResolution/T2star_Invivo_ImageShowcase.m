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

names_to_rule_out = {'20P10C', '20P40_1Month'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
label = labels{5};

dicom_dir = uigetdir;
%% Read dicom INVIVO (T2star mapping)
dicom_glob = glob(cat(2, dicom_dir, '/*'));
sel_array = [501, 153, 242001, 170, 257001, 172, 270, 264, 202, 201];
%sel_array = [501, 153, 242001, 170, 257001, 172, 267, 264, 202, 201];

whatsinit = cell(length(sel_array), 1);
slice_data = cell(length(sel_array), 1);
count = 0;

for i = 1:length(dicom_glob)
    if i ~= 11 && i ~= 8
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

%% Read dicom INVIVO (T2star weighted)
dicom_glob = glob(cat(2, dicom_dir, '/*'));
sel_array = [147, 152, 241001, 169, 256001, 171, 269, 263, 201, 200];
%sel_array = [147, 152, 241001, 169, 256001, 171, 266, 263, 201, 200];

whatsinit_w = cell(length(sel_array), 1);
slice_data_w = cell(length(sel_array), 1);
count = 0;

for i = 1:length(dicom_glob)
    if i ~= 11 && i ~= 8
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
            whatsinit_w{count} = dicomread(yya_glob{6});
            slice_data_w{count} = dicominfo(yya_glob{6});
        end
    end
end

%%
wx = 64;
wy = 64;
%for n = 1:length(Names)
    for n = 8:8
    name = Names{n};
    for t = 1:length(time_label)
        dst_dir = cat(2, base_dir, '/ContourData_Invivo/', name,'/', name, '_', time_label{t}, '/');
        img_save_dir = cat(2, base_dir, '/img/InvivoCases/');
        if ~exist(img_save_dir)
            mkdir(img_save_dir);
        end

        load(cat(2, dst_dir, 'BloodPool/mask_blood.mat'));
        load(cat(2, dst_dir, label, '_SliceLoc.mat'));

        [slc_array_sorted, idx_sorted] = sort(slc_array);
        bloodPool_tight = zeros(size(mask_blood_3D{1},1), size(mask_blood_3D{1},2), length(idx_sorted));

        for ii = 1:length(mask_blood_3D)
            bloodPool_tight(:,:,ii) = mask_blood_3D{idx_sorted(ii)};
        end

        s = regionprops(bloodPool_tight,'centroid');
        centroid = round(s.Centroid);

        figure('Position',[0, 100, 800, 600]);
        im_crop = imcrop(rescaled_t2star{n}, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        imagesc(im_crop); axis image;
        axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        colormap(brewermap([],'*OrRd'));
        caxis([0 50]);
        filename = cat(2, img_save_dir, name, '_Map_OrRd.tif');
        saveas(gcf, filename);

        figure('Position',[0, 100, 800, 600]);
        im_crop = imcrop(whatsinit_w{n}, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        imagesc(im_crop); axis image;
        axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        colormap gray;
        %colormap(brewermap([],'*OrRd'));
        %caxis([0 50]);

        filename = cat(2, img_save_dir, name, '_TE6.tif');
        saveas(gcf, filename);
    end
    end

    %% Save the color bar
    wy = 128;
    wx = 128;
    figure();
    imagesc(zeros(wy,wx));
    colormap(brewermap([],'*OrRd'));

