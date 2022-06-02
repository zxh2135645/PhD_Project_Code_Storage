close all;
clear all;
% AUC analysis across all cases

addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aha_analysis (mat file)
% from T2star_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
label = labels{5};

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name_cell = {'18P90', '18P93', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48', '20P40'};
avg_num_cell = {'Avg0016', 'Invivo'};

%% Read 0.3x0.3x2 T2* value in hemorrhage zone Avg16
% 18P90 - 70
% 18P93 - 73
% 20P03_Exvivo5 - 64
% 20P10_Exvivo7 - 79
% 20P11_Exvivo6 - 71
% 18P92 - 65
% 18P94_Exvivo3 - 94
% 18P95 - 65
% 17P73 - 79
% 20P48 - 85
% 20P40 - 69
[70, 73, 79, 71, 65, 94, 65, 79, 85, 69]
% whatsinit = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    disp(subject_name);
    base_dir = uigetdir;
    folder_glob = glob(cat(2, base_dir, '\*'));
    [list_to_read, order_to_read] = NamePicker(folder_glob);
    f = list_to_read{order_to_read(1)};
    whatsinit{i} = dicom23D(f);
end

% Load masks
mask_cell = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
    mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
end

%% Display image overlay with Mean-2SD mask
mean_t2star_array = zeros(length(subject_name_cell), 1);
sd_t2star_array = zeros(length(subject_name_cell), 1);
mean_t2star_array_remote = zeros(length(subject_name_cell), 1);
sd_t2star_array_remote = zeros(length(subject_name_cell), 1);
%t2star_distribution = cell(length(subject_name_cell), 1);
t2star_array = [];
name_array = char(subject_name_cell);
name_array_full = blanks(length(name_array));
count = 1;
figure();
ks_array = zeros(length(subject_name_cell), 1);
% se = strel('disk',1);

for i = 1:length(whatsinit)
    img2 = whatsinit{i};
    thresh5 = mean(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask)) - 5*std(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask));
    thresh2 = mean(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask)) - 2*std(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask));
    
    hemo_mask2 = (img2 < thresh2) .* mask_cell{i}.mask_struct(1).mi_mask;
    hemo_mask5 = (img2 < thresh5) .* mask_cell{i}.mask_struct(1).mi_mask;
    
    CC = bwconncomp(hemo_mask2);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx_array = find(numPixels == 1);
    for id = 1:length(idx_array)
        hemo_mask2(CC.PixelIdxList{idx_array(id)}) = 0;
    end
    
    mean_t2star_array(i) = mean(nonzeros(img2 .* hemo_mask2));
    sd_t2star_array(i) = std(nonzeros(img2 .* hemo_mask2));
    
    mean_t2star_array_remote(i) = mean(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask));
    sd_t2star_array_remote(i) = std(nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask));
    ks_array(i) = kstest((nonzeros(img2 .* mask_cell{i}.mask_struct(1).remote_mask) - mean_t2star_array_remote(i))/sd_t2star_array_remote(i));
    %t2star_distribution{i} = nonzeros(img2 .* hemo_mask2);
    t2star_array = [t2star_array; nonzeros(img2 .* hemo_mask2)];
    
    for s = 1:length(nonzeros(img2 .* hemo_mask2))
        name_array_full(count,:) = name_array(i,:);
        count = count + 1;
    end
    
    subplot(1,2,1);
    imagesc(img2); axis image; caxis([0 100]);
    title(subject_name_cell{i});
    subplot(1,2,2);
    imagesc(hemo_mask2); axis image;
    pause; 
end

%% Plot T2star values in 0.3x0.3x2 mm3
figure();
x = 1:1:length(subject_name_cell);
errorbar(x, mean_t2star_array, sd_t2star_array, 'LineWidth', 2);
xticklabels(subject_name_cell);

figure();
%boxplot(t2star_array, name_array_full, 'PlotStyle','compact');
boxplot(t2star_array, name_array_full);

%% 10/14/2021
figure();
x = 1:1:length(subject_name_cell);
errorbar(x, mean_t2star_array, sd_t2star_array, '-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2);
hold on;
errorbar(x, mean_t2star_array_remote, sd_t2star_array_remote, '-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2);
%xticklabels(subject_name_cell);
ylim([10 50]);
legend({'Hemorrhage', 'Remote'});
xlim([0, 11]); 
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
%set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[10 20 30 40 50]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');
%%
figure();
rowvals = [1 2 4 8]';
x = bsxfun(@plus,rowvals,randn(4,20));
boxplot(x','ori','horizontal','positions',rowvals,'PlotStyle','compact')
% Then maybe get rid of the group numbers:
set(gca,'ytickmode','auto')
set(gca,'yticklabelmode','auto')