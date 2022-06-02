clear all; close all;

% A new T2* value analysis with
% All resolution in both Ideal and Practical
% At both hemorrahge and remote myocardium

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

%% Read T2* value in hemorrhage zone Avg16
sel_array_avg16_auto = [[70, 68, 66, 64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16];...
    [73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19];...
    [79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25];...
    [71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17];...
    [65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11];...
    [94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 39];...
    [65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11];...
    [79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25];...
    [85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31];...
    [69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15]];
sel_array_invivo_auto = [[167, 165, 163, 161, 159, 157, 155, 153, 151, 149, 147, 145, 143, 141, 139, 137, 135, 133, 131, 129];...
    [170, 168, 166, 164, 162, 160, 158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 138, 136, 134, 132];...
    [170, 172, 174, 176, 178, 160, 162, 164, 166, 168, 150, 152, 154, 156, 158, 140, 142, 144, 146, 148];...
    [162, 164, 166, 168, 170, 152, 154, 156, 158, 160, 142, 144, 146, 148, 150, 132, 134, 136, 138, 140];...
    [162, 160, 158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 138, 136, 134, 132, 130, 128, 126, 124];...
    [191, 189, 187, 185, 183, 181, 179, 177, 175, 173, 171, 169, 167, 165, 163, 161, 159, 157, 155, 153];...
    [162, 160, 158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 138, 136, 134, 132, 130, 128, 126, 124];...
    [176, 174, 172, 170, 168, 166, 164, 162, 160, 158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 138];...
    [118, 120, 122, 124, 126, 108, 110, 112, 114, 116, 97, 100, 102, 104, 106, 87, 89, 91, 93, 95];...
    [110, 108, 106, 104, 102, 100, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71]];
base_dir_auto = {'/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_18P90_EXVIVO2_FORMALIN/BIRI_RESEARCH_JAMES_20200905_221727_180000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_18P93_EXVIVO2_FORMALIN/BIRI_RESEARCH_JAMES_20200913_025501_944000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_20P10_EXVIVO7_FORMALIN/BIRI_RESEARCH_JAMES_20200829_000319_085000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_20P11_EXVIVO6_FORMALIN/BIRI_RESEARCH_JAMES_20200829_233553_298000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_18P92_EXVIVO2_FORMALIN/BIRI_RESEARCH_JAMES_20200920_002138_652000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_18P94_EXVIVO3_FORMALIN/BIRI_RESEARCH_JAMES_20201002_215427_259000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_18P95_EXVIVO2_FORMALIN/BIRI_RESEARCH_JAMES_20200926_231722_771000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/DHARMAKUMAR_17P73_EXVIVO2_FORMALIN/BIRI_RESEARCH_JAMES_20201003_215141_116000',...
    '/Users/jameszhang/Documents/Cedars_pig_2020/JAMES_PHANTOM_20P48_EXVIVO/BIRI_RESEARCH_JAMES_20210220_195324_385000',...
    '/Users/jameszhang/Documents/Cedars_pig_2021/20P40_EXVIVO2_JAMES/BIRI_RESEARCH_JAMES_20210508_171939_321000'};

img_struct = struct;
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    disp(subject_name);
    % base_dir = uigetdir;
    base_dir = base_dir_auto{i};
    folder_glob = glob(cat(2, base_dir, '\*'));
    [list_to_read, order_to_read] = NamePicker(folder_glob, 0, sel_array_avg16_auto(i,:));
    
    for j = 1:length(order_to_read)
        f = list_to_read{order_to_read(j)};
        whatsinit_avg16{j} = dicom23D(f);
    end
    
   [list_to_read, order_to_read] = NamePicker(folder_glob, 0, sel_array_invivo_auto(i,:));

    for j = 1:length(order_to_read)
        f = list_to_read{order_to_read(j)};
        whatsinit_invivo{j} = dicom23D(f);
    end
    
    img_struct(i).whatsinit_avg16 = whatsinit_avg16;
    img_struct(i).whatsinit_invivo = whatsinit_invivo;
end

%% Load masks
mask_cell = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
    mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
end

%% Display image - remote
mean_t2star_array = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
sd_t2star_array = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
mean_t2star_array_remote_avg16 = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
sd_t2star_array_remote_avg16 = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
mean_t2star_array_remote_invivo = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_invivo));
sd_t2star_array_remote_invivo = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_invivo));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];

for i = 1:length(subject_name_cell)
    for j = 1:length(img_struct(1).whatsinit_avg16)
        img_avg16 = img_struct(i).whatsinit_avg16{j};
        remote_mask = mask_cell{i}.mask_struct(j).remote_mask;
        mean_t2star_array_remote_avg16(i,j) = mean(nonzeros(img_avg16 .* remote_mask));
        sd_t2star_array_remote_avg16(i,j) = std(nonzeros(img_avg16 .* remote_mask));
    end
    
    for j = 1:length(img_struct(i).whatsinit_invivo)
        img_invivo = img_struct(i).whatsinit_invivo{j};
        remote_mask = mask_cell{i}.mask_struct(mask_idx_array(j)).remote_mask;
        mean_t2star_array_remote_invivo(i,j) = mean(nonzeros(img_invivo .* remote_mask));
        sd_t2star_array_remote_invivo(i,j) = std(nonzeros(img_invivo .* remote_mask));
    end
end

figure();
img_e = img_struct(1).whatsinit_invivo{1};
% img_e = img_struct(1).whatsinit_avg16{1};
remote_mask = mask_cell{1}.mask_struct(mask_idx_array(1)).remote_mask;
subplot(1,2,1);
imagesc(img_e); axis image; caxis([0 100]);
subplot(1,2,2);
imagesc(img_e.*remote_mask); axis image; caxis([0 100]);
%% Plot (remote)
mean_mean_t2star_array_remote_avg16 = mean(mean_t2star_array_remote_avg16, 1);
sd_mean_t2star_array_remote_avg16 = std(mean_t2star_array_remote_avg16, 0, 1);
mean_mean_t2star_array_remote_invivo = mean(mean_t2star_array_remote_invivo, 1);
sd_mean_t2star_array_remote_invivo = std(mean_t2star_array_remote_invivo, 0, 1);

avg_array = {'Avg0016', 'Invivo'};
d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
hax = axes;
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_mean_t2star_array_remote_avg16;
avg_sd_temp = sd_mean_t2star_array_remote_avg16;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
hE = errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim([0 60]);
ylim_lb = 0; ylim_ub = max(ylim);

patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)
plotHandles(:,3) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;
plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});


avg_temp = mean_mean_t2star_array_remote_invivo;
avg_sd_temp = sd_mean_t2star_array_remote_invivo;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,4) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2}); hold on;
plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});

set(gca, 'FontSize', 18);
%grid on;
%xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
%xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4})
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2})

color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4})
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2})

set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');

%% Statistical Analysis
str_cell = {'030302', '060602', '080802', '101002', '131302', '161602', '212102', '030304', '060604', '080804', '101004', '131304', '161604', '212104','030306', '060606', '080806', '101006', '131306', '161606', '212106','030308', '060608', '080808', '101008', '131308', '161608', '212108'};
str_cell2 = {'Avg0016', 'Invivo'};
N = length(subject_name_cell);

p1 = anova1(mean_t2star_array_remote_avg16); % P1 = 1.0000
p2 = kruskalwallis(mean_t2star_array_remote_avg16); % P2 = 1.0000

p1_invivo = anova1(mean_t2star_array_remote_invivo); % P1_invivo = 0.8266
p2_invivo = kruskalwallis(mean_t2star_array_remote_invivo); % P2_invivo = 0.9508

p3 = mean_t2star_array_remote_avg16(:);
p4 = mean_t2star_array_remote_invivo(:);

strength = [];
alloy = {};

strength = [p3', p4'];
alloy = [alloy, repmat({str_cell2{1}}, [1,length(p3)]), repmat({str_cell2{2}}, [1,length(p4)])];

p5 = anova1(strength, alloy); % 0.0589
p6 = kruskalwallis(strength, alloy); % 0.1047

%% Display image - IMH
mean_t2star_array = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
sd_t2star_array = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
mean_t2star_array_mi_avg16 = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
sd_t2star_array_mi_avg16 = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_avg16));
mean_t2star_array_mi_invivo = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_invivo));
sd_t2star_array_mi_invivo = zeros(length(subject_name_cell), length(img_struct(1).whatsinit_invivo));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];

for i = 1:length(subject_name_cell)
    for j = 1:length(img_struct(1).whatsinit_avg16)
        img_avg16 = img_struct(i).whatsinit_avg16{j};
        remote_mask = mask_cell{i}.mask_struct(j).remote_mask;
        mi_mask = mask_cell{i}.mask_struct(j).mi_mask;
        myo_mask = mask_cell{i}.mask_struct(j).myo_mask;
        thresh = mean(nonzeros(img_avg16 .* remote_mask)) - 5 * std(nonzeros(img_avg16 .* remote_mask));
        imh_mask = (img_avg16 <= thresh) .* mi_mask .* myo_mask;
        
        CC = bwconncomp(imh_mask);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        idx_array = find(numPixels <= 1);
        for id = 1:length(idx_array)
            imh_mask(CC.PixelIdxList{idx_array(id)}) = 0;
        end
    
        mask_cell{i}.mask_struct(j).imh_mask = imh_mask;
        mean_t2star_array_mi_avg16(i,j) = mean(nonzeros(img_avg16 .* imh_mask));
        sd_t2star_array_mi_avg16(i,j) = std(nonzeros(img_avg16 .* imh_mask));
    end
    
    for j = 1:length(img_struct(i).whatsinit_invivo)
        img_invivo = img_struct(i).whatsinit_invivo{j};
        remote_mask = mask_cell{i}.mask_struct(mask_idx_array(j)).remote_mask;
        mi_mask = mask_cell{i}.mask_struct(mask_idx_array(j)).mi_mask;
        myo_mask = mask_cell{i}.mask_struct(mask_idx_array(j)).myo_mask;
        
        thresh = mean(nonzeros(img_invivo .* remote_mask)) - 2 * std(nonzeros(img_invivo .* remote_mask));
        imh_mask = (img_invivo <= thresh) .* mi_mask .* myo_mask;
        mean_t2star_array_mi_invivo(i,j) = mean(nonzeros(img_invivo .* imh_mask));
        sd_t2star_array_mi_invivo(i,j) = std(nonzeros(img_invivo .* imh_mask));
    end
end

figure();
% img_e = img_struct(1).whatsinit_invivo{1};
img_e = img_struct(1).whatsinit_avg16{1};
remote_mask = mask_cell{1}.mask_struct(1).remote_mask;
subplot(1,2,1);
imagesc(img_e); axis image; caxis([0 100]);
subplot(1,2,2);
imagesc(img_e.*remote_mask); axis image; caxis([0 100]);