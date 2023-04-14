clear all;
close all;

%% 1. representitative subjects (20P10_Exvivo7, 18P93 -> 18P95)
addpath('../function/')
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

% disp('Avg 0016 starts here: ');
% avg_num = input('Please type average number here:  ');
% avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));

%% Load Avg16 Data
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Invivo'};
avg_name = avg_num_cell{1};
[70, 73, 69, 79, 71, 65, 94, 65, 79, 85]
whatsinit = cell(length(subject_name_cell), 1);
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
%% 3.1 Barplot of Transmurality
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Invivo'};
avg_name = avg_num_cell{1};

perc_trans16 = cell(length(subject_name_cell), 1);
trans16_avg = zeros(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    % subject_name = input('Please type subject name here:  ', 's');
    subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
    if ~exist(subject_dir, 'dir')
        mkdir(subject_dir)
    end
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    if ~exist(subject_data_dir, 'dir')
        mkdir(subject_data_dir)
    end

    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    
    myo_mask = mask_cell{i}.mask_struct(1).myo_mask;
    remote_mask = mask_cell{i}.mask_struct(1).remote_mask;
    mi_mask = mask_cell{i}.mask_struct(1).mi_mask;
    
    thresh = mean(nonzeros(whatsinit{i} .* myo_mask .* remote_mask)) - 2*std(nonzeros(whatsinit{i} .* myo_mask .* remote_mask));
    hemo_mask = (whatsinit{i} < thresh) .* myo_mask .* mi_mask;
    trans16_avg(i) = sum(hemo_mask(:)) ./ sum(mi_mask(:));
    
    %figure();
    %imagesc(myo_mask + hemo_mask*2 + mi_mask);
%     perc_array_mi = aha16.aha_analysis.perc_array_mi;
%     perc_array_temp = nonzeros(perc_array_mi(1,:))';
%     perc_trans16{i} = perc_array_temp;
%     trans16_avg(i) = mean(perc_array_temp);
end
%%
ax = figure('Position', [100 0 600 500]);
s = 1:length(subject_name_cell);
[trans16_avg_sorted, I] = sort(trans16_avg);
%plot(s, trans16_avg_sorted, 'LineWidth', 2)
hold on;
bar(s, trans16_avg_sorted, 'FaceColor', [35, 47, 233]/255, 'EdgeColor', [35, 47, 233]/255);
%grid on;
bar(4, trans16_avg_sorted(4), 'FaceColor', [222,114,98]/255, 'EdgeColor', [222,114,98]/255); 
bar(10, trans16_avg_sorted(10), 'FaceColor', [222,114,98]/255, 'EdgeColor', [222,114,98]/255); 
ylim([0 0.25]); xlim([0.2 10.8]);
%ylabel('Transmurality'); xlabel('Subject Name');
%set(gca, 'FontSize', 24);

%set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', [4, 10]);
%set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 2,'TickLength',[0.04 0.04]);
set(gca,'TickDir','out'); % The only other option is 'in'
