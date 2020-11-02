% The main function for Longitudinal analysis
% Some scripts needed to be run before excuting this code
clear all;
close all;
addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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

subject_name_cell = {'18P90', '18P93', '20P03_Exvivo5', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73'};
avg_num_cell = {'Avg0016', 'Avg0001', 'Invivo'};


%perc_all16 = zeros(28, length(subject_name_cell));
%perc_all01 = zeros(28, length(subject_name_cell));
%perc_allvivo = zeros(20, length(subject_name_cell));
perc_all16 = [];
perc_all01 = [];
perc_allvivo = [];
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    for j = 1:length(avg_num_cell)
        avg_name = avg_num_cell{j};
        % subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));

        if j == 1
            aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_all16 = [perc_all16, aha16.aha_analysis.perc_array_mi];
        elseif j == 2
            aha01 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_all01 = [perc_all01, aha01.aha_analysis2.perc_array_mi];
%             if size(aha01.aha_analysis2.perc_array_mi,2) ~= size(aha16.aha_analysis.perc_array_mi,2)
%                 disp(subject_name)
%             end
        elseif j == 3
            aha_invivo = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_allvivo = [perc_allvivo, aha_invivo.aha_analysis2.perc_array_mi];
        end
    end
end

gt = perc_all16(1,:) > 0.1;

%% 1. AUC analysis compilation
% Confusion Matrix
% figure('Position', [100 0 1600 1600]);
% mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
% 
% for i = 1:size(perc_all16, 1)
%     [cm, order] = confusionmat(gt,perc_all16(i,:)>0.1);
%     subplot(4,7,i);confusionchart(cm, order); 
%     set(gca, 'FontSize', 18);
% end

% ROC analysis %% Avg0016
auc_all16 = zeros(size(perc_all16, 1),1);
figure('Position', [100 0 1600 1600]);
for i = 1:size(perc_all16, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(gt,perc_all16(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    set(gca, 'FontSize', 16);
    auc_all16(i) = AUC;
end

% ROC analysis %% Avg0001
auc_all01 = zeros(size(perc_all01, 1),1);
figure('Position', [100 0 1600 1600]);
for i = 1:size(perc_all01, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(gt,perc_all01(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    set(gca, 'FontSize', 16);
    auc_all01(i) = AUC;
end

% ROC analysis %% Avg0001
auc_allvivo = zeros(size(perc_allvivo, 1),1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
figure('Position', [100 0 1600 1600]);
for i = 1:size(perc_all01, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(gt,perc_allvivo(i,:), 1);
    subplot(4,5,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    set(gca, 'FontSize', 16);
    auc_allvivo(i) = AUC;
end

%% Save as struct
auc_compilation = struct;
auc_compilation.auc_all16 = auc_all16;
auc_compilation.auc_all01 = auc_all01;
auc_compilation.auc_allvivo = auc_allvivo;
auc_compilation.subject_name_cell = subject_name_cell;

subject_data_dir = GetFullPath(cat(2, data_dir, 'Compilation/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end
save(cat(2, subject_data_dir, 'auc_compilation.mat'), 'auc_compilation');

%% 2. transmurality by only looking at Avg0016
perc_trans16 = cell(length(subject_name_cell), 1);
trans16_avg = zeros(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
        avg_name = avg_num_cell{1};
        aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
        perc_array_mi = aha16.aha_analysis.perc_array_mi;
        perc_array_temp = nonzeros(perc_array_mi(1,:))';
        perc_trans16{i} = perc_array_temp;
        trans16_avg(i) = mean(perc_array_temp);
end
%% Plot and save data
figure();
s = 1:length(subject_name_cell);
[trans16_avg_sorted, I]= sort(trans16_avg);
plot(s, trans16_avg_sorted, 'LineWidth', 2)
hold on;
scatter(s, trans16_avg_sorted, 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]); ylim([0 0.3])
ylabel('Transmurality'); xlabel('Subject Name');
xticklabels(subject_name_cell(I));
%% 3. T2* value in hemorrhage zone
% [70, 73, 64, 79, 71, 65, 94, 65, 79]
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
%% Load masks
mask_cell = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
    mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
end

%% mean - 2sd
hemo_avg = zeros(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    mi = mask_cell{i}.mask_struct(1).mi_mask;
    remote = mask_cell{i}.mask_struct(1).remote_mask;
    thresh = mean(nonzeros(remote.* whatsinit{i})) - 2*std(nonzeros(remote.* whatsinit{i}));
    hemo_f = whatsinit{i} < thresh;
    hemo = hemo_f .* mi .* whatsinit{i};
    hemo_avg(i) = mean(nonzeros(hemo));
end
%% Plot and save data
s = 1:length(subject_name_cell);
figure();
plot(s, hemo_avg)
figure();

[trans16_avg_sorted, I]= sort(trans16_avg);
plot(trans16_avg_sorted, hemo_avg(I), 'LineWidth', 2)
hold on;
%scatter(s, trans16_avg_sorted, 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]); ylim([0 0.3])
xlabel('Transmurality'); ylabel('T2*');
%xticklabels(subject_name_cell(I));

%% 2. transmurality by only looking at Avg0016
perc_trans16 = cell(length(subject_name_cell), 1);
trans16_avg = zeros(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
        avg_name = avg_num_cell{1};
        aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
        perc_array_mi = aha16.aha_analysis.perc_array_mi;
        perc_array_temp = nonzeros(perc_array_mi(1,:))';
        perc_trans16{i} = perc_array_temp;
        trans16_avg(i) = mean(perc_array_temp);
        
        
end