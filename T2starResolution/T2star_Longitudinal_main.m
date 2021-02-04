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
auc_subjects_mean_avg16 = zeros(length(subject_name_cell), 1);
auc_subjects_mean_invivo = zeros(length(subject_name_cell), 1);

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    for j = 1:length(avg_num_cell)
        avg_name = avg_num_cell{j};
        % subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));

        if j == 1
            aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_all16 = [perc_all16, aha16.aha_analysis.perc_array_mi];
            auc_subjects_mean_avg16(i) = mean(aha16.aha_analysis.auc_array_mi);
            
        elseif j == 2
            aha01 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_all01 = [perc_all01, aha01.aha_analysis2.perc_array_mi];
%             if size(aha01.aha_analysis2.perc_array_mi,2) ~= size(aha16.aha_analysis.perc_array_mi,2)
%                 disp(subject_name)
%             end
        elseif j == 3
            aha_invivo = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_allvivo = [perc_allvivo, aha_invivo.aha_analysis2.perc_array_mi];
            auc_subjects_mean_invivo(i) = mean(aha_invivo.aha_analysis2.auc_array_mi);
        end
    end
end

gt = perc_all16(1,:) > 0.1;

%% 1.1 AUC analysis compilation
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

%% 1.2 Save as struct
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

%% 2.1 transmurality by only looking at Avg0016
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
%% 2.2 Plot and save data
figure();
s = 1:length(subject_name_cell);
[trans16_avg_sorted, I]= sort(trans16_avg);
plot(s, trans16_avg_sorted, 'LineWidth', 2)
hold on;
scatter(s, trans16_avg_sorted, 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]); ylim([0 0.3])
ylabel('Transmurality'); xlabel('Subject Name');
xticklabels(subject_name_cell(I));

figure();
s = 1:length(subject_name_cell);
[trans16_avg_sorted, I]= sort(trans16_avg);
%plot(s, trans16_avg_sorted, 'LineWidth', 2)
hold on;
bar(s, trans16_avg_sorted); ylim([0 0.3]);
grid on;
bar(2, trans16_avg_sorted(2), 'FaceColor', [0.8500, 0.3250, 0.0980]); ylim([0 0.3]);
bar(9, trans16_avg_sorted(9), 'FaceColor', [0.8500, 0.3250, 0.0980]); ylim([0 0.3]);

ylabel('Transmurality'); xlabel('Subject Name');
set(gca, 'FontSize', 24)
%xticklabels(subject_name_cell(I));

%% 2b. 
X = [ones(length(trans16_avg_sorted),1), trans16_avg_sorted];
figure();
scatter(trans16_avg_sorted, auc_subjects_mean_avg16(I), 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
hold on;
scatter(trans16_avg_sorted, auc_subjects_mean_invivo(I), 72, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
xlabel('Transmurality'); ylabel('Mean AUC');




b = X \ auc_subjects_mean_avg16(I);
yCalc_avg16 = X*b;
plot(trans16_avg_sorted, yCalc_avg16, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);

b = X \ auc_subjects_mean_invivo(I);
yCalc_invivo = X*b;
plot(trans16_avg_sorted, yCalc_invivo, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
grid on;
xlim([0.03, 0.22]);
legend({'Avg0016', 'Invivo'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 24);
Rsq1 = 1 - sum((auc_subjects_mean_avg16(I) - yCalc_avg16).^2)/sum((auc_subjects_mean_avg16 - mean(auc_subjects_mean_avg16)).^2);
Rsq2 = 1 - sum((auc_subjects_mean_invivo(I) - yCalc_invivo).^2)/sum((auc_subjects_mean_invivo - mean(auc_subjects_mean_invivo)).^2);

%% 3.1 T2* value in hemorrhage zone Avg16
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
%% 3.2 Load masks
mask_cell = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
    mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
end

%% 3.3 mean - 2sd
hemo_avg = zeros(length(subject_name_cell), 1);
remote_avg = zeros(length(subject_name_cell), 1);

for i = 1:length(subject_name_cell)
    mi = mask_cell{i}.mask_struct(1).mi_mask;
    remote = mask_cell{i}.mask_struct(1).remote_mask;
    thresh = mean(nonzeros(remote.* whatsinit{i})) - 2*std(nonzeros(remote.* whatsinit{i}));
    hemo_f = whatsinit{i} < thresh;
    hemo = hemo_f .* mi .* whatsinit{i};
    hemo_avg(i) = mean(nonzeros(hemo));
    remote_myo = remote .* whatsinit{i};
    remote_avg(i) = mean(nonzeros(remote_myo));
end
%% 3.4 Plot and save data
s = 1:length(subject_name_cell);
figure();
plot(s, hemo_avg);
hold on;
plot(s, remote_avg);

figure();

[trans16_avg_sorted, I] = sort(trans16_avg);
plot(trans16_avg_sorted, hemo_avg(I), 'LineWidth', 2)
hold on;
plot(trans16_avg_sorted, remote_avg(I), 'LineWidth', 2)
%scatter(s, trans16_avg_sorted, 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]); ylim([0 0.3])
xlabel('Transmurality'); ylabel('T2*');
%xticklabels(subject_name_cell(I));

%% 3.5 save. transmurality by only looking at Avg0016 (Shall I deprecate it?)
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

%% 4. SNR capped
snr_remote_avg16_3d = zeros(28, 5, length(subject_name_cell));
snr_remote_invivo_3d = zeros(20, 5, length(subject_name_cell));
snr_air_avg16_3d = zeros(28, 5, length(subject_name_cell));
snr_air_invivo_3d = zeros(20, 5, length(subject_name_cell));

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    avg_name16 = avg_num_cell{1};
    snr_avg16 = load(cat(2, subject_data_dir, 'SNR_', avg_name16, '.mat'));
    avg_name_invivo = avg_num_cell{3};
    snr_invivo = load(cat(2, subject_data_dir, 'SNR_', avg_name_invivo, '.mat'));
    snr_remote_avg16_3d(:,:,i) = snr_avg16.SNR.snr_remote;
    snr_remote_invivo_3d(:,:,i) = snr_invivo.SNR.snr_remote;
    snr_air_avg16_3d(:,:,i) = snr_avg16.SNR.snr_air;
    snr_air_invivo_3d(:,:,i) = snr_invivo.SNR.snr_air;
end

mean_snr_remote_avg16 = mean(snr_remote_avg16_3d, 3);
mean_snr_remote_invivo = mean(snr_remote_invivo_3d, 3);
mean_snr_air_avg16 = mean(snr_air_avg16_3d, 3);
mean_snr_air_invivo = mean(snr_air_invivo_3d, 3);
sd_snr_remote_avg16 = std(snr_remote_avg16_3d, 0, 3);
sd_snr_remote_invivo = std(snr_remote_invivo_3d, 0, 3);
sd_snr_air_avg16 = std(snr_air_avg16_3d, 0, 3);
sd_snr_air_invivo = std(snr_air_invivo_3d, 0, 3);
%% 4.1 SNR in remote
avg_array = {'Avg0016', 'Invivo'};
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
TE_avg16 = [4.95, 13.7, 22.5, 31.2, 40.0];
TE_invivo = [2.55, 5.80, 9.90, 15.6, 21.2];
for i = 1:size(mean_snr_remote_avg16, 2)
    
    subplot(2,3,i);
    avg_temp = mean_snr_remote_avg16(:,i);
    avg_sd_temp = sd_snr_remote_avg16(:,i);
    avg_x = reshape(1:length(avg_temp), [], 4).';
    avg_compare = reshape(avg_temp, [], 4).';
    avg_err = reshape(avg_sd_temp, [], 4).';
    % plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});hold on;
    plotHandles(:,1) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;

    avg_temp = mean_snr_remote_invivo(:,i);
    avg_sd_temp = sd_snr_remote_invivo(:,i);
    avg_x = reshape(mask_idx_array, [], 4).';
    avg_compare = reshape(avg_temp, [], 4).';
    avg_err = reshape(avg_sd_temp, [], 4).';
    % plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});hold on;
    plotHandles(:,2) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2});hold on;
    
    set(gca, 'FontSize', 18);
    grid on;
    set(gca,'xticklabel',{[]}); ylabel('SNR_{remote} (A.U.)');
    title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
end
    legend(plotHandles(1,:), avg_array, 'Location', 'southeast');
%% 4.2 SNR in air
avg_array = {'Avg0016', 'Invivo'};
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
TE_avg16 = [4.95, 13.7, 22.5, 31.2, 40.0];
TE_invivo = [2.55, 5.80, 9.90, 15.6, 21.2];
for i = 1:size(mean_snr_air_avg16, 2)
    
    subplot(2,3,i);
    avg_temp = mean_snr_air_avg16(:,i);
    avg_sd_temp = sd_snr_air_avg16(:,i);
    avg_x = reshape(1:length(avg_temp), [], 4).';
    avg_compare = reshape(avg_temp, [], 4).';
    avg_err = reshape(avg_sd_temp, [], 4).';
    plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});hold on;
    %plotHandles(:,1) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;
    
    avg_temp = mean_snr_air_invivo(:,i);
    avg_sd_temp = sd_snr_air_invivo(:,i);
    avg_x = reshape(mask_idx_array, [], 4).';
    avg_compare = reshape(avg_temp, [], 4).';
    avg_err = reshape(avg_sd_temp, [], 4).';
    plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});hold on;
    
    set(gca, 'FontSize', 18);
    grid on;
    set(gca,'xticklabel',{[]}); ylabel('SNR_{air} (A.U.)');
    title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
end
    legend(plotHandles(1,:), avg_array, 'Location', 'northwest');

 %% 4.3 SNR average over 5 TEs
mean_snr_remote_avg16_1d = mean(mean(snr_remote_avg16_3d, 3), 2);
mean_snr_remote_invivo_1d = mean(mean(snr_remote_invivo_3d, 3), 2);
mean_snr_air_avg16_1d = mean(mean(snr_air_avg16_3d, 3), 2);
mean_snr_air_invivo_1d = mean(mean(snr_air_invivo_3d, 3), 2);
sd_snr_remote_avg16_1d = std(reshape(snr_remote_avg16_3d, 28, []), 0, 2);
sd_snr_remote_invivo_1d = std(reshape(snr_remote_invivo_3d, 20, []), 0, 2);
sd_snr_air_avg16_1d = std(reshape(snr_air_avg16_3d, 28, []), 0, 2);
sd_snr_air_invivo_1d = std(reshape(snr_air_invivo_3d, 20, []), 0, 2);

d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_snr_remote_avg16_1d;
avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim_lb = 0; ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)

plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});
%plotHandles(:,1) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;

avg_temp = mean_snr_remote_invivo_1d;
avg_sd_temp = sd_snr_remote_invivo_1d;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});
% plotHandles(:,2) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2});hold on;

set(gca, 'FontSize', 18);
grid on;
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
ylabel('SNR_{remote} (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));

legend(plotHandles(1,:), avg_array, 'Location', 'southeast');


d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(avg_array));

avg_temp = mean_snr_air_avg16_1d;
avg_sd_temp = sd_snr_air_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)

plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});
%plotHandles(:,1) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;

avg_temp = mean_snr_air_invivo_1d;
avg_sd_temp = sd_snr_air_invivo_1d;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});
%plotHandles(:,2) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2});
set(gca, 'FontSize', 18);
grid on;
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
ylabel('SNR_{air} (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
legend(plotHandles(1,:), avg_array, 'Location', 'northwest');

%% 5.1 Fitting Residual

mean_res_avg16 = zeros(28, length(subject_name_cell));
mean_res_invivo = zeros(20, length(subject_name_cell));

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    avg_name16 = avg_num_cell{1};
    FitResults_avg16 = load(cat(2, subject_data_dir, 'FitResults_', avg_name16, '.mat'));
    avg_name_invivo = avg_num_cell{3};
    FitResults_invivo = load(cat(2, subject_data_dir, 'FitResults_', avg_name_invivo, '.mat'));
    
    for j = 1:size(mean_res_avg16, 1)
        mean_res_avg16(j, i) = mean(mean(FitResults_avg16.FitResults_struct(j).FitResults.res(~isnan(FitResults_avg16.FitResults_struct(j).FitResults.res))));
    end
    for j = 1:length(mean_res_invivo)
        mean_res_invivo(j, i) = mean(mean(FitResults_invivo.FitResults_struct(j).FitResults.res(~isnan(FitResults_invivo.FitResults_struct(j).FitResults.res))));        
    end
end

mean_mean_res_avg16 = mean(mean_res_avg16, 2);
mean_mean_res_invivo = mean(mean_res_invivo, 2);

%% 5.2 Plot fitting residual
d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
plot(avg_x.', avg_compare.', 'LineStyle', 'none' );

hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)


avg_temp = mean_mean_res_avg16;
%avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
%avg_err = reshape(avg_sd_temp, [], 4).';

plotHandles(:,1) = plot(avg_x.', avg_compare.', '-o',  'LineWidth', 2, 'Color', color_cell{1}); hold on;
%plot(mean_mean_res_invivo);

avg_temp = mean_mean_res_invivo;
%avg_sd_temp = sd_snr_remote_invivo_1d;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
%avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,2) = plot(avg_x.', avg_compare.', '-o', 'LineWidth', 2, 'Color', color_cell{2}); hold on;

set(gca, 'FontSize', 18);
grid on;
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub]);

text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
ylabel('Fitting Residual (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
legend(plotHandles(1,:), avg_array, 'Location', 'northeast');