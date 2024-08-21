% To pull up images for figures
% Mainly based on T2star_Longitudinal_main.m
% But made for Adobe illustrator (Publication use)
clear all;
close all;

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% aha_anlysis
% T2star_meanSD_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. representitative subjects (20P10_Exvivo7, 18P93 -> 18P95)
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

subject_name = input('Please type subject name here:  ', 's');
subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
if ~exist(subject_dir, 'dir')
    mkdir(subject_dir)
end
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

disp('Avg 0016 starts here: ');
avg_num = input('Please type average number here:  ');
avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));

%% Load data (Longitudinal)
load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
res_mi2 = aha_analysis.perc_array_mi > 0.1;
perc_array_mi = aha_analysis.perc_array_mi;

% ROC analysis 
% 20P10_Exvivo7 & 18P93
figure('Position', [100 0 1600 1600]);
k = ([229, 240, 248] - [0, 113, 188]) / (1 - 0.5); % Light Blue % Dark blue
[243, 209, 215];
[246, 101, 72];
for i = 1:size(res_mi2, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi2(1,:),perc_array_mi(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    %xlabel('FPR');
    %ylabel('TPR');
    text(0.6,0.2,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    %set(h, 'Color', 'k')
    % get rid of the white ticks and tick labels, moving the labels closer to
    % the axes
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    pos = get(gca, 'Position');
    pos(1) = 0.1 + mod(i-1, 7) * 0.09;
    pos(2) = 0.8 - fix((i-1)/7) * 0.165;
    set(gca, 'Position', pos);
    if AUC < 0.5
        set(gca,'color', [229, 240, 248]/255)
    else
        set(gca,'color', ([229, 240, 248]-k*(AUC-0.5))/255)
    end
end

%% ROC Analysis (Longitudinal)
% subject_name_cell = {'18P90', '18P93', '20P03_Exvivo5', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Invivo'};

perc_all16 = [];
perc_all01 = [];
perc_allvivo = [];
auc_subjects_mean_avg16 = zeros(length(subject_name_cell), 1);
auc_subjects_mean_invivo = zeros(length(subject_name_cell), 1);
perc_subjects_avg16_cell = cell(length(subject_name_cell), 1);
perc_subjects_invivo_cell = cell(length(subject_name_cell), 1);

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    for j = 1:length(avg_num_cell)
        avg_name = avg_num_cell{j};

        if j == 1
            aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_all16 = [perc_all16, aha16.aha_analysis.perc_array_mi];
            auc_subjects_mean_avg16(i) = mean(aha16.aha_analysis.auc_array_mi);
            perc_subjects_avg16_cell{i} = aha16.aha_analysis.perc_array_mi;
        elseif j == 2
            aha_invivo = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            perc_allvivo = [perc_allvivo, aha_invivo.aha_analysis2.perc_array_mi];
            auc_subjects_mean_invivo(i) = mean(aha_invivo.aha_analysis2.auc_array_mi);
            perc_subjects_invivo_cell{i} = aha_invivo.aha_analysis2.perc_array_mi;
        end
    end
end

gt = perc_all16(1,:) > 0.1;
%% 2.1. Avg16 (Longitudinal)
lb = 0.5;
figure('Position', [100 0 1600 1600]);
k = ([229, 240, 248] - [0, 113, 188]) / (1 - lb); % Light Blue % Dark blue

for i = 1:size(perc_all16, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(gt,perc_all16(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    %xlabel('FPR');
    %ylabel('TPR');
    text(0.6,0.2,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    %set(h, 'Color', 'k')
    % get rid of the white ticks and tick labels, moving the labels closer to
    % the axes
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    pos = get(gca, 'Position');
    pos(1) = 0.1 + mod(i-1, 7) * 0.09;
    pos(2) = 0.8 - fix((i-1)/7) * 0.165;
    set(gca, 'Position', pos);
    if AUC < lb
        set(gca,'color', [229, 240, 248]/255)
    else
        set(gca,'color', ([229, 240, 248]-k*(AUC-lb))/255)
    end
end

%% 2.2. Invivo (Longitudinal)
lb = 0.5;
figure('Position', [100 0 1000 1600]);
k = ([229, 240, 248] - [0, 113, 188]) / (1 - lb); % Light Blue % Dark blue

for i = 1:size(perc_allvivo, 1)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(gt,perc_allvivo(i,:), 1);
    subplot(4,5,i);
    plot(X,Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    %xlabel('FPR');
    %ylabel('TPR');
    text(0.6,0.2,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    %set(h, 'Color', 'k')
    % get rid of the white ticks and tick labels, moving the labels closer to
    % the axes
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    pos = get(gca, 'Position');
    pos(1) = 0.1 + mod(i-1, 5) * 0.127;
    pos(2) = 0.8 - fix((i-1)/5) * 0.165;
    set(gca, 'Position', pos);
    if AUC < lb
        set(gca,'color', [229, 240, 248]/255)
    else
        set(gca,'color', ([229, 240, 248]-k*(AUC-lb))/255)
    end
    axis tight;
end

%% 2.2.2 Invivo (Longitudinal, averaged AUC is across 10 subjects)
lb = 0.5;
figure('Position', [100 0 1000 1600]);
k = ([229, 240, 248] - [0, 113, 188]) / (1 - lb); % Light Blue % Dark blue
AUC_array = zeros(size(perc_allvivo, 1), size(perc_subjects_invivo_cell, 1));
for j = 1:size(perc_subjects_invivo_cell, 1)
    perc_invivo = perc_subjects_invivo_cell{j};
    gt = perc_subjects_avg16_cell{j} > 0.1;
    for i = 1:size(perc_allvivo, 1)
        [X,Y,T,AUC,OPTROCPT] = perfcurve(gt(1,:),perc_invivo(i,:), 1);
        AUC_array(i,j) = AUC;
    end

    for i = 1:size(perc_allvivo, 1)
        AUC_avg = mean(AUC_array(i,:));
        subplot(4,5,i);
        %plot(X, Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
        %xlabel('FPR');
        %ylabel('TPR');
        text(0.6,0.2,num2str(round(AUC_avg, 2)),'FontSize',24);
        % title('ROC for Classification')
        set(gca, 'FontSize', 16);
        set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
        %set(h, 'Color', 'k')
        % get rid of the white ticks and tick labels, moving the labels closer to
        % the axes
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        pos = get(gca, 'Position');
        pos(1) = 0.1 + mod(i-1, 5) * 0.127;
        pos(2) = 0.8 - fix((i-1)/5) * 0.165;
        set(gca, 'Position', pos);
        if AUC_avg < lb
            set(gca,'color', [229, 240, 248]/255)
        else
            set(gca,'color', ([229, 240, 248]-k*(AUC_avg-lb))/255)
        end
        axis tight;
    end
end

%% Two way ANOVA

%% 3.1 Barplot of Transmurality
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

%% 3.2 Linear regression
%X = [ones(length(trans16_avg_sorted),1), trans16_avg_sorted];

figure('Position', [100 0 600 500]);
scatter(trans16_avg_sorted, auc_subjects_mean_avg16(I), 200, 'filled', 'MarkerFaceColor', [.75 .75 1], 'MarkerEdgeColor', 'none' );
hold on;
scatter(trans16_avg_sorted, auc_subjects_mean_invivo(I), 200, 'filled', 'MarkerFaceColor', [253,190,133]/255);
%xlabel('Transmurality'); ylabel('Mean AUC');


% 
X1 = [trans16_avg_sorted];
mdl1 = fitlm(X1, auc_subjects_mean_avg16(I));
mdl1.Coefficients;
ci = coefCI(mdl1);
b1 = mdl1.Coefficients.Estimate;

fitted_avg16 = b1(1) + X1*b1(2);
fitted_avg16_lb = ci(1,1) + ci(2,1) * X1;
fitted_avg16_ub = ci(1,2) + ci(2,2) * X1;

plot(trans16_avg_sorted, fitted_avg16, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
% b = X \ auc_subjects_mean_avg16(I);
% yCalc_avg16 = X*b;
% plot(trans16_avg_sorted, yCalc_avg16, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
%line(trans16_avg_sorted,fitted_avg16_lb, 'LineStyle', '-.', 'Color', [0 .5 0])
%line(trans16_avg_sorted,fitted_avg16_ub, 'LineStyle', '-.', 'Color', [0 .5 0])


mdl2 = fitlm(X1, auc_subjects_mean_invivo(I));
mdl2.Coefficients;
ci = coefCI(mdl2);
b2 = mdl2.Coefficients.Estimate;
fitted_invivo = b2(1) + X1*b2(2);
plot(trans16_avg_sorted, fitted_invivo, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%b = X \ auc_subjects_mean_invivo(I);
%yCalc_invivo = X*b;

%grid on;
xlim([0.03, 0.22]); ylim([0.4 1]);
%legend({'Avg0016', 'Invivo'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 24);
set(gca, 'FontName', 'Helvetica')
%Rsq1 = 1 - sum((auc_subjects_mean_avg16(I) - yCalc_avg16).^2)/sum((auc_subjects_mean_avg16 - mean(auc_subjects_mean_avg16)).^2);
%Rsq2 = 1 - sum((auc_subjects_mean_invivo(I) - yCalc_invivo).^2)/sum((auc_subjects_mean_invivo - mean(auc_subjects_mean_invivo)).^2);

set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
%set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.04 0.04]);
set(gca,'TickDir','out', 'YGrid', 'on');
%% Try an example
% addpath('../function/Publication_Quality_Graphics/');
% load data xfit yfit xdata_m ydata_m ydata_s xVdata yVdata xmodel ymodel ...
%     ymodelL ymodelU c cint
% %% Create basic plot
% figure
% hold on
% hData = line(xVdata, yVdata);
% hModel = line(xmodel, ymodel);
% hCI(1) = line(xmodel, ymodelL);
% hCI(2) = line(xmodel, ymodelU);
% 
% % Adjust line properties (functional)
% set(hData, 'LineStyle', 'none', 'Marker', '.')
% set(hModel, 'LineStyle', '-', 'Color', 'r')
% set(hCI(1), 'LineStyle', '-.', 'Color', [0 .5 0])
% set(hCI(2), 'LineStyle', '-.', 'Color', [0 .5 0])
% 
% set(hData, 'Marker', 'o', 'MarkerSize', 5, ...
%     'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
% set(hModel, 'LineWidth', 1.5)
% set(hCI(1), 'LineWidth', 1.5)
% set(hCI(2), 'LineWidth', 1.5)
% 
% 
% % Add labels
% hTitle = title('My Publication-Quality Graphics');
% hXLabel = xlabel('Length (m)');
% hYLabel = ylabel('Mass (kg)');
% 
% % Add text
% hText = text(10, 800, ...
%     sprintf('{\\itC = %0.1g \\pm %0.1g (CI)}', c, cint(2)-c));
% 
% % Add legend
% hLegend = legend([hData, hModel, hCI(1)], ...
%     'Validation Data', 'Model (C{\itx}^3)', '95% CI', ...
%     'Location', 'NorthWest');
% 
% % Adjust font
% set(gca, 'FontName', 'Helvetica')
% set([hTitle, hXLabel, hYLabel, hText], 'FontName', 'AvantGarde');
% 
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:500:2500, ...
%     'LineWidth', 1)
% set([hLegend, gca], 'FontSize', 8)
% set([hXLabel, hYLabel, hText], 'FontSize', 10)
% set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

%% Read 0.3x0.3x2 T2* value in hemorrhage zone Avg16
% [70, 73, 64, 79, 71, 65, 94, 65, 79, 85]
% whatsinit = cell(length(subject_name_cell), 1);
% for i = 1:length(subject_name_cell)
%     subject_name = subject_name_cell{i};
%     disp(subject_name);
%     base_dir = uigetdir;
%     folder_glob = glob(cat(2, base_dir, '\*'));
%     [list_to_read, order_to_read] = NamePicker(folder_glob);
%     f = list_to_read{order_to_read(1)};
%     whatsinit{i} = dicom23D(f);
% end
% 
% % Load masks
% mask_cell = cell(length(subject_name_cell), 1);
% for i = 1:length(subject_name_cell)
%     subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
%     mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
% end

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
    avg_name_invivo = avg_num_cell{2};
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
figure('Position', [100 0 800 400]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_snr_remote_avg16_1d;
avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
hE = errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim([0 80]);
ylim_lb = 0; ylim_ub = max(ylim);
% patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
% patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
% patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
% patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)

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
%grid on;
%xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
%xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
ylabel('SNR_{remote} (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', [.3 .3 .3])
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7])

set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', [.4 .4 .4])
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceColor' , [.8 .8 .8])
%legend(plotHandles(1,:), avg_array, 'Location', 'southeast');
%%
% part 2
d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
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
%grid on;
%xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
%xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

% text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
% text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
% text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
% text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabelAvg16,{[]}); 
ylabel('SNR_{air} (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
%legend(plotHandles(1,:), avg_array, 'Location', 'northwest');

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', [.3 .3 .3])
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , [.7 .7 .7])

set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', [.4 .4 .4])
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [.3 .3 .3], 'MarkerFaceColor' , [.8 .8 .8])

%% 4.3 SNR average over 5 TEs (Different color scheme)
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
figure('Position', [100 0 800 400]);
hax = axes;
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_snr_remote_avg16_1d;
avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
hE = errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim([0 80]);
ylim_lb = 0; ylim_ub = max(ylim);
% patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
% patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
% patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
% patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)

patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)
plotHandles(:,3) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1});hold on;
plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});

avg_temp = mean_snr_remote_invivo_1d;
avg_sd_temp = sd_snr_remote_invivo_1d;
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

% text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
% text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
% text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
% text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
%ylabel('SNR_{remote} (A.U.)');

color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4})
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2})

color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4})
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2})
%legend(plotHandles(1,:), avg_array, 'Location', 'southeast');

set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');


%% part 2 Air
d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
plotHandles = zeros(4,length(avg_array));

avg_temp = mean_snr_air_avg16_1d;
avg_sd_temp = sd_snr_air_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
errorbar(avg_x.', avg_compare.', avg_err.', 'LineStyle', 'none' );

hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
% patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
% patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
% patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
% patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)

patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)
plotHandles(:,3) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1}); hold on;
plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{1});

avg_temp = mean_snr_air_invivo_1d;
avg_sd_temp = sd_snr_air_invivo_1d;
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

% text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
% text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
% text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
% text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabelAvg16,{[]}); 
%ylabel('SNR_{air} (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
%legend(plotHandles(1,:), avg_array, 'Location', 'northwest');

color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');

%% 5.1 Fitting Residual and R2
avg_array = {'Avg0016', 'Invivo'};

mean_res_avg16 = zeros(28, length(subject_name_cell));
mean_res_invivo = zeros(20, length(subject_name_cell));

mean_r2_avg16 = zeros(28, length(subject_name_cell));
mean_r2_invivo = zeros(20, length(subject_name_cell));

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    avg_name16 = avg_num_cell{1};
    FitResults_avg16 = load(cat(2, subject_data_dir, 'FitResults_', avg_name16, '.mat'));
    avg_name_invivo = avg_num_cell{2};
    FitResults_invivo = load(cat(2, subject_data_dir, 'FitResults_', avg_name_invivo, '.mat'));
    
    for j = 1:size(mean_res_avg16, 1)
        mean_res_avg16(j, i) = mean(mean(FitResults_avg16.FitResults_struct(j).FitResults.res(~isnan(FitResults_avg16.FitResults_struct(j).FitResults.res))));
    end
    for j = 1:length(mean_res_invivo)
        mean_res_invivo(j, i) = mean(mean(FitResults_invivo.FitResults_struct(j).FitResults.res(~isnan(FitResults_invivo.FitResults_struct(j).FitResults.res))));        
    end
    
    for j = 1:size(mean_r2_avg16, 1)
        temp = FitResults_avg16.FitResults_struct(j).FitResults.R2(~isnan(FitResults_avg16.FitResults_struct(j).FitResults.R2));
        temp(temp<=0) = 1;
        mean_r2_avg16(j, i) = mean(temp);
    end
    for j = 1:length(mean_r2_invivo)
        temp_invivo = FitResults_invivo.FitResults_struct(j).FitResults.R2(~isnan(FitResults_invivo.FitResults_struct(j).FitResults.R2));
        temp_invivo(temp_invivo<=0) = 1;
        mean_r2_invivo(j, i) = mean(temp_invivo);
    end
end

mean_mean_res_avg16 = mean(mean_res_avg16, 2);
mean_mean_res_invivo = mean(mean_res_invivo, 2);
std_res_avg16 = std(mean_res_avg16, 0, 2);
std_res_invivo = std(mean_res_invivo, 0, 2);

mean_mean_r2_avg16 = mean(mean_r2_avg16, 2);
mean_mean_r2_invivo = mean(mean_r2_invivo, 2);
std_r2_avg16 = std(mean_r2_avg16, 0, 2);
std_r2_invivo = std(mean_r2_invivo, 0, 2);
%% 5.2 Plot fitting residual
avg_array2 = {'Ideal', 'Practical'};
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
plot(avg_x.', avg_compare.', 'LineStyle', 'none');

hold on;
ylim([0 60]);
ylim_lb = min(ylim); ylim_ub = max(ylim);
% patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
% patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
% patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
% patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)

avg_temp = mean_mean_res_avg16;
avg_sd_temp = std_res_avg16;
%avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,3) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{1}); hold on;
plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', '-o',  'LineWidth', 2, 'Color', color_cell{1}); 
%plot(mean_mean_res_invivo);

avg_temp = mean_mean_res_invivo;
avg_sd_temp = std_res_invivo;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,4) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2}); hold on;
plotHandles(:,2) = errorbar(avg_x.', avg_compare.', avg_err.', '-o', 'LineWidth', 2, 'Color', color_cell{2}); 


set(gca, 'FontSize', 18);
grid on;
% xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
% xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub]);

% text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
% text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
% text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
% text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
%ylabel('Fitting Residual (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
legend(plotHandles(1,:), avg_array2, 'Location', 'northeast');


color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');

%% 5.3 Plot R2
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_mean_r2_avg16;
avg_sd_temp = std_r2_avg16;
%avg_sd_temp = sd_snr_remote_avg16_1d;
avg_x = reshape(1:length(avg_temp), [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';

avg_array2 = {'Ideal', 'Practical'};
d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
plot(avg_x.', avg_compare.', 'LineStyle', 'none');
%semilogy(avg_x.', avg_compare.')
%set(gca, 'YScale', 'log');

hold on;
ylim([0 1]);
ylim_lb = min(ylim); ylim_ub = max(ylim);
% patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
% patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
% patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
% patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)


xlim([0 x(5)]);ylim([ylim_lb, ylim_ub]);

set(gca, 'Xcolor', 'w', 'Ycolor', 'k')
set(gca, 'XTick', []);
%set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');



figure('Position', [100 0 800 400]);
plotHandles = zeros(4,length(avg_array));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
plot(avg_x.', avg_compare.', 'LineStyle', 'none');
%semilogy(avg_x.', avg_compare.')
%set(gca, 'YScale', 'log');

hold on;
ylim([0 1]);
ylim_lb = min(ylim); ylim_ub = max(ylim);

plotHandles(:,3) = plot(avg_x.', 1-avg_compare.', 'LineWidth', 2, 'Color', color_cell{1}); hold on;
plotHandles(:,1) = errorbar(avg_x.', 1-avg_compare.', avg_err.', '-o',  'LineWidth', 2, 'Color', color_cell{1}); 
%plot(mean_mean_res_invivo);

avg_temp = mean_mean_r2_invivo;
avg_sd_temp = std_r2_invivo;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';
plotHandles(:,4) = plot(avg_x.', 1-avg_compare.', 'LineWidth', 2, 'Color', color_cell{2}); hold on;
plotHandles(:,2) = errorbar(avg_x.', 1-avg_compare.', avg_err.', '-o', 'LineWidth', 2, 'Color', color_cell{2}); 

set(gca, 'FontSize', 18);
%grid on;
% xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
% xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub]);

% text(2,ylim_ub-200, 'Slice Thickness = 2 mm', 'FontSize', 16);
% text(9,ylim_ub-200, 'Slice Thickness = 4 mm', 'FontSize', 16);
% text(16,ylim_ub-200, 'Slice Thickness = 6 mm', 'FontSize', 16);
% text(23,ylim_ub-200, 'Slice Thickness = 8 mm', 'FontSize', 16);
%set(gca,'xticklabel',{[]}); 
%ylabel('Fitting Residual (A.U.)');
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
%legend(plotHandles(1,:), avg_array2, 'Location', 'northeast');


color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
%title(cat(2, 'TE = ', num2str(TE_avg16(i)), '/', num2str(TE_invivo(i)), ' ms'));
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca,'ydir','reverse','yscale','log');

% yticks_temp = get(gca,'yticklabel');
yticks_temp = [-3, -2, -1, 0];
yticks_temp2 = 1-10.^yticks_temp;

set(gca,'yticklabel', yticks_temp2);

set(gca, 'Xcolor', 'w', 'Ycolor', 'k')
set(gca, 'XTick', []);
%set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');

%% 6. Iron content -> after 3.1 
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
iron_remote = [37.524, 44.387, 31.759, 27.296, 32.443, 31.385, 36.850, 37.827, 44.101, 36.543]; % ug/g
iron_mi =     [53.031, 37.556, 30.755, 61.667, 42.751, 33.786, 39.168, 48.570, 52.152, 44.290]; % ug/g
iron_ratio =  iron_mi ./ iron_remote;
trans16_avg

figure('Position', [100 0 700 350]);
plot(trans16_avg, iron_remote, 'o');
hold on; plot(trans16_avg,iron_mi, 'o');

figure('Position', [100 0 700 350]);
plotHandles_iron(:,1) = plot(trans16_avg, iron_ratio, 'o');

X1 = [trans16_avg];
mdl1 = fitlm(X1, iron_ratio);
mdl1.Coefficients;
ci = coefCI(mdl1);
b1 = mdl1.Coefficients.Estimate;

fitted_avg16 = b1(1) + X1*b1(2);
hold on;
plotHandles_iron(:,2) = plot(trans16_avg, fitted_avg16, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
set(plotHandles_iron(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [237, 62, 47]/255, 'MarkerFaceColor' , [237, 62, 47]/255)
set(plotHandles_iron(:,2), 'LineWidth', 2, 'Color', [.3, .3, .3]);
xlim([0.03, 0.22]); ylim([0.6 2.4]);
set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'YTick', [0.6, 1, 1.4, 1.8, 2.2]);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.025 0.025]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off')

%% 6.2 Iron content -> after 3.1 -> plot hemo-width vs iron-content
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
iron_remote = [37.524, 44.387, 31.759, 27.296, 32.443, 31.385, 36.850, 37.827, 44.101, 36.543]; % ug/g
iron_mi =     [53.031, 37.556, 30.755, 61.667, 42.751, 33.786, 39.168, 48.570, 52.152, 44.290]; % ug/g
iron_ratio =  iron_mi ./ iron_remote;

hemo_width = [1.21, 0.53, 0.56, 2.40, 0.91, 0.59, 0.92, 0.59, 1.10, 0.69]';

figure('Position', [100 0 700 350]);
plotHandles_iron(:,1) = plot(hemo_width, iron_ratio, 'o');

X1 = [hemo_width];
mdl1 = fitlm(X1, iron_ratio);
mdl1.Coefficients;
ci = coefCI(mdl1);
b1 = mdl1.Coefficients.Estimate;

fitted_avg16 = b1(1) + X1*b1(2);
hold on;
plotHandles_iron(:,2) = plot(hemo_width, fitted_avg16, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
plotHandles_iron(:,3) = plot(hemo_width([4 8]), iron_ratio([4 8]), 'o');
set(plotHandles_iron(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 14, ...
    'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor' , [0, 0.4470, 0.7410]);
set(plotHandles_iron(:,2), 'LineWidth', 3, 'Color', [.3, .3, .3]);
set(plotHandles_iron(:,3), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 14, ...
    'MarkerEdgeColor', [237, 62, 47]/255, 'MarkerFaceColor' , [237, 62, 47]/255);

xlim([0.5, 2.5]); 
ylim([0.6 2.4]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'YTick', [0.6, 1, 1.4, 1.8, 2.2]);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.025 0.025]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off')
%% T2* values in the hemorrhage (no clear correlation)
% subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
mean_t2star_array = [26.2308, 30.7600, 23.0624, 26.6097, 23.1755, 26.4926, 25.2607, 25.2070, 22.1008, 27.3151];
sd_t2star_array = [4.2677, 2.4475, 6.1379, 3.1481, 4.3568, 4.3842, 2.2723, 4.0519, 3.2397, 5.1619];
figure(); plot(mean_t2star_array, iron_ratio, 'o');
figure(); plot(sd_t2star_array, iron_ratio, 'o');

figure(); plot(trans16_avg, mean_t2star_array, 'o');
figure(); plot(trans16_avg, sd_t2star_array, 'o');

% No clear correlation

%% T2* value in different voxel size
inplane_res = [0.3, 0.6, 0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';
idx_mat = repmat(1:7, [4,1]) + 7*repmat([0:3]', [1,7]);

vol_array = vol_mat(:);
idx_array = idx_mat(:);

[B,I] = sort(vol_array)
BI = idx_array(I)

% 6, 11,12,17,18,23

AUC_avg16 = [1,0.96,0.94,0.92,0.86,0.8,0.71;0.92,0.88,0.88,0.85,0.8,0.76,0.69;0.83,0.81,0.82,0.78,0.77,0.74,0.7;0.76,0.75,0.74,0.74,0.71,0.67,0.67];
AUC_invivo = [nan,nan,0.57,0.68,0.68,0.77,0.64;nan,nan,0.69,0.76,0.76,0.75,0.66;nan,nan,0.76,0.76,0.71,0.69,0.6;nan,nan,0.74,0.69,0.66,0.66,0.58];

AUC_avg16_array = AUC_avg16(:);
AUC_invivo_array = AUC_invivo(:);

figure(); plot(B, AUC_avg16_array(I),'o'); grid on;
figure(); plot(B, AUC_invivo_array(I),'o'); grid on;

% load carbig
% tbl = table(Horsepower,Weight,MPG);
% modelfun = @(b,x)b(1) + b(2)*x(:,1).^b(3) + ...
%     b(4)*x(:,2).^b(5);
% beta0 = [-50 500 -1 500 -1];
% mdl = fitnlm(tbl,modelfun,beta0)

%% Flatten the voxel size, plot vs AUC
% fits the model specified by modelfun to variables in the table or dataset array tbl, and returns the nonlinear model mdl.
% fitnlm estimates model coefficients using an iterative procedure starting from the initial values in beta0.
figure('Position', [100 0 1000 300]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);

ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
plotHandles_auc(:,1) = plot(B, AUC_avg16_array(I),'o'); grid on;
hold on;
% plot(X, Y_pred); %ylim([0.5 1])
% Y_lb = modelfun(ci(:,1), X);
% Y_ub = modelfun(ci(:,2), X);
% plot(X, Y_lb); 
% plot(X, Y_ub); 
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [160 207 221]/255, 'MarkerFaceColor' , [160 207 221]/255)

plotHandles_auc(:,2) = plot(B, AUC_invivo_array(I),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [203 203 203]/255, 'MarkerFaceColor' , [203 203 203]/255)
xlim([0, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%% Flatten the voxel size, plot vs AUC (version 2 for SCMR)

figure('Position', [100 0 1000 300]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);

ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
plotHandles_auc(:,1) = plot(B, AUC_avg16_array(I),'o'); grid on;
hold on;
% plot(X, Y_pred); %ylim([0.5 1])
% Y_lb = modelfun(ci(:,1), X);
% Y_ub = modelfun(ci(:,2), X);
% plot(X, Y_lb); 
% plot(X, Y_ub); 
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor' , [.75 .75 1])

plotHandles_auc(:,2) = plot(B, AUC_invivo_array(I),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [253,190,133]/255, 'MarkerFaceColor' , [253,190,133]/255)
xlim([0, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%% Find optimal range for AUC 
AUC_invivo_array_reorder = AUC_invivo_array(I);
B_new = B(~isnan(AUC_invivo_array_reorder));
AUC_invivo_array_reorder_new = AUC_invivo_array_reorder(~isnan(AUC_invivo_array_reorder));

figure();
plot(B_new, AUC_invivo_array_reorder_new,'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [253,190,133]/255, 'MarkerFaceColor' , [253,190,133]/255);

x = 1:1:20;
for_K_means = [x', AUC_invivo_array_reorder_new];
idx = kmeans(for_K_means, 3);

s_array = zeros(9, 1);
for k = 2:10
    clust = kmeans(for_K_means, k);
    %clust(isnan(clust)) = [];
    s = silhouette(for_K_means,clust);
    s(isnan(s)) = [];
    s_array(k-1) = mean(s);
end

figure();
bar([2:10], s_array);
xlabel('Cluster #');
ylabel('Sihouette Coefficient');