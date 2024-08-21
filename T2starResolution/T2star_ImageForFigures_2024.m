clear all;
close all;

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
%% 5.1 Fitting Residual and R2
avg_array = {'Avg0016', 'Invivo'};
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Invivo'};
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
        % temp(temp<=0) = 1;
        temp(temp<=0) = 0;
        mean_r2_avg16(j, i) = mean(temp);
    end
    for j = 1:length(mean_r2_invivo)
        temp_invivo = FitResults_invivo.FitResults_struct(j).FitResults.R2(~isnan(FitResults_invivo.FitResults_struct(j).FitResults.R2));
        % temp_invivo(temp_invivo<=0) = 1;
        temp_invivo(temp_invivo<=0) = 0;
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

%% Load Real invivo data
FitResults_realinivo = load(cat(2, data_dir, 'Invivo_Fitting/FitResults_Invivo.mat'));
mean_res_realinvivo = zeros(1, length(subject_name_cell));
mean_r2_realinvivo = zeros(1, length(subject_name_cell));
% 17P17, 18P90, 18P92, 18P93, 18P94, 18P95, 20P10C, 20P11, 20P40, 20P48
for j = 1:size(mean_res_realinvivo, 2)
    mean_res_realinvivo(1, j) = mean(mean(FitResults_realinivo.FitResults_struct(j).FitResults.res(~isnan(FitResults_realinivo.FitResults_struct(j).FitResults.res))));
end

for j = 1:size(mean_r2_realinvivo, 2)
    temp = FitResults_realinivo.FitResults_struct(j).FitResults.R2(~isnan(FitResults_realinivo.FitResults_struct(j).FitResults.R2));
    temp(temp<=0) = 0; % why it's not 0???
    mean_r2_realinvivo(1, j) = mean(temp);
end

mean_mean_res_realinvivo = mean(mean_res_realinvivo, 2);
std_res_realinvivo = std(mean_res_realinvivo, 0, 2);

mean_mean_r2_realinvivo = mean(mean_r2_realinvivo, 2);
std_r2_realinvivo = std(mean_r2_realinvivo, 0, 2);
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
plotHandles(:,2) = errorbar(avg_x.', 1-avg_compare.', avg_err.', '-s', 'LineWidth', 2, 'Color', color_cell{2}); 

avg_temp = mean_mean_r2_realinvivo;
avg_sd_temp = std_r2_realinvivo;
avg_x = mask_idx_array(14);
avg_compare = avg_temp;
avg_err = avg_sd_temp;
plotHandles(:,6) = plot(avg_x.', 1-avg_compare.', 'LineWidth', 2, 'Color', color_cell{2}); hold on;
plotHandles(:,5) = errorbar(avg_x.', 1-avg_compare.', avg_err.', '-^', 'LineWidth', 2, 'Color', color_cell{2}); 


set(gca, 'FontSize', 18);
%grid on;
% xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
% xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub]);


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

color_cell_invivo = {[242,240,247]/255, [203,201,226]/255, [158,154,200]/255, [117,107,177]/255, [84,39,143]/255}; % Purple
%color_cell_invivo = {[255,255,204]/255, [194,230,153]/255, [120,198,121]/255, [49,163,84]/255, [0,104,55]/255};

set(plotHandles(:,6), 'LineStyle', 'none', 'Marker', '^', 'Color', color_cell_invivo{4});
set(plotHandles(:,6), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});
set(plotHandles(:,5), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca,'ydir','reverse','yscale','log');
%set(gca,'yscale','log');


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
legend({'Ideal', 'Ex-vivo Practical', 'In-vivo Practical'});