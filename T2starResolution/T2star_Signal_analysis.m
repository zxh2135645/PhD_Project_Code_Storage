close all;
clear all;

% the main body for T2* SNR analysis

addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR_<avg_name>.mat
% from T2star_SNR_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014

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
avg_num_cell = {'Avg0016', 'Avg0001', 'Invivo'};
subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));

sig_avg16 = struct;
sig_avg01 = struct;
sig_invivo = struct;

for i = 1:length(avg_num_cell)
    avg_name = avg_num_cell{i};
    if i == 1
        snr_avg16 = load(cat(2, subject_data_dir, 'SNR_', avg_name, '.mat'));
        sig_avg16.sig_remote_mean = snr_avg16.SNR.sig.sig_remote_mean;
        sig_avg16.sig_remote_sd = snr_avg16.SNR.sig.sig_remote_sd;
        sig_avg16.sig_air_mean = snr_avg16.SNR.sig.sig_air_mean;
        sig_avg16.sig_air_sd = snr_avg16.SNR.sig.sig_air_sd;
    elseif i == 2
        snr_avg01 = load(cat(2, subject_data_dir, 'SNR_', avg_name, '.mat'));
        sig_avg01.sig_remote_mean = snr_avg01.SNR.sig.sig_remote_mean;
        sig_avg01.sig_remote_sd = snr_avg01.SNR.sig.sig_remote_sd;
        sig_avg01.sig_air_mean = snr_avg01.SNR.sig.sig_air_mean;
        sig_avg01.sig_air_sd = snr_avg01.SNR.sig.sig_air_sd;
    elseif i == 3
        snr_invivo = load(cat(2, subject_data_dir, 'SNR_', avg_name, '.mat'));
        sig_invivo.sig_remote_mean = snr_invivo.SNR.sig.sig_remote_mean;
        sig_invivo.sig_remote_sd = snr_invivo.SNR.sig.sig_remote_sd;
        sig_invivo.sig_air_mean = snr_invivo.SNR.sig.sig_air_mean;
        sig_invivo.sig_air_sd = snr_invivo.SNR.sig.sig_air_sd;
    end
end

%% Read SNR metrics
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

vx = zeros(length(slice_data), 1);
for i = 1:length(slice_data)
    x = slice_data{i}.PixelSpacing(1);
    y = slice_data{i}.PixelSpacing(2);
    z = slice_data{i}.SliceThickness;
    vx(i) = x*y*z;
end

%% Plot signals
[vx_sorted, I] = sort(vx);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
res_array = {'03', '06', '08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
lg = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    q = fix((i-1)/7)+1;
    r = mod(i, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    lg{i} = cat(2, slc_array{q}, '__', res_array{real_r});
    
end
sig_invivo_remote_mean = sig_invivo.sig_remote_mean;
sig_invivo_remote_sd = sig_invivo.sig_remote_sd;
sig_invivo_air_mean = sig_invivo.sig_air_mean;
sig_invivo_air_sd = sig_invivo.sig_air_sd;

sig_remote_mean_expand = zeros(length(sig_avg16.sig_remote_mean), 5);
sig_remote_sd_expand = zeros(length(sig_avg16.sig_remote_sd), 5);
sig_air_mean_expand = zeros(length(sig_avg16.sig_air_mean), 5);
sig_air_sd_expand = zeros(length(sig_avg16.sig_air_sd), 5);

for i = 1:length(sig_invivo_remote_mean)
    mask_idx = mask_idx_array(i);
    sig_remote_mean_expand(mask_idx,:) = sig_invivo_remote_mean(i,:);
    sig_remote_sd_expand(mask_idx,:) = sig_invivo_remote_sd(i,:);
    sig_air_mean_expand(mask_idx,:) = sig_invivo_air_mean(i,:);
    sig_air_sd_expand(mask_idx,:) = sig_invivo_air_sd(i,:);
end

sig_remote_mean_expand(sig_remote_mean_expand == 0) = nan;
sig_remote_sd_expand(sig_remote_sd_expand == 0) = nan;
sig_air_mean_expand(sig_air_mean_expand == 0) = nan;
sig_air_sd_expand(sig_air_sd_expand == 0) = nan;

figure('Position', [100 0 1600 1600]);
for i = 1:size(sig_invivo_remote_mean, 2)
    subplot(2,3,i);
    errorbar(sig_avg16.sig_remote_mean(:,i),sig_avg16.sig_remote_sd(:,i), 'LineWidth', 2);
    hold on;
    errorbar(sig_avg01.sig_remote_mean(:,i), sig_avg01.sig_remote_sd(:,i), 'LineWidth', 2);
    errorbar(sig_remote_mean_expand(:,i),sig_remote_sd_expand(:,i), 'LineWidth', 2);
    legend(avg_num_cell, 'Location', 'NorthEast');
    ylim([0 3000]);grid on;
end

figure('Position', [100 0 1600 1600]);
for i = 1:size(sig_invivo_remote_mean, 2)
    subplot(2,3,i);
    plot(sig_avg16.sig_remote_mean(:,i), 'LineWidth', 2);
    hold on;
    plot(sig_avg01.sig_remote_mean(:,i), 'LineWidth', 2);
    plot(sig_remote_mean_expand(:,i), 'LineWidth', 2);
    legend(avg_num_cell, 'Location', 'NorthEast');
    ylim([0 3000]);grid on;
end

figure('Position', [100 0 1600 1600]);
for i = 1:size(sig_invivo_remote_sd, 2)
    subplot(2,3,i);
    plot(sig_avg16.sig_remote_sd(:,i), 'LineWidth', 2);
    hold on;
    plot(sig_avg01.sig_remote_sd(:,i), 'LineWidth', 2);
    plot(sig_remote_sd_expand(:,i), 'LineWidth', 2);
    legend(avg_num_cell, 'Location', 'NorthEast');
    grid on;
end
%% Air
figure('Position', [100 0 1600 1600]);
for i = 1:size(sig_invivo_air_mean, 2)
    subplot(2,3,i);
    plot(sig_avg16.sig_air_mean(:,i), 'LineWidth', 2);
    hold on;
    plot(sig_avg01.sig_air_mean(:,i), 'LineWidth', 2);
    plot(sig_air_mean_expand(:,i), 'LineWidth', 2);
    legend(avg_num_cell, 'Location', 'NorthEast');
    ylim([0 200]);grid on;
end

figure('Position', [100 0 1600 1600]);
for i = 1:size(sig_invivo_air_sd, 2)
    subplot(2,3,i);
    plot(sig_avg16.sig_air_sd(:,i), 'LineWidth', 2);
    hold on;
    plot(sig_avg01.sig_air_sd(:,i), 'LineWidth', 2);
    plot(sig_air_sd_expand(:,i), 'LineWidth', 2);
    legend(avg_num_cell, 'Location', 'NorthEast');
    ylim([0 100]);grid on;
end