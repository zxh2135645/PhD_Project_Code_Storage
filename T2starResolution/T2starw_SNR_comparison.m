clear all;
close all;
addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
% SNR from different avg from T2star_SNR_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

snr_glob = glob(cat(2, subject_data_dir, 'SNR_*'));
[list_to_read, order_to_read] = NamePicker(snr_glob, 1);

%% Read SNR metrics
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = load(f);
end

%% Compare
avg_array = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    strings = strsplit(list_to_read{order_to_read(i)},'/');
    strings2 = strsplit(strings{end},'_');
    strings3 = strsplit(strings2{end}, '.');
    avg_array{i} = strings3{1};
end

%% SNR in the air
figure('Position', [100 0 1600 1600]);
for i = 1:size(whatsinit{1}.SNR.snr_air, 2)
    avg_truth = whatsinit{1}.SNR.snr_air(:,i);
    avg_compare = whatsinit{2}.SNR.snr_air(:,i);
    subplot(2,3,i);
    plot(avg_truth, 'LineWidth', 2);hold on;
    plot(avg_compare, 'LineWidth', 2);
    legend({avg_array{1}, avg_array{2}}, 'Location', 'northwest');
    set(gca, 'FontSize', 18);
    grid on;
end
%% SNR in 3D plot
clear p;
figure('Position', [100 0 1600 1600]);
avg_truth = whatsinit{1}.SNR.snr_air;
%avg_compare = whatsinit{2}.SNR.snr_air;
snr_air_reshape = permute(reshape(avg_truth, 7, 4, []), [2,1,3]);
%snr_remote_reshape = permute(reshape(whatsinit{2}.SNR.snr_air, 7, 4, []), [2,1,3]);
X = 1:size(snr_air_reshape, 2);
color_array = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]];



for i = 1:size(whatsinit{1}.SNR.snr_air, 2)
    subplot(2,3,i);
    
    for k = 1:length(whatsinit)
        snr_air_reshape = permute(reshape(whatsinit{k}.SNR.snr_air, 7, 4, []), [2,1,3]);
        for j = 1:size(snr_air_reshape, 1)
            Y = j*ones(1, size(snr_air_reshape, 2));
            p{j,k} = plot3(X,Y,snr_air_reshape(j,:,i), 'LineWidth', 2, 'Color', color_array(k, :));hold on;
            %plot(avg_compare, 'LineWidth', 2);
            %
            set(gca, 'FontSize', 18);
            grid on;
        end
    end
    %xlabel('In-plane (mm^2)');ylabel('Through-plane (mm)');
    zlabel('SNR (unitless)');
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    h = [p{1,1};p{1,2}];
    legend(h, {avg_array{1},avg_array{2}}, 'Location', 'northwest');
    title(cat(2,'TE', num2str(i)));
end
%% This is not visually appealing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Position', [100 0 1600 1600]);
%num_slice = size(whatsinit{1}.SNR.snr_air, 1);
%num_te = size(whatsinit{1}.SNR.snr_air, 2);
%[X,Y] = meshgrid(1:num_te, 1:num_slice);

%avg_truth = whatsinit{1}.SNR.snr_air;
%avg_compare = whatsinit{2}.SNR.snr_air;
%s1 = surf(X,Y,avg_truth, 'FaceAlpha',0.5);hold on;
%s1.EdgeColor = 'none';
%s2 = surf(X,Y,avg_compare,'FaceAlpha',0.5);
%s2.EdgeColor = 'none';
%legend({avg_array{1}, avg_array{2}}, 'Location', 'northwest');
%set(gca, 'FontSize', 18);
%grid on;
%% SNR in remote
figure('Position', [100 0 1600 1600]);
for i = 1:size(whatsinit{1}.SNR.snr_remote, 2)
    avg_truth = whatsinit{1}.SNR.snr_remote(:,i);
    avg_compare = whatsinit{2}.SNR.snr_remote(:,i);
    subplot(2,3,i);
    plot(avg_truth, 'LineWidth', 2);hold on;
    plot(avg_compare, 'LineWidth', 2);
    legend({avg_array{1}, avg_array{2}}, 'Location', 'northwest');
    set(gca, 'FontSize', 18);
    grid on;
end

%% SNR in 3D plot
clear p;
figure('Position', [100 0 1600 1600]);
avg_truth = whatsinit{1}.SNR.snr_remote;
%avg_compare = whatsinit{2}.SNR.snr_remote;
snr_remote_reshape = permute(reshape(avg_truth, 7, 4, []), [2,1,3]);
%snr_remote_reshape = permute(reshape(whatsinit{2}.SNR.snr_air, 7, 4, []), [2,1,3]);
X = 1:size(snr_remote_reshape, 2);
color_array = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]];



for i = 1:size(whatsinit{1}.SNR.snr_remote, 2)
    subplot(2,3,i);
    
    for k = 1:length(whatsinit)
        snr_remote_reshape = permute(reshape(whatsinit{k}.SNR.snr_remote, 7, 4, []), [2,1,3]);
        for j = 1:size(snr_remote_reshape, 1)
            Y = j*ones(1, size(snr_remote_reshape, 2));
            p{j,k} = plot3(X,Y,snr_remote_reshape(j,:,i), 'LineWidth', 2, 'Color', color_array(k, :));hold on;
            %plot(avg_compare, 'LineWidth', 2);
            %
            set(gca, 'FontSize', 18);
            grid on;
        end
    end
    %xlabel('In-plane (mm^2)');ylabel('Through-plane (mm)');
    zlabel('SNR (unitless)');
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    h = [p{1,1};p{1,2}];
    legend(h, {avg_array{1},avg_array{2}}, 'Location', 'northwest');
    title(cat(2,'TE', num2str(i)));
end