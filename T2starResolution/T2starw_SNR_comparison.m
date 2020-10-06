clear all;
close all;
addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
% SNR from different avg from T2star_SNR_analysis_main.m
% Thus, need to run T2star_SNR_analysis_main.m before this script
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
if length(snr_glob) == 3
    % Avg0016, Avg0001, Invivo
    list_to_read = snr_glob;
    order_to_read = [2, 1, 3];
elseif length(snr_glob) == 2
    % Avg0016, Avg0001
    [list_to_read, order_to_read] = NamePicker(snr_glob, 1);
end

%% Read SNR metrics
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = load(f);
end

% Compare
avg_array = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    strings = strsplit(list_to_read{order_to_read(i)},'/');
    strings2 = strsplit(strings{end},'_');
    strings3 = strsplit(strings2{end}, '.');
    avg_array{i} = strings3{1};
end

%% SNR in the air
figure('Position', [100 0 1600 1600]);
plotHandles = zeros(4,length(whatsinit));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

for i = 1:size(whatsinit{1}.SNR.snr_air, 2)
    for m = 1:length(whatsinit)
        if m == 3
            avg_temp = whatsinit{m}.SNR.snr_air(:,i);
            avg_x = reshape(mask_idx_array, [], 4).';
            avg_compare = reshape(avg_temp, [], 4).';
        else
            avg_temp = whatsinit{m}.SNR.snr_air(:,i);
            avg_x = reshape(1:length(avg_temp), [], 4).';
            avg_compare = reshape(avg_temp, [], 4).';
        end
        
        subplot(2,3,i);
        plotHandles(:,m) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{m});hold on;
        set(gca, 'FontSize', 18);
        grid on;
        set(gca,'xticklabel',{[]}); ylabel('SNR (arb. unit.)');
    end
    legend(plotHandles(1,:), avg_array, 'Location', 'northwest');
end
%% SNR in 3D plot
clear p;
figure('Position', [100 0 1600 1600]);
%avg_truth = whatsinit{1}.SNR.snr_air;
%avg_compare = whatsinit{2}.SNR.snr_air;
%snr_air_reshape = permute(reshape(avg_truth, 7, 4, []), [2,1,3]);
%snr_remote_reshape = permute(reshape(whatsinit{2}.SNR.snr_air, 7, 4, []), [2,1,3]);
%X = 1:size(snr_air_reshape, 2);

color_array = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290, 0.6940, 0.1250]];
row = 4;

for i = 1:size(whatsinit{1}.SNR.snr_air, 2)
    subplot(2,3,i);
    h = [];
    for k = 1:length(whatsinit)
        col = length(whatsinit{k}.SNR.snr_air) / row;
        snr_air_reshape = permute(reshape(whatsinit{k}.SNR.snr_air, col, row, []), [2,1,3]);
        if k == 3
            X = [1:size(snr_air_reshape, 2)]+2;
        else
            X = 1:size(snr_air_reshape, 2);
        end
        for j = 1:size(snr_air_reshape, 1)
            
            
            Y = j*ones(1, size(snr_air_reshape, 2));
            p{j,k} = plot3(X,Y,snr_air_reshape(j,:,i), 'LineWidth', 2, 'Color', color_array(k, :));hold on;
            %plot(avg_compare, 'LineWidth', 2);
            set(gca, 'FontSize', 18);
            grid on;
        end
    h = [h; p{1,k}];
    end
    %xlabel('In-plane (mm^2)');ylabel('Through-plane (mm)');
    zlabel('SNR (unitless)');
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    
    legend(h, avg_array, 'Location', 'northwest');
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
plotHandles = zeros(4,length(whatsinit));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

for i = 1:size(whatsinit{1}.SNR.snr_remote, 2)
    for m = 1:length(whatsinit)
        if m == 3
            avg_temp = whatsinit{m}.SNR.snr_remote(:,i);
            avg_x = reshape(mask_idx_array, [], 4).';
            avg_compare = reshape(avg_temp, [], 4).';
        else
            avg_temp = whatsinit{m}.SNR.snr_remote(:,i);
            avg_x = reshape(1:length(avg_temp), [], 4).';
            avg_compare = reshape(avg_temp, [], 4).';
        end
        
        subplot(2,3,i);
        plotHandles(:,m) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{m});hold on;
        set(gca, 'FontSize', 18);
        grid on;
        set(gca,'xticklabel',{[]}); ylabel('SNR (arb. unit.)');
    end
    legend(plotHandles(1,:), avg_array, 'Location', 'southeast');
end

%% SNR in 3D plot
clear p;
figure('Position', [100 0 1600 1600]);

color_array = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290, 0.6940, 0.1250]];
row = 4;

for i = 1:size(whatsinit{1}.SNR.snr_remote, 2)
    subplot(2,3,i);
    h = [];
    for k = 1:length(whatsinit)
        col = length(whatsinit{k}.SNR.snr_remote) / row;
        snr_remote_reshape = permute(reshape(whatsinit{k}.SNR.snr_remote, col, row, []), [2,1,3]);
        if k == 3
            X = [1:size(snr_remote_reshape, 2)]+2;
        else
            X = 1:size(snr_remote_reshape, 2);
        end
        for j = 1:size(snr_remote_reshape, 1)
            
            
            Y = j*ones(1, size(snr_remote_reshape, 2));
            p{j,k} = plot3(X,Y,snr_remote_reshape(j,:,i), 'LineWidth', 2, 'Color', color_array(k, :));hold on;
            %plot(avg_compare, 'LineWidth', 2);
            set(gca, 'FontSize', 18);
            grid on;
        end
    h = [h; p{1,k}];
    end
    %xlabel('In-plane (mm^2)');ylabel('Through-plane (mm)');
    zlabel('SNR (unitless)');
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    
    legend(h, avg_array, 'Location', 'northwest');
    title(cat(2,'TE', num2str(i)));
end