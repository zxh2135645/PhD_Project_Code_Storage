clear all;
close all;

% the main body for T2* resolution analysis
% integration of ../LineWidth_Analysis.m

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% aha_anlysis
% T2star_meanSD_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014
dicom_dir = uigetdir;
%base_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob);

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

% disp('Avg 0016 starts here: ');
% avg_num = input('Please type average number here:  ');
% avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

%% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end

%%
save_array = 1:1:length(whatsinit);
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    figure();
    img2 = whatsinit{save_idx};
    %thresh = mean(nonzeros(img2 .* mask_struct(save_idx).remote_mask)) - 2*std(nonzeros(img2 .* mask_struct(save_idx).remote_mask));
    %hemo_mask = img2 < thresh;
    %subplot(1,2,1);
    imagesc(img2); caxis([0 100]); axis image; %colormap default; %colorbar;
    colormap(brewermap([],'RdBu'));
    
    mean2sd_dir = cat(2, subject_dir, 'Exvivo_Avg0016/');
    if ~exist(mean2sd_dir, 'dir')
        mkdir(mean2sd_dir)
    end
    saveas(gcf, cat(2, mean2sd_dir, num2str(save_idx), '.png'));
end