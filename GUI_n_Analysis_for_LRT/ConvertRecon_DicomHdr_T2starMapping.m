clear all;
close all;
%% Read DICOM file
addpath('../function/');

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'LRT'};

label = labels{1};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

%% Load LRT recon - T2star mapping
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file));
%% Display image
slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
figure();
for i = 1:size(t2star,3)
    slc = slc_array(i);
    temp = imrotate(t2star(:,:,slc), 90);
    temp = abs(flip(temp,2));
    temp = uint16(temp);
    subplot(4,4,i); imagesc(temp); axis image; axis off;
    colormap gray;
end

%% 
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

save_dir = GetFullPath(cat(2, fid_path, '../DICOM_T2STAR_L/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for i = 1:len
    
    metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
    fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_Slc', num2str(i), '_L', '.dcm'));
    
    slc = slc_array(i);
    temp = imrotate(t2star(:,:,slc), 90);
    temp = abs(flip(temp,2));
    temp = uint16(temp);

    dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
end
