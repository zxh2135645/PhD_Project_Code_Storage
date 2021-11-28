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

%% Load LRT recon
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
%% Display image
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);

temp = Gr\reshape(Phi(:,31,1,1,:), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));

figure();
temp = imrotate(temp, 90);
temp = abs(flip(temp(:,:,1),2)/cw);
temp = uint16(temp*4095);
ax2 = imagesc(temp); axis image; colormap gray;axis off;
%% 6TEs and then next sliceloc
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

save_dir = cat(2, fid_path, 'DICOM_E/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


echo_f_glob = glob(cat(2, fid_path, '*Echo*.mat'));

slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
slc_array = [9 8 7 6 5 4 3 2 1 16 15 14 13 12 11 10];

num_seg = 31;
for te = 1:NEco_old
    f = echo_f_glob{te};
    load(f, 'Gr', 'Phi', 'L', 'U', 'Nx', 'Ny', 'Nz', 'params', 'sizes');    
    for i = 1:len
        slc = slc_array(i);
        dispim = @(x)fftshift(x(:,:,slc,:),1);
        temp = Gr\reshape(Phi(:,num_seg,1,1,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
        cw = max(vec(abs(temp)));
        
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(te), '_Seg', num2str(num_seg), '_Slc', num2str(i), '_E', '.dcm');
        
        temp = imrotate(temp, 90);
        temp = abs(flip(temp(:,:,1),2)/cw);
        temp = uint16(temp*4095);
        
        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    end
end
