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
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%% Display image (Write DICOM)
slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
slc_array = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
num_seg_array = [16, 21, 26, 31, 36, 41];
% num_seg_array = [141, 151, 161, 171, 181, 191];
temtemp_4D = zeros(Ny, Nx, length(num_seg_array), length(slc_array));

% cardiac phase and resp phase needs to be encoded
resp_phase = 1;
card_phase = 18;

for i = 1:length(slc_array)
    slc = slc_array(i);
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    for num = 1:length(num_seg_array)
        temp = Gr\reshape(Phi(:,num_seg_array(num),card_phase,resp_phase,end), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);

        temtemp_4D(:,:,num,i) = temp;
    end
end

temp_4D = permute(temtemp_4D, [1 2 4 3]);

% Write LGE into DICOM
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

save_dir = GetFullPath(cat(2, fid_path, 'DICOM_LGE_CVI/'));
% save_dir = GetFullPath(cat(2, fid_path, '/DICOM_T2star_CVI/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for te = 1:length(num_seg_array)
    num_seg = num_seg_array(te); 
    for i = 1:len
        %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_LGE', '.dcm'));
        %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));
        
        
        slc_array_flip = flip(slc_array);
        slc = slc_array_flip(i);
        temp = imrotate(temp_4D(:,:,slc,te), 90);
        cw = max(vec(abs(temp_4D(:,:,slc,te))));

        % temp = temp_4D(:,:,slc,te);
        temp = abs(flip(temp,2)./cw);
        temp = uint16(temp*4095);
        metadata.WindowCenter = 2048;
        metadata.WindowWidth = 4095;

        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
        metadata.SmallestImagePixelValue = min(temp(:));
        metadata.LargestImagePixelValue = max(temp(:));
        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    end
end