close all;
clear all;
%% ImageJ
% Convert mat to dicom (T1 map) for Image J
[fid_file, fid_path] = uigetfile('*.mat');

load(strcat(fid_path, fid_file(1:19), 'LRT_T1MOLLI_Mappings_Seg15.mat'));
t1_map = map_to_save.t1_map;

figure();
for i = 1:size(t1_map, 3)
    subplot(4,4,i)
    imagesc(t1_map(:,:,i,1).*mask(:,:,i)); caxis([100 1000]); axis image;
end


dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));
label = 'LRT';
idx_array = contains(folder_glob, label);
[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    'TriggerTime',...
    'RepetitionTime',...
    'EchoTime',...
    'EchoTrainLength',
    };

whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end

save_dir = cat(2, fid_path, 'DICOM_T1/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

NumEcho = slice_data(1).EchoTrainLength;
NumSlc = length(slice_data) / NumEcho;
slc = 4;

for m = 1:sizes(5)
    %metadata = dicominfo(slice_data{1}(m).Filename);
    
    i = ceil(m/NumSlc);
    slc_virtual = m - (i-1) * NumSlc;
    
    idx = (slc_virtual - 1)*NumEcho + i;
    metadata = dicominfo(slice_data(idx).Filename);
    
    t1 = uint16(t1_map(:,:,slc,m));
    metadata.WindowCenter = 500;
    metadata.WindowWidth =  800;
    metadata.SmallestImagePixelValue = min(t1(:));
    metadata.LargestImagePixelValue = max(t1(:));
    
    fname = cat(2, save_dir, fid_file(1:19), 'Echo', num2str(1), '_Section', num2str(m), '_Slice', num2str(slc), '_ForImageJ', '.dcm');
    dicomwrite(t1, fname, metadata, 'CreateMode', 'copy');
end
