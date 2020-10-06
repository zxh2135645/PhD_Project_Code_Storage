%% Tools for reading Inversion times from a dataset
close all;
clear all;

addpath('../function/');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '/*'));

Names = ExtractNames(folder_glob);

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

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
    'InversionTime',...
    };

proj_dir = GetFullPath(cat(2, pwd, '/../../Khalid_Patient_TI/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

%% Read T1 MOLLI weighted data
kh_ti = [219.8, 287.5, 909.8, 1007.2, 1600.9, 1733.5, 2460, 3178.2];
rms_cell = cell(length(Names), 1);

%% RMS to every slice for each patient
name_check = 'KIM_KI_BOK';
starting_point = find(strcmp(name_check, Names),1);

for n = starting_point:length(Names)
    f_glob = glob(cat(2, folder_glob{n}, 'T1w_MOCO\*'));
    if isempty(f_glob)
        disp('Folder is empty: ');
        disp(Names{n})
    else
        whatsinit = cell(length(f_glob), 1);
        slice_data = cell(length(f_glob), 1);
        for i = 1:length(f_glob)
            [whatsinit{i}, slice_data{i}] = dicom23D(f_glob{i}, dicom_fields);
        end
        
        if length(slice_data{1}) ~= 8
            disp('Number of Inversion Time is shorter than 8: ');
            disp(Names{n});
        else
            IR_array = zeros(length(kh_ti), length(f_glob));
            for i = 1:size(IR_array, 1)
                for j = 1:length(f_glob)
                    IR_array(i,j) = slice_data{j}(i).InversionTime;
                end
            end
            
            IR_array_sorted = sort(IR_array);
            rms = sqrt(sum((IR_array_sorted - kh_ti').^2));
            rms_cell{n} = rms;
        end
    end
end

%% find smallest rms
thresh = 200;
idx_cell = cell(length(f_glob), 1);
for i = 1:length(rms_cell)
    if ~isempty(rms_cell{i})
        idx_cell{i} = find(rms_cell{i} < thresh);
    else
        disp('empty array');
    end
end

%% Find corresponding name idex
name_idx = [];

for i = 1:length(idx_cell)
   if length(idx_cell{i}) > 1
       name_idx = [name_idx, i];
   end
end

%% name_idx
%Names{name_idx}

%for i = 1:length(name_idx)

f_glob = glob(cat(2, folder_glob{name_idx(5)}, 'T1w_MOCO\*'));
whatsinit = cell(length(f_glob), 1);
slice_data = cell(length(f_glob), 1);
for i = 1:length(f_glob)
    [whatsinit{i}, slice_data{i}] = dicom23D(f_glob{i}, dicom_fields);
end
%end

IR_array = zeros(length(kh_ti), length(f_glob));
for i = 1:size(IR_array, 1)
    for j = 1:length(f_glob)
        IR_array(i,j) = slice_data{j}(i).InversionTime;
    end
end

IR_array_sorted = sort(IR_array);
%% SAVE DATA
TI_encoding = struct;
TI_encoding.Names = Names;
TI_encoding.idx_cell = idx_cell;
TI_encoding.name_idx = name_idx;
save_dir = cat(2, data_dir, 'TI_encoding.mat');
save(save_dir, 'TI_encoding');