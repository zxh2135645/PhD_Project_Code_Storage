clear all;
close all;
%% Read DICOM file
addpath('../function/');

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));


dicom_files = {};
dicom_headers = {};
dicom_images = {};

for i = 1:length(folder_glob)
    files = dir(fullfile(folder_glob{i}, '*.dcm'));
    for j = 1:length(files)
        file_path = fullfile(folder_glob{i}, files(j).name);
        dicom_files{end+1} = file_path;
        dicom_headers{end+1} = dicominfo(file_path);
        dicom_images{end+1} = dicomread(file_path);
    end
end

%%
figure;
num_images = length(dicom_images);
cols = ceil(sqrt(num_images));
rows = ceil(num_images / cols);

for k = 1:num_images
    subplot(rows, cols, k);
    imshow(dicom_images{k}, []);
    title(sprintf('Image %d', k));
end

%%
% Extract SliceLocation from dicom_headers
slice_locations = zeros(1, length(dicom_headers));
for idx = 1:length(dicom_headers)
    if isfield(dicom_headers{idx}, 'SliceLocation')
        slice_locations(idx) = dicom_headers{idx}.SliceLocation;
    else
        slice_locations(idx) = NaN; % If SliceLocation is missing
    end
end
disp('SliceLocations:');
disp(slice_locations);

%% CINE
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));


dicom_files_cine = {};
dicom_headers_cine = {};
dicom_images_cine = {};

for i = 1:length(folder_glob)
    files = dir(fullfile(folder_glob{i}, '*.dcm'));
    for j = 1:length(files)
        file_path = fullfile(folder_glob{i}, files(j).name);
        dicom_files_cine{end+1} = file_path;
        dicom_headers_cine{end+1} = dicominfo(file_path);
        dicom_images_cine{end+1} = dicomread(file_path);
    end
end

%% Compare the difference between dicom_headers_cine{1} and dicom_headers{1}
fields_cine = fieldnames(dicom_headers_cine{1});
fields_static = fieldnames(dicom_headers{1});
all_fields = unique([fields_cine; fields_static]);

fprintf('Differences between dicom_headers_cine{1} and dicom_headers{1}:\n');
for i = 1:length(all_fields)
    field = all_fields{i};
    has_cine = isfield(dicom_headers_cine{1}, field);
    has_static = isfield(dicom_headers{1}, field);
    if has_cine && has_static
        val_cine = dicom_headers_cine{1}.(field);
        val_static = dicom_headers{1}.(field);
        if ~isequal(val_cine, val_static)
            %fprintf('Field "%s" differs:\n', field);
            %disp('  cine:');
            %disp(val_cine);
            %disp('  static:');
            %disp(val_static);
        end
    elseif has_cine
        fprintf('Field "%s" only in dicom_headers_cine{1}\n', field);
    elseif has_static
        fprintf('Field "%s" only in dicom_headers{1}\n', field);
    end
end

%% Load LRT recon 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%%

% Need to assign 
TriggerTime = 'TriggerTime';
NominalInterval = 'NominalInterval';
NumberOfTemporalPositions = 'NumberOfTemporalPositions';
InstanceNumber = 'InstanceNumber';
AcquisitionNumber = 'AcquisitionNumber';
SpacingBetweenSlices = 'SpacingBetweenSlices';

% Can just copy
% 'IntervalsAcquired'
% 'IntervalsRejected'
% 'LowRRValue'
% 'HighRRValue'
% 'ImageComments'

% Copy selected fields from dicom_headers_cine{1} to dicom_headers{1}
fields_to_copy = {'IntervalsAcquired', 'IntervalsRejected', 'LowRRValue', 'HighRRValue', 'ImageComments'};
dicom_headers_to_copy = cell(1, length(dicom_headers)*size(Phi,3));

for slc = 1:length(dicom_headers)
    for frame = 1:size(Phi,3)
        idx = (slc-1)*size(Phi,3) + frame;
        dicom_headers_to_copy{idx} = dicom_headers{slc};
        dicom_headers_to_copy{idx}.(TriggerTime) = (frame-1) * RR_int / size(Phi,3);
        dicom_headers_to_copy{idx}.(NominalInterval) = RR_int;
        dicom_headers_to_copy{idx}.(NumberOfTemporalPositions) = size(Phi,3);
        dicom_headers_to_copy{idx}.(InstanceNumber) = idx;
        dicom_headers_to_copy{idx}.(AcquisitionNumber) = idx;
        dicom_headers_to_copy{idx}.(SpacingBetweenSlices) = 6;

        for i = 1:length(fields_to_copy)
            field = fields_to_copy{i};
            if isfield(dicom_headers_cine{1}, field)
                dicom_headers_to_copy{idx}.(field) = dicom_headers_cine{1}.(field);
            end
        end
    end
end
%% Display image (Write DICOM)
slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
slc_array = fftshift(1:Nz);
num_seg_array = [61];
% num_seg_array = [141, 151, 161, 171, 181, 191];
temtemp_4D = zeros(Ny, Nx, size(Phi,3), length(slc_array));

% cardiac phase and resp phase needs to be encoded
resp_phase = 1;
card_phase = 18;

for i = 1:length(slc_array)
    slc = slc_array(i);
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    for num = 1:length(num_seg_array)
        temp = Gr\reshape(Phi(:,num_seg_array(num),:,resp_phase,1), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);

        temtemp_4D(:,:,:,i) = temp;
    end
end

% Crop first dimension of temtemp_4D to dicom_headers_to_copy{1}.Rows, centered
target_rows = dicom_headers_to_copy{1}.Rows;
current_rows = size(temtemp_4D, 1);
row_start = floor((current_rows - target_rows)/2) + 1;
row_end = row_start + target_rows - 1;
temtemp_4D = temtemp_4D(row_start:row_end, :, :, :);
%temp_4D = permute(temtemp_4D, [1 2 4 3]);

%% Write LGE into DICOM
X = dicom_headers{1}.Filename;
NEco_old = params.NEco_old;
%len = length(slice_data{1})/NEco_old;

save_dir = GetFullPath(cat(2, fid_path, 'DICOM_CINE_V2/'));
% save_dir = GetFullPath(cat(2, fid_path, '/DICOM_T2star_CVI/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for te = 1:length(num_seg_array)
    num_seg = num_seg_array(te);
    for i = 1:Nz
        %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
        %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);

        f_dir = GetFullPath(cat(2, save_dir, 'Slice_', num2str(i,'%02.0f'), '/'));

        if ~exist(f_dir, 'dir')
            mkdir(f_dir);
        end

        for j = 1:size(Phi,3)
            fname = GetFullPath(cat(2, f_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i,'%02.0f'), '_CINE', num2str(j,'%02.0f'), '.dcm'));
            %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));

            metadata = dicom_headers_to_copy{(i-1)*size(Phi,3) + j};

            %slc_array_flip = flip(slc_array);
            %slc = slc_array_flip(i);
            temp = imrotate(temtemp_4D(:,:,j,i), 90);
            cw = 0.8*max(vec(abs(temtemp_4D(:,:,j,i))));

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
end