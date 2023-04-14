% This code is specifically for fresh heart scans
% Jesse
% Quantitative maps including: T2* CMR, T2-SE, T1-IR-SE
% First part of pipeline for 
% Fresh heart analysis:
% 1. FreshHeart_08072022.m
% 2. 
clear all;
close all;

% T2-SE
addpath('../function/')
qMRinfo('mono_t2'); % set it up first in qMRLab
dicom_dir = uigetdir;
dicom_glob = glob(cat(2, dicom_dir, '/T1_SE_TE*'));
dicom_glob_reorder = dicom_glob([2,3,4,6,7,1]); % JESSE
dicom_glob_reorder = dicom_glob([2,3,4,5,6,7,1]); % CHILI, NUTMEG, PAPRIKA, CINNAMON

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

whatsinit = cell(length(dicom_glob_reorder), 1);
slice_data = cell(length(dicom_glob_reorder), 1);
TE_array = zeros(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
    TE_array(i) = slice_data{i}.EchoTime;
end



%% Draw Mask (rect)
% sizes(2) -> T1
% sizes(3) -> Cardiac
% sizes(4) -> resp

mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_rect.mat'));
mask = zeros(size(whatsinit{1},1), size(whatsinit{1},2), size(whatsinit{1},3));
if ~exist(mask_f)
    for i = 1:size(whatsinit{1},3)

        figure();
        imagesc(whatsinit{1}); axis image;
        roi = drawpolygon;
        mask(:,:,i) = createMask(roi);
    end
    save(mask_f, 'mask');
else
    load(mask_f);
end

% mask = ones(size(whatsinit{1},1), size(whatsinit{1},2), size(whatsinit{1},3));

%% T2 mapping
NEco = length(TE_array); % 6
% for i = 1:Nz
addpath('../function/');
Nx = size(whatsinit{1},1);
Ny = size(whatsinit{1},2);
Nz = size(whatsinit{1},3);

t2_map = zeros([Nx, Ny, Nz]);

% i = input(sprintf('Select Slice of Interest [%d]: ', 3));
%
for slc = 1:Nz
    mask_temp = mask(:,:,slc);
    if any(mask_temp(:))
        for neco = 1:NEco
            if neco == 1
                temp_4D = zeros([Nx, Ny, Nz, NEco]);
            end
            temp_4D(:,:,:,neco) = whatsinit{neco};
        end


        % Reshape matrix as [Width x Height x #Slice x #TE]
        ipt = abs(temp_4D(:,:,:,:));

        Model = mono_t2;  % Create class from model
        %Model = Custom_OptionsGUI(Model);
        Model.Prot.SEdata.Mat = TE_array; %
        Model.st = [50 2000];
        Model.lb = [1 2000];
        Model.fx = [0 0];
        Model.voxelwise = 1;
        %Model.options.FitType = 'Linear';
        data = struct;  % Create data structure
        data.SEdata = ipt;
        data.Mask = squeeze(mask(:,:,slc));
        FitResults = FitData(data, Model); %fit data
        t2_map(:,:,slc) = FitResults.T2;
    end
end

figure();
imagesc(t2_map); axis image; caxis([0 100]); axis off;
colormap(brewermap([],'RdYlBu')); colorbar;

%% T1-IR-SE
dicom_glob = glob(cat(2, dicom_dir, '/IR_T1_TSE_TI*'));
dicom_glob_reorder = dicom_glob([7,3,5,6,1,2,4]); % JESSE, GINGER, PAPRIKA
%dicom_glob_reorder = dicom_glob([9,3,5,6,1,2,4]); % CHILI
%dicom_glob_reorder = dicom_glob([3,5,6,1,2,4]); % NUTMEG -> first TI has different resolution


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
    'EchoTrainLength',...
    'InversionTime'
    };

whatsinit = cell(length(dicom_glob_reorder), 1);
slice_data = cell(length(dicom_glob_reorder), 1);
IR_array = zeros(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
    IR_array(i) = slice_data{i}.InversionTime;
end

Nx = size(whatsinit{1},1);
Ny = size(whatsinit{1},2);
Nz = size(whatsinit{1},3);
t1_map = zeros(Nx, Ny, Nz);
%i = input(sprintf('Select Slice of Interest [%d]: ', 3));
%for i = 1:Nz

N_inv = length(IR_array);
for i = 1:Nz
    mask_temp = mask(:,:,i);
    if any(mask_temp(:))
        for neco = 1:N_inv
            if neco == 1
                temp_4D = zeros([Nx, Ny, Nz, N_inv]);
            end
            temp_4D(:,:,:,neco) = whatsinit{neco};
        end
        % a - create object
        Model = inversion_recovery;


        data = struct;
        data.IRData= double(temp_4D);
        data.Mask= double(mask(:,:,i));

        Model.Prot.IRData.Mat = IR_array;
        Model.voxelwise = 1;

        % b- fit dataset
        FitResults = FitData(data,Model,0);

        mask_temp = double(mask(:,:,i));

        % t1_map(:,:,i,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
        t1_map(:,:,i) = FitResults.T1;
    end
end

figure();
imagesc(t1_map); axis image; caxis([800 1600]); axis off;
colormap(brewermap([],'RdBu')); colorbar;
%% T2* mapping
dicom_glob = glob(cat(2, dicom_dir, '/T2_MULTIECHO_2D_SAX*'));
%dicom_glob_reorder = dicom_glob([3,4,5,6,7,8,9,10,11,1,2]); % JESSE
%dicom_glob_reorder = dicom_glob([2,3,4,5,6,7,8,9,11,1]);    % CHILI
%dicom_glob_reorder = dicom_glob([8,10,12,14,16,18,20,22,24,1,5]);    % NUTMEG
%dicom_glob_reorder = dicom_glob([1,6,8,10,12,14,16,18,20,22,4]); % GINGER
%dicom_glob_reorder = dicom_glob([3,5,7,9,11,13,15,17,19,1]); % PAPRIKA
dicom_glob_reorder = dicom_glob([3,6,8,11,13,15,17,19,22,1]); % CINNAMON

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

whatsinit = cell(length(dicom_glob_reorder), 1);
slice_data = cell(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end

TE_array = zeros(length(slice_data{1}), 1);
for i = 1:length(slice_data{1})
    TE_array(i) = slice_data{1}(i).EchoTime;
end

mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_rect_multiecho.mat'));
mask = zeros(size(whatsinit{1},1), size(whatsinit{1},2), length(dicom_glob_reorder));
if ~exist(mask_f)
    for i = 1:length(dicom_glob_reorder)

        figure();
        temp = whatsinit{i};
        imagesc(temp(:,:,4)); axis image;
        roi = drawpolygon;
        mask(:,:,i) = createMask(roi);
    end
    save(mask_f, 'mask');
else
    load(mask_f);
end


NEco = length(TE_array); % 8
Nx = size(whatsinit{1},1);
Ny = size(whatsinit{1},2);
Nz = length(whatsinit);

t2star_map = zeros([Nx, Ny, Nz]);

for slc = 1:Nz
    mask_temp = mask(:,:,slc);
    if any(mask_temp(:))

        temp_3D = zeros([Nx, Ny, NEco]);
        temp_3D = whatsinit{slc};
        temp_4D = reshape(temp_3D, Nx, Ny, 1, NEco);

        % Reshape matrix as [Width x Height x #Slice x #TE]
        ipt = abs(temp_4D(:,:,:,:));

        Model = mono_t2;  % Create class from model
        %Model = Custom_OptionsGUI(Model);
        Model.Prot.SEdata.Mat = TE_array; %
        Model.st = [50 2000];
        Model.lb = [1 2000];
        Model.fx = [0 0];
        Model.voxelwise = 1;
        %Model.options.FitType = 'Linear';
        data = struct;  % Create data structure
        data.SEdata = ipt;
        data.Mask = squeeze(mask(:,:,slc));
        FitResults = FitData(data, Model); %fit data
        t2star_map(:,:,slc) = FitResults.T2;
    end
end

%%
figure();
imagesc(t2star_map(:,:,8)); axis image; caxis([0 80]); colorbar;
colormap(brewermap([],'RdBu')); axis off;

%% Save
map_to_save = struct;
map_to_save.t1_map = t1_map;
map_to_save.t2_map = t2_map;
map_to_save.t2star_map = t2star_map;

save_dir = GetFullPath(cat(2, dicom_dir, '/../Result/QuantitativeMap.mat'));
save(save_dir, '-struct', 'map_to_save');


