% This code is specifically for ex-vivo heart scans in formalin
% Carlos, George, Tony, Roger
% Quantitative maps including: T2* mGRE, T2-TSE, T1-IRTSE
% First part of pipeline for 
% Fresh heart analysis:
% 1. FreshHeart_08072022.m
% 1. QuantMap_R2_Exvivo.m
% 2. FFmap_R2_Phase_Exvivo.m


clear all;
close all;

% T2-SE
addpath('../function/')
qMRinfo('mono_t2'); % set it up first in qMRLab
dicom_dir = uigetdir;
dicom_glob = glob(cat(2, dicom_dir, '/t2_tse_*'));
dicom_glob_reorder = dicom_glob([2,3,4,5,6,1]); % Carlos, George

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
    f_true = glob(cat(2, f));
    [whatsinit{i} slice_data{i}] = dicom23D(f_true{1}, dicom_fields);
    TE_array(i) = slice_data{i}.EchoTime;
end

%% Draw Mask (rect)
save_dir = GetFullPath(cat(2, dicom_dir, '/Results/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

mask_f = GetFullPath(cat(2, dicom_dir, '/Results/mask_rect.mat'));
mask = zeros(size(whatsinit{1},1), size(whatsinit{1},2), size(whatsinit{1},3));
if ~exist(mask_f)
    for i = 1:size(whatsinit{1},3)
        img = whatsinit{1};
        figure();
        imagesc(img(:,:,i)); axis image;
        roi = drawpolygon;
        mask(:,:,i) = createMask(roi);
    end
    save(mask_f, 'mask');
else
    load(mask_f);
end

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
        ipt = abs(temp_4D(:,:,slc,:));

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
imagesc(t2_map(:,:)); axis image; caxis([0 100]); axis off;
colormap(brewermap([],'RdYlBu')); colorbar;

%% T1-IR-SE
dicom_glob = glob(cat(2, dicom_dir, '/t1_irtse_*'));
dicom_glob_reorder = dicom_glob([1,4,5,6,2,3]); % Carlos
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
    f_true = glob(cat(2, f));
    [whatsinit{i} slice_data{i}] = dicom23D(f_true{1}, dicom_fields);
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
        data.IRData= double(temp_4D(:,:,i,:));
        data.Mask= double(mask(:,:,i));

        Model.Prot.IRData.Mat = IR_array;
        Model.voxelwise = 1;

        % b- fit dataset
        FitResults = FitData(data,Model,0);

        % mask_temp = double(mask(:,:,i));

        % t1_map(:,:,i,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
        t1_map(:,:,i) = FitResults.T1;
    end
end

figure();
imagesc(t1_map(:,:)); axis image; caxis([0 1000]); axis off;
colormap(brewermap([],'RdBu')); colorbar;

%% Save
map_to_save = struct;
map_to_save.t1_map = t1_map;
map_to_save.t2_map = t2_map;

save_dir = GetFullPath(cat(2, save_dir, 'QuantitativeMap.mat'));
save(save_dir, '-struct', 'map_to_save');

