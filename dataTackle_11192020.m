% This code is for parsing 3 specific LGE examples from Khalid
clear all;
close all;
clc;
current_dir = pwd;
% 3 Lines needs to be modified accordingly

addpath('function/');
addpath('T1NFF/')
base_dir = uigetdir;
popu_label = {'Canine', 'Patient'};
popu = popu_label{2};
name_glob = glob(cat(2, base_dir, '/', popu_label{2}, '/*')); % Modify popu_label here
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'MAG', 'PSIR'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

%name_check = {'RYN'}; % This line needs to be modified accordingly
name_check = {'KIM_BONG_KI'};
name_idx_list = linspace(1, length(Names), length(Names)); % initialize with incremental add

if length(name_check) == 1
    starting_point = find(strcmp(name_check, Names),1);
else
    name_idx_list = zeros(1, length(name_check));
    for n = 1:length(name_check)
        % Check an array of names
        name_idxo = find(strcmp(name_check(n), Names),1);
        name_idx_list(n) = name_idxo;
    end
end

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
    };
output_label = {'LGE', 'T1'};
OutputPath = GetFullPath(cat(2, base_dir, '/ContourData/'));


%% Parsing XML file
% Read DICOM and contours
% This line needs to be modified accordingly

for n = starting_point:length(Names)
    %for n = 1:length(name_idx_list)
    
    name = Names{name_idx_list(n)};
    
    switch name
        case {'RYN', 'SUNNY'}
            old_freeROI_label = 1;
        otherwise
            old_freeROI_label = 0;
    end
    % XML file is independent on Labels
    xml_glob = glob(cat(2, base_dir,'/', popu, '/', name, '/*.cvi42wsx'));
    cvi42wsx = char(xml_glob);
    con_cell = cell(0);
    for xml_ind = 1:size(cvi42wsx, 1)
        % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
        % Thus, there are two different files
        con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
    end
    % Iterate through MAG, PSIR and LGE
    for con_idx = 1:length(con_cell)
        con = con_cell{con_idx};
        for ll = 1:length(sequence_label)
            label = sequence_label{ll};
            
            switch label
                case {'MAG', 'PSIR'}
                    labelo = output_label{1};
                otherwise
                    labelo = output_label{2};
            end
            
            dicom_glob = glob(cat(2, base_dir,'/', popu, '/', name, '/', label, '/*.dcm'));
            
            
            dstFolder = cat(2, OutputPath, popu, '/',  name, '/', labelo, '/');
            
            
            dicom_dir_glob = glob(cat(2, base_dir, '/', popu, '/', name, '/', label, '/'));
            
            ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_dir_glob, dstFolder, dicom_fields, old_freeROI_label);
        end
    end
end


%% after freeROI is drawn config MI to different label
% This line needs to be modified accordingly
for n = starting_point:length(Names)
    %for n = 1:length(name_idx_list)
    name = Names{name_idx_list(n)};
    freeROI_glob = glob(cat(2, OutputPath, popu, '/', name, '/LGE/freeROI/freeROI.mat'));
    load(freeROI_glob{1});
    infarct_masked_3D = freeROIMask_3D;
    
    mi_dir = cat(2, OutputPath, popu, '/', name, '/LGE/LGE_MI/');
    if ~exist(mi_dir, 'dir')
        mkdir(mi_dir);
    end
    
    save(cat(2, mi_dir, 'Infarct_Whole.mat'), 'infarct_masked_3D');
end

