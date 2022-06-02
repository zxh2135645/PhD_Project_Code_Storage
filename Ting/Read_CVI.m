clear all;
close all;
clc;
current_dir = pwd;
% Dog data configuration for Ting 06/18/2021
%
addpath('../function/');
addpath('../T1NFF/');
base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'T1Map'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

% name_check = {'AXEL_Day0'};
name_check = {'Evelyn_Day0'};
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

t1w_dicom_fields = {...
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

output_label = {'T1'};

%%
% Make it not always overwrite
for n = starting_point:length(Names)
    %for n = starting_point:starting_point+4
    name = Names{n};
    %for tp = 1:length(time_points)
    
    
    %time_point = time_points{end-tp+1};
    % XML file is independent on Labels
    xml_glob = glob(cat(2, base_dir, '/',  name, '/Data', '/*.cvi42wsx'));
    if isempty(xml_glob)
        disp(cat(2, 'Missing CVI XML: ', name));
    else
        cvi42wsx = char(xml_glob);
        con_cell = cell(0);
        for xml_ind = 1:size(cvi42wsx, 1)
            % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
            % Thus, there are two different files
            con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
        end
        % Iterate through MAG, PSIR and LGE
        for con_idx = 1:length(con_cell)
            % A different label for Exvivo
            con = con_cell{con_idx};
            
            for ll = 1:length(sequence_label)
                label = sequence_label{ll};
                label_out = output_label{ll};
                % Check if the dstFolder
                dstFolder = cat(2, base_dir, '/', name, '/ContourData/', label_out, '/');
                
                dicom_glob = glob(cat(2, base_dir, '/', name, '/Data/', label, '/*'));
                
                ReadCVI_Workflow_Cell_Func(con, dicom_glob, dstFolder, dicom_fields);
            end
            
            
        end
    end
    %end
end