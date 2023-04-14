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
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

sequence_label = {'PSIR', 'MAG'};
% sequence_label = {'DICOM_E', 'DICOM_L', 'PSIR', 'MAG'};
% sequence_label = {'LRT_LGE'};
time_label = {'D6', 'D8', 'WK8', 'WK8+2'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
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
output_label = {'PSIR', 'MAG'};
% output_label = {'LRT_LGE'};
%% Main Body
%for n = 1:length(Names)
for n = 4:4
    name = Names{n};
    for t = 1:length(time_label)
    %for t = 1:1
        cvi_glob = glob(cat(2, name_glob{n}, name, '_', time_label{t}, '/', name, '*.cvi42wsx.xml'));
        if isempty(cvi_glob)
            disp(cat(2, 'NO CVI XML: ', name));
        else
            cvi42wsx = char(cvi_glob);
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
                    dstFolder = cat(2, base_dir, '/', name, '/', name, '_', time_label{t}, '/ContourData/', label_out, '/');

                    if strcmp(label, 'DICOM_E') || strcmp(label, 'DICOM_L')
                        dicom_glob = glob(cat(2, base_dir, '/', name, '/Data/LRT_recon_Seg15/', label, '/'));
                    else
                        % dicom_glob = glob(cat(2, base_dir, '/', name, '/Data/OrigData/*', label, '_[0123456789]*/'));
                        dicom_glob = glob(cat(2, name_glob{n}, name, '_', time_label{t}, '/*/*', label, '_[0123456789]*/'));
                    end

                    ReadCVI_Workflow_Cell_Func(con, dicom_glob, dstFolder, dicom_fields);
                end

            end
        end
    end
end

