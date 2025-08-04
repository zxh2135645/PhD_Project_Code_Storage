% Read CVI longitudinal main script
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../function/');
addpath('../GUI_n_Analysis_for_LRT/');
addpath('../T1NFF/');

% ================================ Identify your major folder, in this case
% ROC_analysis/
base_dir = uigetdir; % more generic -> % ROC_analysis
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, '/DICOM/*'));

Names = ExtractNames(folder_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
% time_points = {'D6', 'D8'};
% time_points = {'D5', 'D7'};
time_points = {'WK8', 'WK8+2'};
methods = {'LRT', 'ConvCMR'};

OutputPath = GetFullPath(cat(2, base_dir, '/ContourData/'));
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
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

sequence_label = {'LGE', 'EGE', 'EGE1', 'EGE2',  'T2_MULTIECHO'};
% sequence_label = {'temp'};
%% Name check
name_check = 'LATTE_chronic';
starting_point = find(strcmp(name_check, Names),1);

%
% Make it not a lways overwrite
%for n = starting_point:length(Names)
% for n = starting_point:starting_point
for n = starting_point:starting_point
    name = Names{n};
    real_name_temp = strsplit(name, '_');
    real_name = real_name_temp{1};
    for tp = 1:length(time_points)
    %for tp = 1:1

        time_point = time_points{end-tp+1};


        % XML file is independent on Labels
        % Redirect XML file to XML_Data

        xml_glob = glob(cat(2, base_dir, '/XML_Data/',  real_name, '_', time_point,  '*.cvi42wsx'));
        xml_glob2 = glob(cat(2, base_dir, '/XML_Data/',  real_name, '_', time_point, '*.cvi42wsx.xml'));
        if isempty(xml_glob) && isempty(xml_glob2)
            % Try to glob again based on .cvi42wsx.xml
            disp(cat(2, 'Missing CVI XML: ', name, ' ', time_point));
        else
            if ~isempty(xml_glob)
                cvi42wsx = char(xml_glob);
            elseif ~isempty(xml_glob2)
                cvi42wsx = char(xml_glob2);
            end
            con_cell = cell(0);
            for xml_ind = 1:size(cvi42wsx, 1)
                % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
                % Thus, there are two different files
                con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
            end
            % Iterate through MAG, PSIR and LGE
            % for con_idx = 1:length(con_cell)
            for con_idx = 1:1
                % A different label for Exvivo
                con = con_cell{con_idx};
                % for mth = 1:length(methods)
                for mth = 2:2
                    method = methods{mth};
                    % for ll = 1:length(sequence_label)
                    % for ll = 1:1
                    % for ll = 4:4
                    for ll = 5:5
                        label = sequence_label{ll};
                        % Check if the dstFolder
                        dstFolder = cat(2, OutputPath, name, '/', time_point, '/', method, '/', label, '/');

                        dicom_folder = cat(2, base_dir, '/DICOM/', name, '/', time_point, '/', method, '/', label);

                        %ReadCVI_Workflow_ROC_Analysis_Func(con, dicom_folder, dstFolder, dicom_fields);
                        
                        dicom_glob = glob(cat(2, base_dir, '/DICOM/', name, '/', time_point, '/', method, '/', label, '/*'));
                        
                        
                        if (~strcmp(label, 'T2_MULTIECHO') && strcmp(method, 'ConvCMR'))
                            if ~isempty(dicom_glob)
                                ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_glob, dstFolder, dicom_fields);
                            end
                        else
                            if ~isempty(dicom_glob)
                                if strcmp(label, 'T2_MULTIECHO')
                                    % echo_idx_te = 6;
                                    % echo_idx_te = 5; % for CHILI_acute, for DAVE_acute
                                    echo_idx_te = 7; % for Latte_chronic
                                else
                                    echo_idx_te = 1;
                                end
                                old_freeROI_label = 0;
                                ReadCVI_Workflow_ROC_Analysis_Func(con, dicom_folder, dstFolder, dicom_fields, old_freeROI_label, echo_idx_te);
                            end
                        end
                    end
                end

            end

        end
    end
end

% Error in Ryn_0D_baseline T2 resolution is different in base slice