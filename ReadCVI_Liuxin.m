% Read CVI longitudinal main script
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./function/');
addpath('./GUI_n_Analysis_for_LRT/');
addpath('./T1NFF/');
addpath('./PatientFatAnalysis/');

% ================================ Identify your major folder, in this case
% ROC_analysis/
base_dir = uigetdir; % more generic -> % ROC_analysis
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, '/DICOM/*'));

Names = ExtractNames(folder_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
% time_points = {'D6', 'D8'};
% time_points = {'D5', 'D7'};
time_points = {''};

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

sequence_label = {'BOOST', 'LGE'};
sequence_label = {'BOOST2', 'LGE', 'T2'};
sequence_label = {'BOOST2'};
sequence_label = {'LGE', 'T2'};
% sequence_label = {'temp'};
%% Name check
name_check = 'BAI_YAN_SHENG';
%name_check = 'LIANG_JI_CHUN';
%name_check = 'WANG_SHU_LI';
%name_check = 'CHENG_YONG';
%name_check = 'ZHENG_YI';
name_check = 'SHU_WEN_LI';
name_check =   'XU_YONG';
name_check = 'ZHENG_YI';
starting_point = find(strcmp(name_check, Names),1);

% Issues in CHENG_YONG, CHEN_JING, CHEN_PENG, FENG_LEI, FENG_SHU_FANG,
% JIANG_XIN_HUA, LIANG_JI_CHUN, 'LIU_CHEN_HAO', 'LI_JIN', 'LI_JING_HUI',
% 'LU_HAI_HUI', 

% Missing LGE: JING_AN_KUN

% Something wrong with 'PEI_LI_YAN' - LGE - excludeContour
% 'XU_JI_JIE', - LGE - excludeContour
% 'ZHAO_ZHONG_FU' - T2 - excludeContour

%
% Make it not always overwrite
%for n = starting_point:length(Names)
for n = starting_point:length(Names)
%for n = starting_point:starting_point
    name = Names{n};
    real_name_temp = strsplit(name, '_');
    real_name = real_name_temp{1};
    for tp = 1:length(time_points)
        %for tp = 1:1

        time_point = time_points{end-tp+1};


        % XML file is independent on Labels
        % Redirect XML file to XML_Data

        xml_glob = glob(cat(2, base_dir, '/XML_Data/',  name, '/', time_point,  '/*.cvi42wsx'));
        xml_glob2 = glob(cat(2, base_dir, '/XML_Data/',  name, '/', time_point, '/*.cvi42wsx.xml'));
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
                xml_temp = strrep(cvi42wsx(xml_ind,:), ' ', '');
                con_cell{end+1} = CMR42ContourReader(xml_temp);
            end
            % Iterate through MAG, PSIR and LGE
            for con_idx = 1:length(con_cell)
                % for con_idx = 1:1
                % A different label for Exvivo
                con = con_cell{con_idx};

                for ll = 1:length(sequence_label)
                    label = sequence_label{ll};
                    % Check if the dstFolder
                    dstFolder = cat(2, OutputPath, name, '/', time_point, '/', label, '/');

                    dicom_folder = cat(2, base_dir, '/DICOM/', name, '/', time_point, '/', label);

                    %ReadCVI_Workflow_ROC_Analysis_Func(con, dicom_folder, dstFolder, dicom_fields);

                    dicom_glob = glob(cat(2, base_dir, '/DICOM/', name, '/', time_point, '/', label, '/*'));

                    old_freeROI_label = 0;
                    echo_idx_te = 1;
                    ReadCVI_Workflow_Liuxin_Analysis_Func(con, dicom_folder, dstFolder, dicom_fields, old_freeROI_label, echo_idx_te);
                    
                    % if ~isempty(dicom_glob)
                    %     for dd = 1:length(dicom_glob)
                    %         dicom = dicom_glob{dd};
                    % 
                    %         old_freeROI_label = 0;
                    %         ReadCVI_Workflow_IndianPatient_Analysis_Func(con, dicom, dstFolder, dicom_fields, old_freeROI_label, echo_idx_te);
                    %     end
                    % end

                end

            end

        end
    end
end

        % Error in Ryn_0D_baseline T2 resolution is different in base slice