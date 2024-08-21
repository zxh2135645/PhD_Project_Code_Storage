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

% ================================ Identify your major folder, in this case
% Example/
base_dir = uigetdir; % more generic
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
% folder_glob = glob(cat(2, base_dir, '/T1FP_Data_02242022/Data/*'));
folder_glob = glob(cat(2, base_dir, '/HV_Data/*'));

Names = ExtractNames(folder_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
%time_points = {'0D_baseline', '0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR', 'Exvivo'};
% time_points = {'6MO', '9MO', '1YR', '15YR', 'Exvivo'};

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


sequence_label = {'T1MOLLI', 'T2star', 'T2'};
sequence_label_exvivo = {'mGRE_3D', 'T1_SIEMENS', 'T2_SIEMENS', 'T2star_SIEMENS', 'T2_CMR', 'T1_CMR', 'T2star_CMR'};
%% Name check
name_check = '2795';
starting_point = find(strcmp(name_check, Names),1);

%
% Make it not always overwrite
for n = starting_point:length(Names)
%for n = starting_point:starting_point
    name = Names{n};
    for tp = 1:length(time_points)
    %for tp = 7:7
    
        time_point = time_points{end-tp+1};
        % XML file is independent on Labels
        % Redirect XML file to XML_Data
        % xml_glob = glob(cat(2, base_dir, '/T1FP_Data_02242022/Data/',  name, '/', name, '_', time_point,  '/*.cvi42wsx'));
        % xml_glob2 = glob(cat(2, base_dir, '/T1FP_Data_02242022/Data/',  name, '/', name, '_', time_point,  '/*.cvi42wsx.xml'));
        xml_glob = glob(cat(2, base_dir, '/Updated_XML/',  name,   '/*.cvi42wsx'));
        xml_glob2 = glob(cat(2, base_dir, '/Updated_XML/',  name,  '/*.cvi42wsx.xml'));
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
            for con_idx = 1:length(con_cell)
                % A different label for Exvivo
                con = con_cell{con_idx};

                for ll = 1:length(sequence_label)
                %for ll = 1:1
                    label = sequence_label{ll};
                    % Check if the dstFolder
                    dstFolder = cat(2, OutputPath, name, '/', label, '/');

                    dicom_glob = glob(cat(2, base_dir, '/HV_Data/', name,  '/DICOM/', label, '/*'));

                    % ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_glob, dstFolder, dicom_fields);
                    ReadCVI_Workflow_Longitudinal_Study_Cell_Func(con, dicom_glob, dstFolder, dicom_fields);


                end
            end
        end
    end
end

% Error in Ryn_0D_baseline T2 resolution is different in base slice