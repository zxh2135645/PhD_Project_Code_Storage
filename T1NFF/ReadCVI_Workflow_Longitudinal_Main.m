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
base_dir = uigetdir; % more generic
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, '/Data/*'));

Names = ExtractNames(folder_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
time_points = {'0D_baseline', '0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR', 'Exvivo'};

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


sequence_label = {'LGE', 'T1', 'T2star'};
sequence_label_exvivo = {'mGRE_3D', 'T1_SIEMENS', 'T2_SIEMENS', 'T2star_SIEMENS', 'T2_CMR', 'T1_CMR', 'T2star_CMR'};
%% Name check
name_check = '11D05';
starting_point = find(strcmp(name_check, Names),1);

%
% Make it not always overwrite
for n = starting_point:length(Names)
%for n = starting_point:starting_point+4
    name = Names{n};
    %for tp = 1:length(time_points)
    for tp = 12:12
        
    
        time_point = time_points{end-tp+1};
        % XML file is independent on Labels
        xml_glob = glob(cat(2, base_dir, '/Data/',  name, '/', name, '_', time_point,  '/*.cvi42wsx'));
        if isempty(xml_glob)
            disp(cat(2, 'Missing CVI XML: ', name, ' ', time_point));
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
                con = con_cell{con_idx};
                if tp == 1
                    for ll = 1:length(sequence_label_exvivo)
                    %for ll = 7:7
                        label = sequence_label_exvivo{ll};
                        % Check if the dstFolder
                        dstFolder = cat(2, OutputPath, name, '/',  name, '_', time_point, '/', label, '/');
                        
                        dicom_glob = glob(cat(2, base_dir, '/Data/', name, '/', name, '_', time_point, '/', label, '/*'));
                        
                        if tp == 13 && (ll == 1)% D3 do not have 2D LGE, and sax5 and sax6 of D3 has different resolution
                            % Needed to be modified here for a more generic use
                            disp(['Skipped ', time_point, '  ', label])
                            
                        else
                            
                            ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_glob, dstFolder, dicom_fields);
                        end
                    end
                else
                    for ll = 1:length(sequence_label)
                        %for ll = 2:2
                            label = sequence_label{ll};
                            % Check if the dstFolder
                            dstFolder = cat(2, OutputPath, name, '/',  name, '_', time_point, '/', label, '/');
                            
                            dicom_glob = glob(cat(2, base_dir, '/Data/', name, '/', name, '_', time_point, '/', label, '/*'));
                            
                        if tp == 13 && (ll == 1)% D3 do not have 2D LGE, and sax5 and sax6 of D3 has different resolution
                            % Needed to be modified here for a more generic use
                            disp(['Skipped ', time_point, '  ', label])
                            
                        else
                            
                            ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_glob, dstFolder, dicom_fields);
                        end
                    end
                end
                
            end
            % break; 
            % Not looking at most chronic time point, this line needs to be
            % deleted for future dev.
        end
    end
end