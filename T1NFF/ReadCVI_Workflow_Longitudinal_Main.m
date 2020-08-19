% Read CVI longitudinal main script
clear all;
close all;

cd 'D:\src\Longitudinal'
addpath('D:\src\function');
base_dir = 'D:\T1_Fat_Project\';
folder_glob = glob(cat(2, base_dir, 'Data\*'));

Names = ExtractNames(folder_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
time_points = {'0D', '0D_occl', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};

OutputPath = GetFullPath(cat(2, base_dir, 'ContourData\'));
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

name_check = 'Merry';
starting_point = find(strcmp(name_check, Names),1);

%%
for n = 1:starting_point
    name = Names{n};
    %for tp = 1:length(time_points)
    for tp = 4:4
        time_point = time_points{tp};
        % XML file is independent on Labels
        xml_glob = glob(cat(2, base_dir, 'Data\',  name, '\', name, '_', time_point,  '\*.cvi42wsx'));
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
                
                dicom_glob = glob(cat(2, base_dir, 'Data\', name, '\', name, '_', time_point,'\', label, '\*'));
                
                dstFolder = cat(2, OutputPath, name, '\',  name, '_', time_point, '\', label, '\');
                
                
                if tp == 4 && (ll == 1)% D3 do not have 2D LGE, and sax5 and sax6 of D3 has different resolution
                    disp(['Skipped', time_point, '  ', label])
                else
                    
                    ReadCVI_Workflow_Longitudinal_Study_Func(con, dicom_glob, dstFolder, dicom_fields);
                end
            end
        end
    end
end