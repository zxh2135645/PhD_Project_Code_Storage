clear all;
close all;

clc;
current_dir = pwd;

%% 
addpath('..\function\');
base_dir = GetFullPath(cat(2, current_dir, '\..\..\Data\Diane\ResultsFolder_180718\'));
name_glob = glob(cat(2, base_dir, '*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'MultiEcho'};


name_check = 'Queenie_Final_28May2018';
starting_point = find(strcmp(name_check, Names),1);

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

output_label = {'LGE', 'T1', 'MultiEcho'};

%% Parsing XML file
% Read DICOM and contours
% for n = starting_point:starting_point
for n = starting_point:starting_point
    name = Names{n};
    % XML file is independent on Labels
    xml_glob = glob(cat(2, base_dir, name, '/*.cvi42wsx'));
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
                 
             case {'T1Map'}
                 labelo = output_label{2};
                 
             otherwise
                 labelo = output_label{3};
         end
         
    dicom_glob = glob(cat(2, base_dir, name, '/', label, '/*'));    
    
    
    OutputPath = GetFullPath(cat(2, base_dir, '../ContourData/'));
    dstFolder = cat(2, OutputPath, name, '\', labelo);
    
    dicom = char(dicom_glob);
    id_cell = cell(size(dicom, 1), 1);
    total_match = 0;
    
    slc_start = 1;
    slc_end = 1;
    clear vol_img_4D mask_heart_4D mask_myocardium_4D mask_blood_4D excludeMask_4D myoRefMask_4D noReflowMask_4D freeROIMask_4D volume_image
    for i = 1:size(dicom, 1)
        dicom_file = glob(cat(2, dicom(i,:), '*.dcm'));        
        dicom_f = dicom_file{1};
        info = dicominfo(dicom_f);
        
        if contains(info.InstitutionName, 'Vida') && ~isfield(info, 'SliceLocation')
            
            [volume_image, slice_data] = convert_vida_header(dicom_f, info);
        else
            [volume_image, slice_data, image_meta_data] = dicom23D(dicom(i,:), dicom_fields);
        end
        id_cell{i} = slice_data.MediaStorageSOPInstanceUID; % Didn't do anything to it
        [mask_heart, mask_myocardium, mask_blood, excludeContour, myoRefCell, noReflowCell, freeROICell, match_count] = ...
            CMR42ContourMatrixGenerator3(con, volume_image, slice_data, dstFolder);
        
        % get all contours from excludeContour
        excludeMask_3D = zeros(size(volume_image));
        if ~isempty(excludeContour)
            keys = fieldnames(excludeContour);
            % The code below assumes 2D slice of image; can be improved for more
            % generic use.
            for j = 1:length(keys)
                temp_cell = getfield(excludeContour, keys{j});
                temp_mat = zeros(size(volume_image));
                if size(temp_cell{1}, 3) > 1
                    for k = 1:size(temp_cell{1}, 3)
                        temp_mat = temp_mat + temp_cell{1}(:,:,k);
                    end
                end
                excludeMask_3D = temp_mat + excludeMask_3D;
            end
        end
        
        % Get all contours from NoReFlowArea
        noReflowMask_3D = zeros(size(volume_image));
        if ~isempty(noReflowCell)
            temp_mat = noReflowCell{1};
            for j = 1:size(noReflowCell{1}, 3)
                noReflowMask_3D = temp_mat(:,:,j) + noReflowMask_3D;
            end
        end
        
        % Get all contours from myoRefCell
        % Why there is a empty RefMat
        myoRefMask_3D = zeros(size(volume_image));
        if ~isempty(myoRefCell)
            temp_mat = myoRefCell{1};
            for j = 1:size(myoRefCell{1}, 3)
                myoRefMask_3D = temp_mat(:,:,j) + myoRefMask_3D;
            end
        end
        
        freeROIMask_3D = zeros(size(volume_image));
        if ~isempty(freeROICell)
            temp_mat = freeROICell{1};
            for j = 1:size(freeROICell{1}, 3)
                freeROIMask_3D = temp_mat(:,:,j) + freeROIMask_3D;
            end
        end
        
        if match_count > 0
            vol_img_4D(:,:,:,slc_start:slc_end+match_count-1) = volume_image;
            
            %% TODO
            %te_array = zeros(size(vol_img_4D, 3), 1);
            %for inv = 1:size(vol_img_4D, 3)
            %    te_array(inv) = slice_data(inv).InversionTime;
            %end
            %invt_cell{end+1} = invt_array;
            %%
            mask_heart_4D(:,:,:,slc_start:slc_end+match_count-1) = mask_heart;
            mask_myocardium_4D(:,:,:,slc_start:slc_end+match_count-1) = mask_myocardium;
            mask_blood_4D(:,:,:,slc_start:slc_end+match_count-1) = mask_blood;
            excludeMask_4D(:,:,:,slc_start:slc_end+match_count-1) = excludeMask_3D;
            myoRefMask_4D(:,:,:,slc_start:slc_end+match_count-1)  = myoRefMask_3D;
            noReflowMask_4D(:,:,:,slc_start:slc_end+match_count-1)  = noReflowMask_3D;
            freeROIMask_4D(:,:,:,slc_start:slc_end+match_count-1)  = freeROIMask_3D;

            slc_start = slc_start + match_count;
            slc_end = slc_end + match_count;
            total_match = total_match + match_count;
        end
        
    end
    
        
    
    %% Need to match MAG with PSIR or vice versa
    
    %% Save all as mat file
    %  dsts = {'Heart', 'Myocardium', 'excludeContour', 'myoReference', 'noReflowAreaContour', 'BloodPool'};
    dsts = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'BloodPool', 'freeROI'};
    
    % TODO
    if total_match ~= 0
        dstPath = cat(2, dstFolder, '/', label, '_vol_img_4D.mat');
        save(dstPath, 'vol_img_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{1});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/mask_heart.mat'), 'mask_heart_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{2});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2,dstPath, '/mask_myocardium.mat'), 'mask_myocardium_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{6});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/mask_blood.mat'), 'mask_blood_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{3});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/excludeArea.mat'), 'excludeMask_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{4});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/myoRef.mat'), 'myoRefMask_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{5});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/noReflow.mat'), 'noReflowMask_4D');
        
        dstPath = cat(2, dstFolder, '/', dsts{7});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/freeROI.mat'), 'freeROIMask_4D');
        
        disp(name)
        disp('Done!')
    end
    end
    end
end