clear all;
close all;
clc;
current_dir = pwd;
% Patient data configuration for Khalid
%% 
addpath('function\');
base_dir = uigetdir;
% base_dir needs to be changed

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'MAG', 'PSIR', 'T1Map', 'T1w_MOCO'};


name_check = 'CHOI_TAE_SIK';
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

    output_label = {'LGE', 'T1'};

%% Parsing XML file
% Read DICOM and contours
%for n = starting_point:starting_point
for n = starting_point:length(Names)
    name = Names{n};
    % XML file is independent on Labels
    xml_glob = glob(cat(2, base_dir, '/', name, '/*.cvi42wsx'));
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
     for ll = 1:(length(sequence_label)-1)
     %for ll = 3:3
            label = sequence_label{ll};
            
            switch label
                case {'MAG', 'PSIR'}
                    labelo = output_label{1};
                otherwise
                    labelo = output_label{2};
            end
    
    dicom_glob = glob(cat(2, base_dir, '/', name, '/', label, '/*'));
    % MAG and PSIR is mutually exclusive
    % Only contours in T1Map, but also need to export T1-weighted images
    if strcmp(label,sequence_label{1})  % equals 1 or 2
        sigOtherIdx = 2;
    elseif strcmp(label, sequence_label{2})
        sigOtherIdx = 1;
    else
        sigOtherIdx = 4;
    end
        
    sigOtherLabel = sequence_label{sigOtherIdx};
    sig_dicom_glob = glob(cat(2, base_dir, '/', name, '/', sigOtherLabel, '/*'));
    sig_dicom = char(sig_dicom_glob);
    
    
    OutputPath = GetFullPath(cat(2, base_dir, '/../ContourData/'));
    dstFolder = cat(2, OutputPath, name, '\', labelo);
    
    dicom = char(dicom_glob);
    id_cell = cell(size(dicom, 1), 1);
    total_match = 0;
    
    slc_start = 1;
    slc_end = 1;
    clear vol_img_3D mask_heart_3D mask_myocardium_3D mask_blood_3D excludeMask_3D myoRefMask_3D noReflowMask_3D sig_vol_img_3D sig_vol_img_4D volume_image
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
        excludeMask_2D = zeros(size(volume_image));
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
                excludeMask_2D = temp_mat + excludeMask_2D;
            end
        end
        
        % Get all contours from NoReFlowArea
        noReflowMask_2D = zeros(size(volume_image));
        if ~isempty(noReflowCell)
            temp_mat = noReflowCell{1};
            for j = 1:size(noReflowCell{1}, 3)
                noReflowMask_2D = temp_mat(:,:,j) + noReflowMask_2D;
            end
        end
        
        % Get all contours from myoRefCell
        % Why there is a empty RefMat
        myoRefMask_2D = zeros(size(volume_image));
        if ~isempty(myoRefCell)
            temp_mat = myoRefCell{1};
            for j = 1:size(myoRefCell{1}, 3)
                myoRefMask_2D = temp_mat(:,:,j) + myoRefMask_2D;
            end
        end
        
        if match_count > 0
            % All files of LGE and T1 are under T1 folder
            vol_img_3D = volume_image;
            mask_heart_3D = mask_heart;
            mask_myocardium_3D = mask_myocardium;
            mask_blood_3D = mask_blood;
            excludeMask_3D = excludeMask_2D;
            myoRefMask_3D  = myoRefMask_2D;
            noReflowMask_3D  = noReflowMask_2D;
            
            
            % Find sigOther in matched slices
            % MAG and PSIR is mutually exclusive
            if strcmp(label, sequence_label{1}) || strcmp(label, sequence_label{2}) % equals 1 or 2
                for s = 1:size(sig_dicom,1)
                    sig_dicom_file = glob(cat(2, sig_dicom(s,:), '*.dcm'));
                    sig_dicom_f = sig_dicom_file{1};
                    sig_info = dicominfo(sig_dicom_f);
                    
                    if contains(sig_info.InstitutionName, 'Vida') && ~isfield(sig_info, 'SliceLocation')
                        [sig_volume_image, sig_slice_data] = convert_vida_header(sig_dicom_f, sig_info);
                    else
                        [sig_volume_image, sig_slice_data, sig_image_meta_data] = dicom23D(sig_dicom(s,:), dicom_fields);
                    end

                    sig_vol_img_3D = sig_volume_image;
                    
                end
            else
                % T1w_MOCO will be 3D matrix
                invt_cell = {};
                for s = 1:size(sig_dicom,1)
                    sig_dicom_file = glob(cat(2, sig_dicom(s,:), '*.dcm'));
                    sig_dicom_f = sig_dicom_file{1};
                    sig_info = dicominfo(sig_dicom_f);
                    
                    if contains(sig_info.InstitutionName, 'Vida') && ~isfield(sig_info, 'SliceLocation')
                        sig_volume_image = zeros(size(volume_image, 1), size(volume_image, 2), length(sig_dicom_file));
                        
                        for m = 1:length(sig_dicom_file)
                            
                            [sig_slice_image, sig_slice_info] = convert_vida_header(sig_dicom_file{m}, dicominfo(sig_dicom_file{m}));
                            sig_volume_image(:,:,m) = sig_slice_image;
                            sig_slice_data = setfield(sig_slice_data, {m}, 'SliceLocation', sig_slice_info.SliceLocation);
                            
                            vida_invt = sig_slice_info.PerFrameFunctionalGroupsSequence.Item_1.Private_0021_11fe.Item_1.Private_0021_1189;
                            sig_slice_data = setfield(sig_slice_data, {m}, 'InversionTime', vida_invt);
                        end
                    else
                        [sig_volume_image, sig_slice_data, sig_image_meta_data] = dicom23D(sig_dicom(s,:), t1w_dicom_fields);
                    end
                    
                    for ss = 1:size(sig_dicom, 1)
                        if slice_data(ss).SliceLocation == sig_slice_data(1).SliceLocation
                            sig_vol_img_4D(:,:,:,slc_start:slc_end+match_count-1) = sig_volume_image;
                            invt_array = zeros(size(sig_vol_img_4D, 3), 1);
                            for inv = 1:size(sig_vol_img_4D, 3)
                                invt_array(inv) = sig_slice_data(inv).InversionTime;
                            end
                            invt_cell{end+1} = invt_array;
                            slc_start = slc_start + match_count;
                            slc_end = slc_end + match_count;
                            total_match = total_match + match_count;
                        end
                    end
                end
            end

        end
        
    end
    
        
    
    %% Need to match MAG with PSIR or vice versa
    
    %% Save all as mat file
    %  dsts = {'Heart', 'Myocardium', 'excludeContour', 'myoReference', 'noReflowAreaContour', 'BloodPool'};
    dsts = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'BloodPool'};
    
    % TODO
    if total_match ~= 0
        dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
        save(dstPath, 'vol_img_3D');
        
        if strcmp(label, sequence_label{1}) || strcmp(label, sequence_label{2})
            dstPath = cat(2, dstFolder, '/', sigOtherLabel, '_vol_img_3D.mat');
            save(dstPath, 'sig_vol_img_3D');
        else
            sig_vol_img_4D = permute(sig_vol_img_4D, [1,2,4,3]);
            T1w_MOCO_struct = struct;
            T1w_MOCO_struct.vol_img_4D = sig_vol_img_4D;
            T1w_MOCO_struct.invt_cell = invt_cell;
            dstPath = cat(2, dstFolder, '/', sigOtherLabel, '_vol_img_4D.mat');
            save(dstPath, 'sig_vol_img_4D');
        
            % Save T1w images and inversion time
        end
        
        dstPath = cat(2, dstFolder, '/', dsts{1});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/mask_heart.mat'), 'mask_heart_3D');
        
        dstPath = cat(2, dstFolder, '/', dsts{2});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2,dstPath, '/mask_myocardium.mat'), 'mask_myocardium_3D');
        
        dstPath = cat(2, dstFolder, '/', dsts{6});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/mask_blood.mat'), 'mask_blood_3D');
        
        dstPath = cat(2, dstFolder, '/', dsts{3});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/excludeArea.mat'), 'excludeMask_3D');
        
        dstPath = cat(2, dstFolder, '/', dsts{4});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/myoRef.mat'), 'myoRefMask_3D');
        
        dstPath = cat(2, dstFolder, '/', dsts{5});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/noReflow.mat'), 'noReflowMask_3D');
        
        disp(name)
        disp('Done!')
    end
    end
    end
end
   
%% Check the files
data_glob = glob(cat(2, base_dir, '/*'));
contour_glob = glob(cat(2, OutputPath, '*'));

Names_data = cell(length(data_glob), 1);
for i = 1:length(data_glob)
    strings = strsplit(data_glob{i},'\');
    Names_data{i} = strings{end-1};
end

Names_contour = cell(length(contour_glob), 1);
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'\');
    Names_contour{i} = strings{end-1};
end

disp(length(data_glob));
disp(length(contour_glob));
contour_idx = [];
data_idx = linspace(1, length(Names_data), length(Names_data));
for j = 1:length(Names_contour)
    name_contour = Names_contour{j};
    idx =  find(strcmp(Names_data, name_contour));
    contour_idx(end+1) = idx;
end

diff = setdiff(data_idx, contour_idx);
disp(Names_data(diff))






