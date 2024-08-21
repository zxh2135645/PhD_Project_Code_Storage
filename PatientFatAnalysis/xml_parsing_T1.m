clear all;
close all;
clc;
current_dir = pwd;
% Patient data configuration for Khalid
%% 
addpath('../function/');
base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/Yinyin_Patient_data/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

% sequence_label = {'MAG', 'PSIR', 'T2star', 'T1MOLLI'};
sequence_label = {'T1MOLLI'};
%sequence_label = {'T2star'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

% What's going on with 27?
name_check = {'484060000001'};
%name_check = {'KIM_BONG_KI', 'HAN_BONG_SANG'};
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

% output_label = {'LGE', 'T2star'};
output_label = {'T1'};
%output_label = {'T2star'};
OutputPath = GetFullPath(cat(2, base_dir, '/ContourData/'));

%timepoints = {'BL'};
timepoints = {'FU'};

%% Parsing XML file
% Read DICOM and contours
%for n = 28:28
for n = starting_point:length(Names)
%for n = starting_point:starting_point
    name = Names{name_idx_list(n)};
for tp = 1:length(timepoints)
%for tp = 1:3
    time_point = timepoints{tp};
    
    % XML file is independent on Labels
    xml_glob = glob(cat(2, base_dir,'/Updated_XML/', name, '/', time_point, '/*.cvi42wsx'));
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

        %for ll = 3:3
    % for ll = 1:2
    label = sequence_label{ll};
    labelo = output_label{ll};
    sigOtherLabel = sequence_label{ll};

    dicom_glob = glob(cat(2, base_dir, '/Yinyin_Patient_data/', name, '/', time_point, '/', label, '/*'));
    % MAG and PSIR is mutually exclusive
    % Only contours in T1Map, but also need to export T1-weighted images
    
    sig_dicom_glob = glob(cat(2, base_dir, '/Yinyin_Patient_data/', name, '/', time_point, '/', sigOtherLabel, '/*'));
    sig_dicom = char(sig_dicom_glob);
    
   
    dstFolder = cat(2, OutputPath, name, '/', time_point, '/', labelo);
    
    dicom = char(dicom_glob);
    id_cell = cell(size(dicom, 1), 1);
    total_match = 0;
    % SliceLocation_array = [];
    slc_array = [];

    % For T1
    % echotime_cell = {};
    % sliceloc_array = [];
    
    slc_start = 1;
    slc_end = 1;
    clear  mask_heart_3D mask_myocardium_3D mask_blood_3D excludeMask_3D myoRefMask_3D noReflowMask_3D vol_img_3D sig_vol_img_3D sig_vol_img_4D volume_image freeROIMask_3D
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
        match_count = double(match_count > 0);


        
        
        % get all contours from excludeContour
        excludeMask_3Ds = zeros(size(volume_image));
        if ~isempty(excludeContour)
            keys = fieldnames(excludeContour);
 
            for j = 1:length(keys)
                temp_cell = getfield(excludeContour, keys{j});
                temp_mat = zeros(size(volume_image));
                temp_idx = temp_cell{2};
                if size(temp_cell{1}, 3) > 1
                    if length(temp_idx) == size(temp_cell{1}, 3)
                        temtem_idx = temp_idx;
                        for k = 1:size(temp_cell{1}, 3)
                            temp_mat(:,:,temp_idx(k)) = temp_cell{1}(:,:,k);
                        end
                    else
                        for k = 1:length(temp_idx)
                            innd = find(temp_idx(k) == temtem_idx);
                            temp_mat(:,:,temp_idx(k)) = temp_cell{1}(:,:,innd);
                        end
                    end
                    
                else
                    temp_mat = temp_cell{1};
                end
                excludeMask_3Ds = temp_mat + excludeMask_3Ds;
            end
        end
        
        excludeMask_3Ds(excludeMask_3Ds > 0) = 1; 
        % I realized overlapping exclusion contour between
        % excludeEnhancementArea and excludeEnhancementArea0001
        % Make sure the mask is binary
        
        % Get all contours from NoReFlowArea
        noReflowMask_3Ds = zeros(size(volume_image));
        if ~isempty(noReflowCell)
            if ~isempty(noReflowCell{1})
                temp_mat = zeros(size(volume_image));
                temp_idx = noReflowCell{2};
                %if size(noReflowCell{1}, 3) > 1
                for j = 1:size(noReflowCell{1}, 3)
                    temp_mat(:,:,temp_idx(j)) = noReflowCell{1}(:,:,j);
                end
                %else
                %    temp_mat = noReflowCell{1};
                %end
                noReflowMask_3Ds = temp_mat + noReflowMask_3Ds;
            end
        end
        
        myoRefMask_3Ds = zeros(size(volume_image));
        if ~isempty(myoRefCell)
            if ~isempty(myoRefCell{1})
                temp_mat = zeros(size(volume_image));
                temp_idx = myoRefCell{2};
                %if size(myoRefCell{1}, 3) > 1
                for j = 1:size(myoRefCell{1}, 3)
                    temp_mat(:,:,temp_idx(j)) = myoRefCell{1}(:,:,j);
                end
                %else
                %    temp_mat = myoRefCell{1};
                %end
                myoRefMask_3Ds = temp_mat + myoRefMask_3Ds;
            end
        end
        
        freeROIMask_3Ds = zeros(size(volume_image));
        if ~isempty(freeROICell)
            if ~isempty(freeROICell{1})
                temp_mat = zeros(size(volume_image));
                temp_idx = freeROICell{2};
                %if size(freeROICell{1}, 3) > 1
                for j = 1:size(freeROICell{1}, 3)
                    temp_mat(:,:,temp_idx(j)) = freeROICell{1}(:,:,j);
                    
                end
                %else
                %    temp_mat = freeROICell{1};
                %end
                freeROIMask_3Ds = temp_mat + freeROIMask_3Ds;
            end
        end

        if match_count > 0

            % T2star will be the first echo
            vol_img_3D(:,:,total_match + match_count) = volume_image(:,:,1);
            mask_heart_3D(:,:,total_match + match_count) = mask_heart(:,:,1);
            mask_myocardium_3D(:,:,total_match + match_count) = mask_myocardium(:,:,1);
            mask_blood_3D(:,:,total_match + match_count) = mask_blood(:,:,1);
            excludeMask_3D(:,:,total_match + match_count) = excludeMask_3Ds(:,:,1);
            myoRefMask_3D(:,:,total_match + match_count)  = myoRefMask_3Ds(:,:,1);
            noReflowMask_3D(:,:,total_match + match_count)  = noReflowMask_3Ds(:,:,1);
            freeROIMask_3D(:,:,total_match + match_count)  = freeROIMask_3Ds(:,:,1);

            % total_match = 0;
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

                        % This needs to be modified
                        vida_invt = sig_slice_info.PerFrameFunctionalGroupsSequence.Item_1.Private_0021_11fe.Item_1.Private_0021_1189;
                        sig_slice_data = setfield(sig_slice_data, {m}, 'InversionTime', vida_invt);
                    end
                else
                    [sig_volume_image, sig_slice_data, sig_image_meta_data] = dicom23D(sig_dicom(s,:), dicom_fields);
                end

            end

            total_match = total_match + match_count;


            if strcmp(label, 'T2star')
                slc_array = [slc_array, slice_data(1).SliceLocation];
                glob_idx = [glob_idx, i];
            else
                slc_array = [slc_array, slice_data.SliceLocation];
            end
        end
    end



        
    
    %% Need to match MAG with PSIR or vice versa
    
    %% Save all as mat file
    %  dsts = {'Heart', 'Myocardium', 'excludeContour', 'myoReference', 'noReflowAreaContour', 'BloodPool'};
    dsts = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'BloodPool', 'freeROI'};

    if total_match ~= 0


        T1_struct = struct;
        T1_struct.vol_img_3D = vol_img_3D;
        dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
        save(dstPath, 'T1_struct');


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

        dstPath = cat(2, dstFolder, '/', dsts{7});
        if ~ exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        save(cat(2, dstPath, '/freeROI.mat'), 'freeROIMask_3D');

        dstPath = cat(2, dstFolder, '/', labelo, '_SliceLoc.mat');
        save(dstPath, 'slc_array');

        if strcmp(label, 'T2star')
            glob_names = cell(1, length(glob_idx));
            for i = 1:length(glob_names)
                dicom = dicom_glob{glob_idx(i)};
                strings = strsplit(dicom, '/');
                glob_names{i} = strings{end-1};
            end
            dstPath = cat(2, dstFolder, '/', label, '_Index.mat');
            save(dstPath, 'glob_names');
        end

        disp(name)
        disp('Done!')
    end
    end
    end
end
end
   