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

sequence_label = {'MAG', 'PSIR', 'T2star'};
%sequence_label = {'T1MOLLI'};
%sequence_label = {'T2star'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

% What's going on with 27?
name_check = {'484060000018'};
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

t2star_dicom_fields = {...
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
    'EchoTime',...
    };

output_label = {'LGE', 'T2star'};
%output_label = {'T2star'};
%output_label = {'T1'};
OutputPath = GetFullPath(cat(2, base_dir, '/ContourData/'));

timepoints = {'BL'};
% timepoints = {'FU'};
%% Parsing XML file
% Read DICOM and contours
%for n = 28:28
%for n = starting_point:length(Names)
for n = starting_point:starting_point
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
        
    %for ll = 1:length(sequence_label)
    for ll = 3:3
    %for ll = 1:2
            label = sequence_label{ll};
            
            switch label
                case {'MAG', 'PSIR'}
                    labelo = output_label{1};
                otherwise
                    labelo = output_label{2};
            end
    
    dicom_glob = glob(cat(2, base_dir,'/Yinyin_Patient_data/', name, '/', time_point, '/', label, '/*'));
    % MAG and PSIR is mutually exclusive
    % Only contours in T1Map, but also need to export T1-weighted images
    if strcmp(label,sequence_label{1})  % equals 1 or 2
        sigOtherIdx = 2;
    elseif strcmp(label, sequence_label{2})
        sigOtherIdx = 1;
    else
        sigOtherIdx = 3;
        glob_idx = [];
    end
        
    sigOtherLabel = sequence_label{sigOtherIdx};
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

        if strcmp(label, 'T2star')
            slice_data = slice_data(1);
            volume_image = volume_image(:,:,1);
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
                        
            if strcmp(label, sequence_label{1}) || strcmp(label, sequence_label{2}) % equals 1 or 2
                vol_img_3D = volume_image;
                mask_heart_3D = mask_heart;
                mask_myocardium_3D = mask_myocardium;
                mask_blood_3D = mask_blood;
                excludeMask_3D = excludeMask_3Ds;
                myoRefMask_3D  = myoRefMask_3Ds;
                noReflowMask_3D  = noReflowMask_3Ds;
                freeROIMask_3D  = freeROIMask_3Ds;
                % Find sigOther in matched slices
                % MAG and PSIR is mutually exclusive
            
            % Find sigOther in matched slices
            % MAG and PSIR is mutually exclusive
            % All files of LGE and T1 are under T1 folder
           
                % Need another array to encode info of which slice is
                % accpected if there is multiple version of same slice
                % i is the index of dicom (has the same order as its counterpart)
                
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
                 % T2star will be the first echo
                 vol_img_3D(:,:,total_match + match_count) = volume_image(:,:,1);
                 mask_heart_3D(:,:,total_match + match_count) = mask_heart(:,:,1);
                 mask_myocardium_3D(:,:,total_match + match_count) = mask_myocardium(:,:,1);
                 mask_blood_3D(:,:,total_match + match_count) = mask_blood(:,:,1);
                 excludeMask_3D(:,:,total_match + match_count) = excludeMask_3Ds(:,:,1);
                 myoRefMask_3D(:,:,total_match + match_count)  = myoRefMask_3Ds(:,:,1);
                 noReflowMask_3D(:,:,total_match + match_count)  = noReflowMask_3Ds(:,:,1);
                 freeROIMask_3D(:,:,total_match + match_count)  = freeROIMask_3Ds(:,:,1);
                 
                echotime_cell = {};
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
                        [sig_volume_image, sig_slice_data, sig_image_meta_data] = dicom23D(sig_dicom(s,:), t2star_dicom_fields);
                    end
                
                
                    % for ss = 1:size(sig_dicom, 1)
                    if slice_data(1).SliceLocation == sig_slice_data(1).SliceLocation
                        sig_vol_img_4D(:,:,:,slc_start:slc_end+match_count-1) = sig_volume_image;
                        echotime_array = zeros(size(sig_vol_img_4D, 3), 1);
                        for inv = 1:size(sig_vol_img_4D, 3)
                            echotime_array(inv) = sig_slice_data(inv).EchoTime;
                        end
                        echotime_cell{end+1} = echotime_array;
                        slc_start = slc_start + match_count;
                        slc_end = slc_end + match_count;
                        % total_match = total_match + match_count;
                    end
                    % end
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
        
        if strcmp(label, sequence_label{1}) || strcmp(label, sequence_label{2})
            vi3 = struct;
            vi3.vol_img_3D = vol_img_3D;
            dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
            save(dstPath, 'vi3');
            
            svi3 = struct;
            svi3.sig_vol_img_3D = sig_vol_img_3D;
            
            dstPath = cat(2, dstFolder, '/', sigOtherLabel, '_vol_img_3D.mat');
            save(dstPath, 'svi3');
        else
            vi3 = struct;
            vi3.vol_img_3D = vol_img_3D;
            dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
            save(dstPath, 'vi3');
            
            sig_vol_img_4D = permute(sig_vol_img_4D, [1,2,4,3]);
            T2star_struct = struct;
            T2star_struct.vol_img_4D = sig_vol_img_4D;
            T2star_struct.echotime_cell = echotime_cell;
            dstPath = cat(2, dstFolder, '/', sigOtherLabel, '_vol_img_4D.mat');
            save(dstPath, 'T2star_struct');
        
            % Save T2star images and echo time
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
   
%% Check the files
data_glob = glob(cat(2, base_dir, '/Updated_XML/*'));
contour_glob = glob(cat(2, OutputPath, '*'));

Names_data = cell(length(data_glob), 1);
for i = 1:length(data_glob)
    strings = strsplit(data_glob{i},'/');
    Names_data{i} = strings{end-1};
end

Names_contour = cell(length(contour_glob), 1);
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'/');
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

%% Check img dimension

contour_glob = glob(cat(2, OutputPath, '*'));
Names_contour = cell(length(contour_glob), 1);
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'/');
    Names_contour{i} = strings{end-1};
end

for i = 1:length(Names_contour)
    name = Names_contour{i};
       for ll = 1:(length(sequence_label)-1)
            clear vol_img_3D sig_vol_img_3D vi3 svi3
            label = sequence_label{ll};
            
            switch label
                case {'MAG', 'PSIR'}
                    labelo = output_label{1};
                otherwise
                    labelo = output_label{2};
            end
            
            load(cat(2, OutputPath, name, '/', time_point,  '/', labelo, '/', label, '_vol_img_3D.mat'));
            
            if exist('vi3', 'var')
                img = vi3.vol_img_3D;
            else
                img = svi3.sig_vol_img_3D;
            end


            disp(name)
            disp(label)
            disp(size(img))
            
       end
end










