clear all;
close all;
% So this goes for AHA segmented analysis for invivo data, because simply
% comparing volumes is biased (As its name indicated)

addpath('../function/');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

dicom_dir = uigetdir;
dicom_glob = glob(cat(2, dicom_dir, '/*'));

name_glob = glob(cat(2, dicom_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

names_to_rule_out = {'20P40_1Month'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);
%% Read DICOM
label = labels{5};
count = 0;
% 17P73, 18P90, 18P92, 18P93, 18P94, 18P95, 20P10, 20P10C,  20P11, 20P40, 20P48
sel_array = [501, 153, 242001, 170, 257001, 172, 270, 267, 264, 202, 201];

whatsinit = cell(length(sel_array), 1);
slice_data = cell(length(sel_array), 1);

for i = 1:length(dicom_glob)

    if i ~= 11
        count = count + 1;
        idx_array = contains(dicom_glob, label);

        ya_glob = glob(cat(2, dicom_glob{i}, 'DICOM/*'));
        dst_names = ExtractNames(ya_glob);
        sub_name = num2str(sel_array(count), '%04.f');
        
        ind_array = zeros(size(dst_names, 1), 1);
        ind_array = ind_array + contains(dst_names, sub_name);
        ind_array2 = find(ind_array == 1);
        list_to_read = ya_glob(ind_array2);

        for j = 1:length(list_to_read)
            % [whatsinit{count}, slice_data{count}] = dicom23D(list_to_read{j});
            yya_glob = glob(cat(2, list_to_read{j}, '*'));
            whatsinit{count} = dicomread(yya_glob{1});
            slice_data{count} = dicominfo(yya_glob{1});
        end
    end
end

rescaled_t2star = cell(length(sel_array), 1);
for i = 1:length(sel_array)
    rescaled_t2star{i} = double(whatsinit{i} .* slice_data{i}.RescaleSlope);
end

%% Main Body (1)
time_label = {'8WK'};
%for n = 1:length(Names)
for n = 8:8
    name = Names{n};
    for t = 1:length(time_label)
        %for t = 1:1
        cvi_glob = glob(cat(2, name_glob{n}, '*.cvi42wsx.xml'));
        if isempty(cvi_glob)
            cvi_glob = glob(cat(2, name_glob{n}, '*.cvi42wsx'));
        end
        if isempty(cvi_glob)
            disp(cat(2, 'NO CVI XML: ', name));
        else
            cvi42wsx = char(cvi_glob);
            con_cell = cell(0);
            for xml_ind = 1:size(cvi42wsx, 1)
                % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
                % Thus, there can be two different files
                con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
            end


            for con_idx = 1:length(con_cell)
                % A different label for Exvivo
                con = con_cell{con_idx};

                for ll = 1:length(sel_array)
                    clear vol_img_3D mask_heart_3D mask_myocardium_3D mask_blood_3D excludeMask_3D myoRefMask_3D noReflowMask_3D freeROIMask_3D volume_image;
                    volume_image = double(whatsinit{ll});
                    slice_data_sg = slice_data{ll};
                    dstFolder = cat(2, base_dir, '/ContourData_Invivo/', name, '/', name, '_', time_label{t}, '/');
                    
                    [mask_heart, mask_myocardium, mask_blood, excludeContour, myoRefCell, noReflowCell, freeROICell, match_count] = ...
                        CMR42ContourMatrixGenerator3(con, volume_image, slice_data_sg, dstFolder);
                    total_match = 0;
                    slc_start = 1;
                    slc_end = 1;
                    slc_array = [];
                    % The following is the same as "ReadCVI_Workflow_Cell_Func.m" in
                    % GUI_for_LRT folder
                    % get all contours from excludeContour
                    excludeMask_3Ds = zeros(size(volume_image));
                    if ~isempty(excludeContour)
                        keys = fieldnames(excludeContour);
                        % The code below assumes 2D slice of image; can be
                        % improved for more      []'
                        % Should be

                        % generic use.
                        temp_mat = zeros(size(volume_image));
                        for j = 1:length(keys)
                            temp_cell = getfield(excludeContour, keys{j});

                            % Add Contour001, Contour002 ...
                            for k = 1:size(temp_cell{1}, 3)
                                temp_mat(:,:,temp_cell{2}(k)) = temp_mat(:,:,temp_cell{2}(k)) + temp_cell{1}(:,:,k);
                            end
                        end
                        excludeMask_3Ds = temp_mat + excludeMask_3Ds;
                    end

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
                            %temp_mat = [];
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
                        if size(volume_image, 3) == 1
                            %slc_end = slc_end + size(volume_image, 3) - 1;
                            % This is for when matrix size doesn't match
                            if i == 1
                                vol_img_3D = {};
                                mask_heart_3D = {};
                                mask_myocardium_3D = {};
                                mask_blood_3D = {};
                                excludeMask_3D = {};
                                myoRefMask_3D = {};
                                noReflowMask_3D = {};
                                freeROIMask_3D = {};
                            end

                            vol_img_3D{slc_start} = volume_image;
                            mask_heart_3D{slc_start} = mask_heart;
                            mask_myocardium_3D{slc_start} = mask_myocardium;
                            mask_blood_3D{slc_start} = mask_blood;
                            excludeMask_3D{slc_start} = excludeMask_3Ds > 0;
                            myoRefMask_3D{slc_start}  = myoRefMask_3Ds;
                            noReflowMask_3D{slc_start}  = noReflowMask_3Ds;
                            freeROIMask_3D{slc_start}  = freeROIMask_3Ds;

                        elseif size(volume_image, 3) > 1 % If read 3D slices, directly apply
                            vol_img_3D = volume_image;
                            mask_heart_3D = mask_heart;
                            mask_myocardium_3D = mask_myocardium;
                            mask_blood_3D = mask_blood;
                            excludeMask_3D = excludeMask_3Ds > 0;
                            myoRefMask_3D  = myoRefMask_3Ds;
                            noReflowMask_3D  = noReflowMask_3Ds;
                            freeROIMask_3D  = freeROIMask_3Ds;
                        end
                        slc_start = slc_start + match_count;
                        slc_end = slc_end + match_count;
                        total_match = total_match + match_count;

                        slc_array = [slc_array, slice_data_sg.SliceLocation];
                    end

                    % Save all as mat file
                    dsts = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'BloodPool', 'freeROI'};

                    % TODO
                    if total_match ~= 0
                        dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
                        save(dstPath, 'vol_img_3D');

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

                        dstPath = cat(2, dstFolder, '/', label, '_SliceLoc.mat');
                        save(dstPath, 'slc_array');

                        disp(cat(2, name, ':   ', label));
                        disp('Done!')
                    end
                end
            end
        end
    end
end

%% Main Body (2)
% Part 1 For varying resolutions (20P40 and 20P48)

name_glob = glob(cat(2, dicom_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

names_to_rule_out = {'17P73', '18P90', '18P92', '18P93', '18P94', '18P95', '20P10', '20P11', '20P40'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

label = labels{5};
% Read DICOM
count = 0;
sel_cell = cell(size(Names, 1), 1); 
sel_cell{1} = [204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240];
sel_cell{2} = [216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248];

whatsinit = cell(size(sel_cell, 1), max(length(sel_cell{1}), length(sel_cell{2})));
slice_data = cell(size(sel_cell, 1), max(length(sel_cell{1}), length(sel_cell{2})));

for i = 1:length(name_glob)
    idx_array = contains(name_glob, label);

    ya_glob = glob(cat(2, name_glob{i}, 'DICOM/*'));
    dst_names = ExtractNames(ya_glob);
    
    temp_array = sel_cell{i};
    for j = 1:length(temp_array)
        
        sub_name = num2str(temp_array(j), '%04.f');

        ind_array = zeros(size(dst_names, 1), 1);
        ind_array = ind_array + contains(dst_names, sub_name);
        ind_array2 = find(ind_array == 1);
        list_to_read = ya_glob(ind_array2);

        % [whatsinit{count}, slice_data{count}] = dicom23D(list_to_read{j});
        yya_glob = glob(cat(2, list_to_read{1}, '*'));
        whatsinit{i, j} = dicomread(yya_glob{1});
        slice_data{i, j} = dicominfo(yya_glob{1});
    end
end

rescaled_t2star_cell = cell(size(sel_cell, 1), max(length(sel_cell{1}), length(sel_cell{2})));
for i = 1:size(rescaled_t2star_cell, 1)
    for j = 1:size(rescaled_t2star_cell, 2)
        if ~(i == 2 && j == 18) && ~(i == 2 && j == 19)
            rescaled_t2star_cell{i,j} = double(whatsinit{i,j} .* slice_data{i,j}.Private_0019_1017);
        end
    end
end

%% Main Body (2)
% part 2

time_label = {'8WK_VaryingRes'};
for n = 1:length(Names)
    name = Names{n};
    for t = 1:length(time_label)
        %for t = 1:1
        cvi_glob = glob(cat(2, name_glob{n}, '*.cvi42wsx.xml'));
        if isempty(cvi_glob)
            cvi_glob = glob(cat(2, name_glob{n}, '*.cvi42wsx'));
        end
        if isempty(cvi_glob)
            disp(cat(2, 'NO CVI XML: ', name));
        else
            cvi42wsx = char(cvi_glob);
            con_cell = cell(0);
            for xml_ind = 1:size(cvi42wsx, 1)
                % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
                % Thus, there can be two different files
                con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
            end


            for con_idx = 1:length(con_cell)
                % A different label for Exvivo
                con = con_cell{con_idx};
                temp_array = sel_cell{n};

                for ll = 1:length(temp_array)
                    clear vol_img_3D mask_heart_3D mask_myocardium_3D mask_blood_3D excludeMask_3D myoRefMask_3D noReflowMask_3D freeROIMask_3D volume_image;
                    volume_image = double(whatsinit{n, ll});
                    slice_data_sg = slice_data{n, ll};

                    strings = strsplit(slice_data_sg.Filename, '/');
                    series_name = strings{end-1};
                    dstFolder = cat(2, base_dir, '/ContourData_Invivo/', name, '/', name, '_', time_label{t}, '/', series_name, '/');
                    
                    [mask_heart, mask_myocardium, mask_blood, excludeContour, myoRefCell, noReflowCell, freeROICell, match_count] = ...
                        CMR42ContourMatrixGenerator3(con, volume_image, slice_data_sg, dstFolder);
                    total_match = 0;
                    slc_start = 1;
                    slc_end = 1;
                    slc_array = [];
                    % The following is the same as "ReadCVI_Workflow_Cell_Func.m" in
                    % GUI_for_LRT folder
                    % get all contours from excludeContour
                    excludeMask_3Ds = zeros(size(volume_image));
                    if ~isempty(excludeContour)
                        keys = fieldnames(excludeContour);
                        % The code below assumes 2D slice of image; can be
                        % improved for more      []'
                        % Should be

                        % generic use.
                        temp_mat = zeros(size(volume_image));
                        for j = 1:length(keys)
                            temp_cell = getfield(excludeContour, keys{j});

                            % Add Contour001, Contour002 ...
                            for k = 1:size(temp_cell{1}, 3)
                                temp_mat(:,:,temp_cell{2}(k)) = temp_mat(:,:,temp_cell{2}(k)) + temp_cell{1}(:,:,k);
                            end
                        end
                        excludeMask_3Ds = temp_mat + excludeMask_3Ds;
                    end

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
                            %temp_mat = [];
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
                        if size(volume_image, 3) == 1
                            %slc_end = slc_end + size(volume_image, 3) - 1;
                            % This is for when matrix size doesn't match
                            if i == 1
                                vol_img_3D = {};
                                mask_heart_3D = {};
                                mask_myocardium_3D = {};
                                mask_blood_3D = {};
                                excludeMask_3D = {};
                                myoRefMask_3D = {};
                                noReflowMask_3D = {};
                                freeROIMask_3D = {};
                            end

                            vol_img_3D{slc_start} = volume_image;
                            mask_heart_3D{slc_start} = mask_heart;
                            mask_myocardium_3D{slc_start} = mask_myocardium;
                            mask_blood_3D{slc_start} = mask_blood;
                            excludeMask_3D{slc_start} = excludeMask_3Ds > 0;
                            myoRefMask_3D{slc_start}  = myoRefMask_3Ds;
                            noReflowMask_3D{slc_start}  = noReflowMask_3Ds;
                            freeROIMask_3D{slc_start}  = freeROIMask_3Ds;

                        elseif size(volume_image, 3) > 1 % If read 3D slices, directly apply
                            vol_img_3D = volume_image;
                            mask_heart_3D = mask_heart;
                            mask_myocardium_3D = mask_myocardium;
                            mask_blood_3D = mask_blood;
                            excludeMask_3D = excludeMask_3Ds > 0;
                            myoRefMask_3D  = myoRefMask_3Ds;
                            noReflowMask_3D  = noReflowMask_3Ds;
                            freeROIMask_3D  = freeROIMask_3Ds;
                        end
                        slc_start = slc_start + match_count;
                        slc_end = slc_end + match_count;
                        total_match = total_match + match_count;

                        slc_array = [slc_array, slice_data_sg.SliceLocation];
                    end

                    % Save all as mat file
                    dsts = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'BloodPool', 'freeROI'};

                    % TODO
                    if total_match ~= 0
                        dstPath = cat(2, dstFolder, '/', label, '_vol_img_3D.mat');
                        save(dstPath, 'vol_img_3D');

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

                        dstPath = cat(2, dstFolder, '/', label, '_SliceLoc.mat');
                        save(dstPath, 'slc_array');

                        disp(cat(2, name, ':   ', label));
                        disp('Done!')
                    end
                end
            end
        end
    end
end

