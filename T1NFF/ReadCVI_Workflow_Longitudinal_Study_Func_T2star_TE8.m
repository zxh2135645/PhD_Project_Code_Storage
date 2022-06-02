% A generic function for reading and parsing cvi42wsx files 
% And generate all its masks
% For - T1_fat_Project

% The size of 90min does not match the rest
%% This should be deprecated (02/25/2022)
function ReadCVI_Workflow_Longitudinal_Study_Func_T2star_TE8(con, dicom_glob, dstFolder, dicom_fields, old_freeROI_label)

if nargin == 4
    old_freeROI_label = 0; % What does this mean?
    % This is hard-coded to solve freeROI coordinate record
end
% if ispc
%     strings = strsplit(dicom_glob{1}, '\');
% else
%     strings = strsplit(dicom_glob{1}, '/');
% end
strings = strsplit(dstFolder, '/');
name = strings{end-2};
label = strings{end-1};
t2star_labels = {'T2star', 'T2star_CMR', 'T2star_SIEMENS'};
if strcmp(label, 'MAG') || strcmp(label, 'PSIR')
    labelo = 'LGE';
else
    labelo = label;
end
            dicom = char(dicom_glob); % Strings should be the same length, otherwise recategorize
            id_cell = cell(size(dicom, 1), 1);
            total_match = 0;
            slc_array = [];
            if any(strcmp(label, t2star_labels))
                glob_idx = [];
            end
            slc_start = 1;
            slc_end = 1;
            clear vol_img_3D mask_heart_3D mask_myocardium_3D mask_blood_3D excludeMask_3D myoRefMask_3D noReflowMask_3D freeROIMask_3D volume_image
            for i = 1:size(dicom, 1)
                dicom_file = glob(cat(2, dicom(i,:), '*.dcm'));
                if isempty(dicom_file)
                    
                    dicom_file = glob(cat(2, dicom(i,:), '*.IMA'));
                    
                    if isempty(dicom_file)
                        
                        dicom_file = {dicom(i,:)};
                    end

                end
                
                dicom_f = dicom_file{1};
                info = dicominfo(dicom_f);
                
                if contains(info.InstitutionName, 'Vida')
                    
                    [volume_image, slice_data] = convert_vida_header(dicom_f, info);
                else
                    dicom_string = dicom(i,:);
                    dicom_string = dicom_string(~isspace(dicom_string));
                    [volume_image, slice_data, image_meta_data] = dicom23D(dicom_string, dicom_fields);
                end
                
                
                if any(strcmp(label, t2star_labels))
                    sloc_array = zeros(length(slice_data), 1);
                    for ss = 1:length(slice_data)
                        sloc_array(ss) = slice_data(ss).SliceLocation;
                    end
                    % Reshape volume_image
                    sz = size(volume_image);
                    slc = length(unique(sloc_array));
                    nte = length(slice_data)/slc;
                    volume_image = reshape(volume_image, sz(1), sz(2), nte, slc);
                    
                    echo_idx = 1;
                    % echo_idx_te8 = 8;
                    volume_image = squeeze(volume_image(:,:,echo_idx, :));
                    slice_data = slice_data(echo_idx:nte:end);
                    % slice_data = slice_data(echo_idx);
                end
                % id_cell{i} = slice_data.MediaStorageSOPInstanceUID; % Didn't do anything to it
                [mask_heart, mask_myocardium, mask_blood, excludeContour, myoRefCell, noReflowCell, freeROICell, match_count] = ...
                    CMR42ContourMatrixGenerator3(con, volume_image, slice_data, dstFolder, old_freeROI_label);
                
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
                    
                    vol_img_3D(:,:,slc_start:slc_end+match_count-1) = volume_image;
                    mask_heart_3D(:,:,slc_start:slc_end+match_count-1) = mask_heart;
                    mask_myocardium_3D(:,:,slc_start:slc_end+match_count-1) = mask_myocardium;
                    mask_blood_3D(:,:,slc_start:slc_end+match_count-1) = mask_blood;
                    excludeMask_3D(:,:,slc_start:slc_end+match_count-1) = excludeMask_3Ds > 0;
                    myoRefMask_3D(:,:,slc_start:slc_end+match_count-1)  = myoRefMask_3Ds;
                    noReflowMask_3D(:,:,slc_start:slc_end+match_count-1)  = noReflowMask_3Ds;
                    freeROIMask_3D(:,:,slc_start:slc_end+match_count-1)  = freeROIMask_3Ds;
                    
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
                    
                    slc_array = [slc_array, slice_data.SliceLocation];
                    if any(strcmp(label, t2star_labels))
                        glob_idx = [glob_idx, i];
                    end
                end
                
            end
            
            %% Save all as mat file
            %  dsts = {'Heart', 'Myocardium', 'excludeContour', 'myoReference', 'noReflowAreaContour', 'BloodPool'};
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
                
                dstPath = cat(2, dstFolder, '/', labelo, '_SliceLoc.mat');
                save(dstPath, 'slc_array');
                
                if any(strcmp(label, t2star_labels))
                    glob_names = cell(1, length(glob_idx));
                    for i = 1:length(glob_names)
                        dicom = dicom_glob{glob_idx(i)};
                        strings = strsplit(dicom, '/');
                        glob_names{i} = strings{end-1};
                    end
                    dstPath = cat(2, dstFolder, '/', label, '_Index.mat');
                    save(dstPath, 'glob_names');
                end
            
                disp(cat(2, name, ':   ', label));
                disp('Done!')
            end

end