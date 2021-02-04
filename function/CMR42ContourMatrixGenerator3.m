function [shifted_heart, shifted_myo, shifted_blood, excludeContour, myoRefCell, noReflowCell, freeROICell, match_count] ...
    = CMR42ContourMatrixGenerator3(con, volume_image, slice_data, dstFolder, old_freeROI_label)
% Second version, improved performance. Initially used for CNN
% segmentation.
% Need to work on this to make it 3D compatible? 02/03/2021
% This is hard-coded to solve freeROI coordinate record
if nargin == 4
    old_freeROI_label = 0;
end
% Find Epi- and Endo-cardium from the contours
num_slice = length(slice_data);
num_contours = size(con.contours, 1);
contour_idx = zeros(1, num_slice);

for i = 1:num_slice
    for j = 1:num_contours
        if strcmp(slice_data(i).MediaStorageSOPInstanceUID, con.contours(j).iuid)
            contour_idx(i) = j;
        end
    end
end

match_array = contour_idx > 0;
unmatch_count = 0; % See if slice matches
epi_match_array = zeros(1, num_slice);
endo_match_array = zeros(1, num_slice);

if any(contour_idx(:))
    epi_flow = cell(1, num_slice);
    endo_flow = cell(1, num_slice);
    excludeContour = cell(1, num_slice);
    myoRef = cell(1, num_slice);
    NoReFlow = cell(1, num_slice);
    freeROI = cell(1, num_slice);
    excludeCtr_struct = struct;
    
    %for i = 1:num_contours
    for i = 1:num_slice
        ctype = [];
        if contour_idx(i) ~= 0
            ctype = con.contours(contour_idx(i)).ctype;
        end
        epi_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
        endo_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
        excludeContour{i} = zeros(size(volume_image, 2), size(volume_image,1));
        myoRef{i} = zeros(size(volume_image, 2), size(volume_image,1));
        NoReFlow{i} = zeros(size(volume_image, 2), size(volume_image, 1));
        freeROI{i} = zeros(size(volume_image, 2), size(volume_image, 1));
        for j = 1:length(ctype)
            contour_type = ctype{j};
            if strcmp(contour_type, 'saendocardialContour')
                for c = 1: length(con.contours(contour_idx(i)).pts{j})
                    endo_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                end
            elseif strcmp(contour_type, 'saepicardialContour')
                for c = 1: length(con.contours(contour_idx(i)).pts{j})
                    epi_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                end
            elseif contains(contour_type, 'excludeEnhancementAreaContour')
                for c = 1: length(con.contours(contour_idx(i)).pts{j})
                    excludeContour{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                end
                excludeCtr_struct.(contour_type) = excludeContour;
            elseif contains(contour_type, 'saReferenceMyoContour')
                for c = 1: length(con.contours(contour_idx(i)).pts{j})
                    myoRef{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                end
            elseif contains(contour_type, 'noReflowAreaContour')
                for c = 1: length(con.contours(contour_idx(i)).pts{j})
                    NoReFlow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                end
            elseif contains(contour_type, 'freeDrawRoiContour')
                if old_freeROI_label == 0
                    for c = 1: length(con.contours(contour_idx(i)).pts{j})
                        freeROI{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
                    end
                else
                    for c = 1: length(con.contours(contour_idx(i)).pts{j})
                        freeROI{i}([round(con.contours(contour_idx(i)).pts{j}(c,1)/2)], [round(con.contours(contour_idx(i)).pts{j}(c,2)/2)]) = 1;
                    end
                end
            end
        end
    end
    
    % Remove zero matrix
    count = 1;
    index_array = [];
    
    for i = 1 : num_slice
        if any(epi_flow{i}(:)) && any(endo_flow{i}(:))
            epi_flow_edit(:,:,count) = epi_flow{i}(:,:);
            endo_flow_edit(:,:,count) = endo_flow{i}(:,:);
            index_array = [index_array; i];
            count = count + 1;
            
            epi_match_array(i) = 1;
            endo_match_array(i) = 1;
        end
    end
    
    if exist('epi_flow_edit', 'var') == 1 || exist('endo_flow_edit', 'var') == 1
        %%%%% Cases with no epi/endo myocardium such as there's is only reference point
        %%%%% using OR rather than and because I still want error pops up
        %%%%% when there is missing epi or endo - cardium contour
        %%%%% (Different situation)
        
        
        % For excludeContour (Convert to matrix)
        clear excludeContour_edit
        excludeContour_edit = {};
        fname = fieldnames(excludeCtr_struct); % Considering multiple exclusion contours
        
        for i = 1 : length(fname)
            count = 1;
            
            ctr_index_array = []; %%%%% Not necessary
            for j = 1 : num_slice
                if any(excludeCtr_struct.(fname{i}){j}(:))
                    excludeCtr_mat(:,:,count) =  excludeCtr_struct.(fname{i}){j};
                    ctr_index_array = [ctr_index_array; j];
                    count = count + 1;
                end
            end
            excludeContour_edit(i, 1:2) = {fname{i}, excludeCtr_mat};
            excludeContour_edit(i, 3) = {ctr_index_array};
        end
        
        % For myoRef (Convert to matrix)
        count = 1;
        myo_index_array = []; %%%%% NOT necessary
        myoRef_edit = [];
        for i = 1 : num_slice
            if any(myoRef{i}(:))
                myoRef_edit(:,:,count) = myoRef{i}(:,:);
                myo_index_array = [myo_index_array; i];
                count = count + 1;
            end
        end
                
        
        % For NoReFlow (Convert to matrix)
        count = 1;
        noRef_index_array = [];
        NoReFlow_edit = [];
        for i = 1:num_slice
            if any(NoReFlow{i}(:))
                NoReFlow_edit(:,:,count) = NoReFlow{i}(:,:);
                noRef_index_array = [noRef_index_array; i];
                count = count + 1;
            end
        end
        
        count = 1;
        freeROI_index_array = [];
        freeROI_edit = [];
        for i = 1 : num_slice
            if any(freeROI{i}(:))
                freeROI_edit(:,:,count) = freeROI{i}(:,:);
                freeROI_index_array = [freeROI_index_array; i];
                count = count + 1;
            end
        end
        
        
        
        
        % Apply heart to the image and save it to OutputPath
        if ~ exist(dstFolder, 'dir')
            mkdir(dstFolder);
        end
        
        % Get center coordinate and region growing
        contours = epi_flow_edit | endo_flow_edit;
        num = size(contours, 3);
        myocardium = zeros(size(volume_image,2), size(volume_image,1), num);
        heart = zeros(size(volume_image,2), size(volume_image,1), num);
        blood_pool = zeros(size(volume_image,2), size(volume_image,1), num);
        center_array = zeros(num, 2);
        for i = 1: num
            [x_epi, y_epi] = find(epi_flow_edit(:,:,i) ~= 0);
            x_epi_center = round(mean(x_epi));
            y_epi_center = round(mean(y_epi));
            
            epi_J = regiongrowing(epi_flow_edit(:,:,i),x_epi_center,y_epi_center);
            epi_J = epi_J | epi_flow_edit(:,:,i); % Modified
            
            [x_endo, y_endo] = find(endo_flow_edit(:,:,i) ~= 0);
            x_endo_center = round(mean(x_endo));
            y_endo_center = round(mean(y_endo));
            
            endo_J = regiongrowing(endo_flow_edit(:,:,i),x_endo_center,y_endo_center);
            
            myocardium(:,:,i) = xor(epi_J, endo_J);
            heart(:,:,i) = epi_J;
            blood_pool(:,:,i) = endo_J;
            center_array(i, :) = [(x_epi_center + x_endo_center)/2, (y_epi_center + y_endo_center)/2];
        end
        
        shifted_heart = zeros(size(volume_image));
        shifted_myo = zeros(size(volume_image));
        shifted_blood = zeros(size(volume_image));
        
        for i = 1:num
            shifted_heart(:,:,i) = flipud(rot90(heart(:,:,i),1));
            shifted_myo(:,:,i) = flipud(rot90(myocardium(:,:,i),1));
            shifted_heart(:,:,i) = imfill(shifted_heart(:,:,i), 'holes');
            mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,index_array(i));
            mask_myocardium = shifted_myo(:,:,i) .* volume_image(:,:,index_array(i));
            
            shifted_blood(:,:,i) = flipud(rot90(blood_pool(:,:,i),1));
            shifted_blood(:,:,i) = imfill(shifted_blood(:,:,i), 'holes');
            mask_blood = shifted_blood(:,:,i) .* volume_image(:,:,index_array(i));
            
        end
        
        % Export Exclude Contour
        excludeContour = struct;
        if ~isempty(excludeContour_edit)
            for i = 1:length(fname)
                excludeContour.(fname{i}) = cell(1, 2); % No index array needed
                for j = 1:size(excludeContour_edit{i, 2}, 3)
                    excludeContour_edit{i, 2}(:,:,j) = imfill(excludeContour_edit{i, 2}(:,:,j), 'holes');
                    excludeContour_edit{i, 4}(:,:,j) = flipud(rot90(excludeContour_edit{i, 2}(:,:,j),1));
                end
                excludeContour.(fname{i}){1} = excludeContour_edit{i, 4};
                excludeContour.(fname{i}){2} = excludeContour_edit{i, 3}; % No
                % need for index array
            end
        end
        
        
        % Export Refernce area
        myoRefCell = cell(1, 2);
        if ~isempty(myoRef_edit)
            for i = 1:size(myoRef_edit, 3)
                myoRef_edit(:,:,i) = imfill(myoRef_edit(:,:,i), 'holes');
                shifted_myoRef(:,:,i) = flipud(rot90(myoRef_edit(:,:,i),1));
            end
            myoRefCell{1} = shifted_myoRef;
            myoRefCell{2} = myo_index_array;
            % Can quickly add myo_index_array back
        end
        
        
        
        % Export No Reflow
        % if strcmp(dstLabel, 'noReflowAreaContour') && ~isempty(noRef_index_array)
        noReflowCell = cell(1, 2);
        if ~isempty(NoReFlow_edit)
            for i = 1:size(NoReFlow_edit, 3)
                NoReFlow_edit(:,:,i) = imfill(NoReFlow_edit(:,:,i), 'holes');
                shifted_NoReFlow(:,:,i) = flipud(rot90(NoReFlow_edit(:,:,i),1));
            end
            noReflowCell{1} = shifted_NoReFlow;
            noReflowCell{2} = noRef_index_array;
        end
        
        % Export freeROI
        freeROICell = cell(1, 2);
        if ~isempty(freeROI_edit)
            for i = 1:size(freeROI_edit, 3)
                freeROI_edit(:,:,i) = imfill(freeROI_edit(:,:,i), 'holes');
                shifted_freeROI(:,:,i) = flipud(rot90(freeROI_edit(:,:,i),1));
            end
            freeROICell{1} = shifted_freeROI;
            freeROICell{2} = freeROI_index_array;
        end
    else
        % unmatch_count = unmatch_count + 1;
    end
    
else
   % unmatch_count = unmatch_count + 1;
    
    % disp('No Match');
end

match_count = sum(epi_match_array & endo_match_array & match_array);

% match_count = (num_slice>0) - unmatch_count;
if match_count == 0
    mask_heart = [];
    mask_myocardium = [];
    mask_blood = [];
    excludeContour = [];
    myoRefCell = [];
    noReflowCell = [];
    freeROICell = [];
    shifted_heart = [];
    shifted_myo = [];
    shifted_blood = [];
end
end