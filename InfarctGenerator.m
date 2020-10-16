clear all;
close all;
clc;
%% 
addpath('function\');
% I'm referring ContourData in ShareDrive/Yinyin
%base_dir = GetFullPath(cat(2, current_dir, '\..\ContourData\'));
base_dir = uigetdir;
img_dir = GetFullPath(cat(2, base_dir, '/Data_Patient/'));
contour_dir = GetFullPath(cat(2, base_dir, '/ContourData/'));

sequence_label = {'MAG', 'PSIR', 'T1Map', 'T1w_MOCO'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea'};
output_label = {'LGE', 'T1'};

name_glob = glob(cat(2, img_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

names_to_rule_out = {'RUTHIE', '11D27'};

RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);


infarct_perc_struct = struct;

infarct_perc_array = zeros(size(Names));


%label = sequence_label{3}; %%%
for la = 1:length(sequence_label)-1
%for la = 1:1
for i = 1:length(Names)
%for i = 1:1
    name = Names{i};
    disp(name)
    label = sequence_label{la};
    switch label
        case {'MAG', 'PSIR'}
            labelo = output_label{1};
        otherwise
            labelo = output_label{2};
    end
    
    clear vi3 svi3
    % Assume images must be present
    load(cat(2, contour_dir, name, '/', labelo, '/', label, '_vol_img_3D.mat'));
    if exist('vi3', 'var')
        img = vi3.vol_img_3D;
    else
        img = svi3.sig_vol_img_3D;
    end

    mask_myocardium_3D = [];
    excludeMask_3D = [];
    myoRefMask_3D = [];
    noReflowMask_3D = [];
    
    if length(ls(cat(2, contour_dir,'/', name, '/', labelo, '/', anatomy_label{2}))) > 2
        load(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{2}, '/', 'mask_myocardium.mat'));
    end
    if length(ls(cat(2, contour_dir,'/', name, '/', labelo, '/', anatomy_label{3}))) > 2
        load(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{3}, '/', 'excludeArea.mat'));
    end
    if length(ls(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{4}))) > 2
        load(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{4}, '/', 'myoRef.mat'));
    end
    if length(ls(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{5}))) > 2
        load(cat(2, contour_dir, name, '/', labelo, '/', anatomy_label{5}, '/', 'noReflow.mat'));
    end
   
    % Loaded are all binary masks with double floating numer type with same
    % size (should be)
    
    % Load image...
    
    compositeIm = mask_myocardium_3D + myoRefMask_3D + (-5) * excludeMask_3D + 3 * noReflowMask_3D;
    Outpath = cat(2, contour_dir, name, '/', labelo, '/', 'compositeLabel.mat');
    % Always overwrite
    save(Outpath, 'compositeIm');
    
    
    refIm = (compositeIm == -3) | (compositeIm == 2);
    exIm = compositeIm < 0;
    noReflowIm = compositeIm == 4;
    % Keep them 3D
    
    %figure();
    slice_array = linspace(1, size(mask_myocardium_3D, 3), size(mask_myocardium_3D, 3));
    infarct_raw_3D = zeros(size(img));
    for slc = 1:length(slice_array)
        %subplot(2,2,slc);
        %imagesc(compositeIm(:,:,slc))
                                   
        if any(refIm(:,:,slc))
            slcr = slc;
        else
            dist_array = abs(slice_array - slc);
            start_dist = 1;
            flag = 0;
            for s = 1:length(unique(dist_array)) - 1
                idx = find(dist_array == start_dist);
                for dx = 1:length(idx)
                    refIm_slc = refIm(:,:,idx(dx));
                    if any(refIm_slc(:))
                        slcr = idx(dx);
                        flag = 1;
                        break;
                    end
                end
                
                if flag == 1
                    break;
                end
                start_dist = start_dist + 1;
            end
        end
        refMasked = refIm(:,:,slcr) .* img(:,:,slcr); 
        ref_mean = mean(nonzeros(refMasked(:)));
        ref_sd = std(nonzeros(refMasked(:)));
        thresh = ref_mean + 5*ref_sd;
        
        infarct_raw = img(:,:,slc) > thresh;
        
        infarct_raw_3D(:,:,slc) = infarct_raw;
    end
    
    % Include excludeArea (Without Hemorrhage)
    infarct_masked = infarct_raw_3D .* mask_myocardium_3D;
    infarct_ex_masked = infarct_masked .* (1-exIm);
    
    % Apply noReflowArea (non-Hemorrhage + Hemorrhage)
    infarct_whole_masked = noReflowIm + infarct_ex_masked;
    infarct_whole_masked(infarct_whole_masked > 1) = 1; % Make sure infarct mask is binary
    
    infarct_ex_masked = ImPostProc(infarct_ex_masked);
    infarct_whole_masked = ImPostProc(infarct_whole_masked);
    
    infarct_noReflow_masked = infarct_whole_masked - infarct_ex_masked;
    
    infarct_out = cat(2, contour_dir, name, '/', labelo, '/', label, '_MI/');
    if ~ exist(infarct_out, 'dir')
        mkdir(infarct_out);
    end
    
    infarct_out_path = cat(2, infarct_out, 'Infarct_Whole.mat');
    save(infarct_out_path, 'infarct_whole_masked');
    
    if strcmp(label, sequence_label{1}) || strcmp(label, sequence_label{2})
        infarct_out_path = cat(2, infarct_out, 'Infarct_BeforeNoReflowArea.mat');
        save(infarct_out_path, 'infarct_ex_masked');
        
        infarct_out_path = cat(2, infarct_out, 'Infarct_AfterNoReflowArea.mat');
        save(infarct_out_path, 'infarct_noReflow_masked');
    else
        infarct_out_path = cat(2, infarct_out, 'Infarct_NonHemo.mat');
        save(infarct_out_path, 'infarct_ex_masked');
        
        infarct_out_path = cat(2, infarct_out, 'Infarct_Hemo.mat');
        save(infarct_out_path, 'infarct_noReflow_masked');
    end
end
end

%%
figure();
for i = 1:size(infarct_whole_masked, 3)
    subplot(3,3,i)
    imagesc(infarct_whole_masked(:,:,i)); axis equal;
end
%%
figure();
for i = 1:size(infarct_whole_masked, 3)
    subplot(3,3,i)
    imagesc(compositeIm(:,:,i)); axis equal;
end
%%
figure();
for i = 1:size(infarct_whole_masked, 3)
    subplot(3,3,i)
    imagesc(excludeMask_3D(:,:,i)); axis equal;
end
%%
figure();
for i = 1:size(infarct_whole_masked, 3)
    subplot(3,3,i)
    imagesc(infarct_ex_masked(:,:,i)); axis equal;
end

%%
figure();
for i = 1:size(infarct_whole_masked, 3)
    subplot(3,3,i)
    imagesc(myoRefMask_3D(:,:,i)); axis equal;
end