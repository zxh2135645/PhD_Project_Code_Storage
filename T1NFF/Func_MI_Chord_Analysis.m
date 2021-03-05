function [MI_Chord_Analysis, center_mask_ff] = ...
    Func_MI_Chord_Analysis(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, ...
    MI_Chord_Analysis_fname)
Mipix = cell(Segn, size(t1, 3));
Mipix_epi = cell(Segn, size(t1, 3));
Mipix_endo = cell(Segn, size(t1, 3));
Mipix2 = cell(Segn, size(t1, 3));
Mipix2_epi = cell(Segn, size(t1, 3));
Mipix2_endo = cell(Segn, size(t1, 3));
Mipix3 = cell(Segn, size(t1, 3));
Mipix3_epi = cell(Segn, size(t1, 3));
Mipix3_endo = cell(Segn, size(t1, 3));

Mipix_mean = [];
Mipix_mean2 = [];
Mipix_mean3 = [];
Mipix_mean_epi = [];
Mipix_mean2_epi = [];
Mipix_mean3_epi = [];
Mipix_mean_endo = [];
Mipix_mean2_endo = [];
Mipix_mean3_endo = [];

se = strel('disk', 1);
[optimizer, metric] = imregconfig('multimodal');
MI_Chord_Analysis = struct;
center_mask_ff = zeros(size(roi_in_myo_t1));
center_mask_t1 = zeros(size(roi_in_myo_t1));
for i = 1:size(t1, 3)
    %for i = 4:4
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    img3 = r2star(:,:,i);
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    img2 = imwarp(img2,tform,'OutputView',imref2d(size(fixed)));
    img3 = imwarp(img3,tform,'OutputView',imref2d(size(fixed)));

    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
    
    fixed_eroded = imerode(fixed, se);
    BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
    center_mask_t1(:,:,i) = imfill(BW_skel, 'hole');
    center_fixed = center_mask_t1(:,:,i);
    fixedRegistered_epi = fixed_eroded - center_fixed > 0;
    fixedRegistered_endo = center_fixed + fixed_eroded > 1;
    
    movingRegistered_eroded = imerode(movingRegistered, se);
    BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
    center_mask_ff(:,:,i) = imfill(BW_skel, 'hole');
    center_moving = center_mask_ff(:,:,i);
    movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
    movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
    
    movingRegistered_roi_ff = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    %movingRegistered_epi = imwarp(myo_epi_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    %movingRegistered_endo = imwarp(myo_endo_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(img, fixed_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, movingRegistered_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, movingRegistered_eroded, Segn, Groove);
    
    [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(img, fixedRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(img, fixedRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_epi] = AHASegmentation(img2, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_endo] = AHASegmentation(img2, movingRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_epi] = AHASegmentation(img3, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_endo] = AHASegmentation(img3, movingRegistered_endo, Segn, Groove);
        
    for j = 1:Segn
        Mipix{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        Mipix_epi{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        Mipix_endo{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j);
        
        Mipix2{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        Mipix2_epi{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        Mipix2_endo{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded == j);
        
        Mipix3{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
        Mipix3_epi{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
        Mipix3_endo{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded == j);
    end
    
    Mi_idx1 = [];
    Mi_idx1_epi = [];
    Mi_idx1_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix{j,i})
            Mi_idx1 = [Mi_idx1, j];
        end
        if ~isempty(Mipix_epi{j,i})
            Mi_idx1_epi = [Mi_idx1_epi, j];
        end
        if ~isempty(Mipix_endo{j,i})
            Mi_idx1_endo = [Mi_idx1_endo, j];
        end
    end
    
    Mi_idx2 = [];
    Mi_idx2_epi = [];
    Mi_idx2_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix2{j,i})
            Mi_idx2 = [Mi_idx2, j];
        end
        if ~isempty(Mipix2_epi{j,i})
            Mi_idx2_epi = [Mi_idx2_epi, j];
        end
        if ~isempty(Mipix_endo{j,i})
            Mi_idx2_endo = [Mi_idx2_endo, j];
        end
    end
    
    Mi_idx3 = [];
    Mi_idx3_epi = [];
    Mi_idx3_endo = [];
    for j = 1:Segn
        if ~isempty(Mipix3{j,i})
            Mi_idx3 = [Mi_idx3, j];
        end
        if ~isempty(Mipix3_epi{j,i})
            Mi_idx3_epi = [Mi_idx3_epi, j];
        end
        if ~isempty(Mipix3_endo{j,i})
            Mi_idx3_endo = [Mi_idx3_endo, j];
        end
    end
    
     
    Mi_idx_intersect = intersect(Mi_idx1, Mi_idx2);
    Mi_idx_intersect_epi = intersect(Mi_idx1_epi, Mi_idx2_epi);
    Mi_idx_intersect_endo = intersect(Mi_idx1_endo, Mi_idx2_endo);
    
    Mipix_flat = {};
    Mipix_flat2 = {};
    Mipix_flat3 = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect)
        ind = Mi_idx_intersect(j);
        Mipix_flat{idx} = Mipix{ind,i};
        Mipix_flat2{idx} = Mipix2{ind,i};
        Mipix_flat2{idx}(Mipix_flat2{idx} > 100) = 100;
        Mipix_flat2{idx}(Mipix_flat2{idx} < 0) = 0;
        
        Mipix_flat3{idx} = Mipix3{ind,i};
        Mipix_flat3{idx}(Mipix_flat3{idx} < 0) = 0;
        
        Mipix_mean(idx, i) = mean(Mipix_flat{idx});
        Mipix_mean2(idx, i) = mean(Mipix_flat2{idx});
        Mipix_mean3(idx, i) = mean(Mipix_flat3{idx});
        idx = idx + 1;
    end
    
    Mipix_flat_epi = {};
    Mipix_flat2_epi = {};
    Mipix_flat3_epi = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect_epi)
        ind = Mi_idx_intersect_epi(j);
        Mipix_flat_epi{idx} = Mipix_epi{ind,i};
        Mipix_flat2_epi{idx} = Mipix2_epi{ind,i};
        Mipix_flat2_epi{idx}(Mipix_flat2_epi{idx} > 100) = 100;
        Mipix_flat2_epi{idx}(Mipix_flat2_epi{idx} < 0) = 0;
        
        Mipix_flat3_epi{idx} = Mipix3_epi{ind,i};
        Mipix_flat3_epi{idx}(Mipix_flat3_epi{idx} < 0) = 0;
        
        Mipix_mean_epi(idx, i) = mean(Mipix_flat_epi{idx});
        Mipix_mean2_epi(idx, i) = mean(Mipix_flat2_epi{idx});
        Mipix_mean3_epi(idx, i) = mean(Mipix_flat3_epi{idx});
        idx = idx + 1;
    end
    
    Mipix_flat_endo = {};
    Mipix_flat2_endo = {};
    Mipix_flat3_endo = {};
    
    idx = 1;
    for j = 1:length(Mi_idx_intersect_endo)
        ind = Mi_idx_intersect_endo(j);
        Mipix_flat_endo{idx} = Mipix_endo{ind,i};
        Mipix_flat2_endo{idx} = Mipix2_endo{ind,i};
        Mipix_flat2_endo{idx}(Mipix_flat2_endo{idx} > 100) = 100;
        Mipix_flat2_endo{idx}(Mipix_flat2_endo{idx} < 0) = 0;
        
        Mipix_flat3_endo{idx} = Mipix3_endo{ind,i};
        Mipix_flat3_endo{idx}(Mipix_flat3_endo{idx} < 0) = 0;
        
        Mipix_mean_endo(idx, i) = mean(Mipix_flat_endo{idx});
        Mipix_mean2_endo(idx, i) = mean(Mipix_flat2_endo{idx});
        Mipix_mean3_endo(idx, i) = mean(Mipix_flat3_endo{idx});
        idx = idx + 1;
    end
    
    MI_Chord_Analysis(i).Mipix = Mipix; 
    MI_Chord_Analysis(i).Mipix2 = Mipix2; 
    MI_Chord_Analysis(i).Mipix3 = Mipix3;
    MI_Chord_Analysis(i).Mipix_mean = Mipix_mean; 
    MI_Chord_Analysis(i).Mipix_mean2 = Mipix_mean2; 
    MI_Chord_Analysis(i).Mipix_mean3 = Mipix_mean3;
    MI_Chord_Analysis(i).Mask_Segn = Mask_Segn; 
    MI_Chord_Analysis(i).Mask_Segn2 = Mask_Segn2; 
    MI_Chord_Analysis(i).Mask_Segn3 = Mask_Segn3;
    
    MI_Chord_Analysis(i).Mipix_epi = Mipix_epi; 
    MI_Chord_Analysis(i).Mipix2_epi = Mipix2_epi; 
    MI_Chord_Analysis(i).Mipix3_epi = Mipix3_epi;
    MI_Chord_Analysis(i).Mipix_mean_epi = Mipix_mean_epi; 
    MI_Chord_Analysis(i).Mipix_mean2_epi = Mipix_mean2_epi; 
    MI_Chord_Analysis(i).Mipix_mean3_epi = Mipix_mean3_epi;
    MI_Chord_Analysis(i).Mask_Segn_epi = Mask_Segn_epi; 
    MI_Chord_Analysis(i).Mask_Segn2_epi = Mask_Segn2_epi; 
    MI_Chord_Analysis(i).Mask_Segn3_epi = Mask_Segn3_epi;
    
    MI_Chord_Analysis(i).Mipix_endo = Mipix_endo; 
    MI_Chord_Analysis(i).Mipix2_endo = Mipix2_endo; 
    MI_Chord_Analysis(i).Mipix3_endo = Mipix3_endo;
    MI_Chord_Analysis(i).Mipix_mean_endo = Mipix_mean_endo; 
    MI_Chord_Analysis(i).Mipix_mean2_endo = Mipix_mean2_endo; 
    MI_Chord_Analysis(i).Mipix_mean3_endo = Mipix_mean3_endo;
    MI_Chord_Analysis(i).Mask_Segn_endo = Mask_Segn_endo; 
    MI_Chord_Analysis(i).Mask_Segn2_endo = Mask_Segn2_endo; 
    MI_Chord_Analysis(i).Mask_Segn3_endo = Mask_Segn3_endo;
    
    MI_Chord_Analysis(i).fixed_eroded = fixed_eroded;
    MI_Chord_Analysis(i).movingRegistered_eroded = movingRegistered_eroded;
end
    save(MI_Chord_Analysis_fname, 'MI_Chord_Analysis');
end