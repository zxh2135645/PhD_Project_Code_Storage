function [MI_Chord_Analysis2, center_mask_ff] = ...
    Func_MI_Chord_Analysis2(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, ...
    MI_Chord_Analysis_fname)

% Making ff10 masks and R2*80 masks
Mipix_hemo_p = cell(Segn, size(t1, 3));
Mipix_hemo_n = cell(Segn, size(t1, 3));
Mipix_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix_epi_hemo_p = cell(Segn, size(t1, 3));
Mipix_epi_hemo_n = cell(Segn, size(t1, 3));
Mipix_epi_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix_epi_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix_endo_hemo_p = cell(Segn, size(t1, 3));
Mipix_endo_hemo_n = cell(Segn, size(t1, 3));
Mipix_endo_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix_endo_hemo_n_ff_n = cell(Segn, size(t1, 3));

Mipix2_hemo_p = cell(Segn, size(t1, 3));
Mipix2_hemo_n = cell(Segn, size(t1, 3));
Mipix2_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix2_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix2_epi_hemo_p = cell(Segn, size(t1, 3));
Mipix2_epi_hemo_n = cell(Segn, size(t1, 3));
Mipix2_epi_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix2_epi_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix2_endo_hemo_p = cell(Segn, size(t1, 3));
Mipix2_endo_hemo_n = cell(Segn, size(t1, 3));
Mipix2_endo_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix2_endo_hemo_n_ff_n = cell(Segn, size(t1, 3));

Mipix3_hemo_p = cell(Segn, size(t1, 3));
Mipix3_hemo_n = cell(Segn, size(t1, 3));
Mipix3_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix3_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix3_epi_hemo_p = cell(Segn, size(t1, 3));
Mipix3_epi_hemo_n = cell(Segn, size(t1, 3));
Mipix3_epi_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix3_epi_hemo_n_ff_n = cell(Segn, size(t1, 3));
Mipix3_endo_hemo_p = cell(Segn, size(t1, 3));
Mipix3_endo_hemo_n = cell(Segn, size(t1, 3));
Mipix3_endo_hemo_n_ff_p = cell(Segn, size(t1, 3));
Mipix3_endo_hemo_n_ff_n = cell(Segn, size(t1, 3));

Mipix_mean_hemo_p = [];
Mipix_mean2_hemo_p = [];
Mipix_mean3_hemo_p = [];
Mipix_mean_hemo_n = [];
Mipix_mean2_hemo_n = [];
Mipix_mean3_hemo_n = [];
Mipix_mean_hemo_n_ff_p = [];
Mipix_mean2_hemo_n_ff_p = [];
Mipix_mean3_hemo_n_ff_p = [];
Mipix_mean_hemo_n_ff_n = [];
Mipix_mean2_hemo_n_ff_n = [];
Mipix_mean3_hemo_n_ff_n = [];

Mipix_mean_endo_hemo_p = [];
Mipix_mean2_endo_hemo_p = [];
Mipix_mean3_endo_hemo_p = [];
Mipix_mean_endo_hemo_n = [];
Mipix_mean2_endo_hemo_n = [];
Mipix_mean3_endo_hemo_n = [];
Mipix_mean_endo_hemo_n_ff_p = [];
Mipix_mean2_endo_hemo_n_ff_p = [];
Mipix_mean3_endo_hemo_n_ff_p = [];
Mipix_mean_endo_hemo_n_ff_n = [];
Mipix_mean2_endo_hemo_n_ff_n = [];
Mipix_mean3_endo_hemo_n_ff_n = [];

Mipix_mean_epi_hemo_p = [];
Mipix_mean2_epi_hemo_p = [];
Mipix_mean3_epi_hemo_p = [];
Mipix_mean_epi_hemo_n = [];
Mipix_mean2_epi_hemo_n = [];
Mipix_mean3_epi_hemo_n = [];
Mipix_mean_epi_hemo_n_ff_p = [];
Mipix_mean2_epi_hemo_n_ff_p = [];
Mipix_mean3_epi_hemo_n_ff_p = [];
Mipix_mean_epi_hemo_n_ff_n = [];
Mipix_mean2_epi_hemo_n_ff_n = [];
Mipix_mean3_epi_hemo_n_ff_n = [];

se = strel('disk', 1);
[optimizer, metric] = imregconfig('multimodal');
MI_Chord_Analysis2 = struct;
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
        
    % FF >< 10%
    % r2star >< 80 s-1
    % Registration issue
    movingRegistered_ff10 = img2 .* movingRegistered_roi_ff >= 10;
    movingRegistered_r2star80 = img3 .* movingRegistered_roi_r2star >= 80;
    movingRegistered_ffsub10 = movingRegistered_roi_ff - movingRegistered_ff10;
    movingRegistered_r2starsub80 = movingRegistered_roi_r2star - movingRegistered_r2star80;
    
    % movingRegistered_r2stars80 - as hemo+ mask for whatever the ff is
    % movingRegistered_r2starsub80 .* movingRegistered_ff10 hemo- + ff >= 10
    % movingRegistered_r2star80 .* movingRegistered_ffsub10 hemo- + ff < 10
    hemo_p = movingRegistered_r2star80;
    hemo_n = movingRegistered_r2starsub80;
    hemo_n_ff_p = movingRegistered_r2starsub80 .* movingRegistered_ff10;
    hemo_n_ff_n = movingRegistered_r2starsub80 .* movingRegistered_ffsub10;
    
    for j = 1:Segn
        Mipix_hemo_p{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_p == j);
        Mipix_epi_hemo_p{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_p == j);
        Mipix_endo_hemo_p{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_p == j);
        
        Mipix_hemo_n{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n == j);
        Mipix_epi_hemo_n{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n == j);
        Mipix_endo_hemo_n{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n == j);
        
        Mipix_hemo_n_ff_p{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_p == j);
        Mipix_epi_hemo_n_ff_p{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_p == j);
        Mipix_endo_hemo_n_ff_p{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_p == j);
        
        Mipix_hemo_n_ff_n{j,i} = img(Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_n == j);
        Mipix_epi_hemo_n_ff_n{j,i} = img(Mask_Segn_epi .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_n == j);
        Mipix_endo_hemo_n_ff_n{j,i} = img(Mask_Segn_endo .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_n == j);
        
        Mipix2_hemo_p{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_p == j);
        Mipix2_epi_hemo_p{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_p == j);
        Mipix2_endo_hemo_p{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_p == j);
        
        Mipix2_hemo_n{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n == j);
        Mipix2_epi_hemo_n{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n == j);
        Mipix2_endo_hemo_n{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n == j);
        
        Mipix2_hemo_n_ff_p{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_p == j);
        Mipix2_epi_hemo_n_ff_p{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_p == j);
        Mipix2_endo_hemo_n_ff_p{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_p == j);
        
        Mipix2_hemo_n_ff_n{j,i} = img2(Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_n == j);
        Mipix2_epi_hemo_n_ff_n{j,i} = img2(Mask_Segn2_epi .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_n == j);
        Mipix2_endo_hemo_n_ff_n{j,i} = img2(Mask_Segn2_endo .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_n == j);
        
        Mipix3_hemo_p{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_p == j);
        Mipix3_epi_hemo_p{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_p == j);
        Mipix3_endo_hemo_p{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_p == j);
        
        Mipix3_hemo_n{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n == j);
        Mipix3_epi_hemo_n{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n == j);
        Mipix3_endo_hemo_n{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n == j);
        
        Mipix3_hemo_n_ff_p{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_p == j);
        Mipix3_epi_hemo_n_ff_p{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_p == j);
        Mipix3_endo_hemo_n_ff_p{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_p == j);
        
        Mipix3_hemo_n_ff_n{j,i} = img3(Mask_Segn3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_n == j);
        Mipix3_epi_hemo_n_ff_n{j,i} = img3(Mask_Segn3_epi .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_n == j);
        Mipix3_endo_hemo_n_ff_n{j,i} = img3(Mask_Segn3_endo .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_n == j);
    end
    
    Mi_idx1_hemo_p = [];
    Mi_idx1_epi_hemo_p = [];
    Mi_idx1_endo_hemo_p = [];
    Mi_idx1_hemo_n = [];
    Mi_idx1_epi_hemo_n = [];
    Mi_idx1_endo_hemo_n = [];
    Mi_idx1_hemo_n_ff_p = [];
    Mi_idx1_epi_hemo_n_ff_p = [];
    Mi_idx1_endo_hemo_n_ff_p = [];
    Mi_idx1_hemo_n_ff_n = [];
    Mi_idx1_epi_hemo_n_ff_n = [];
    Mi_idx1_endo_hemo_n_ff_n = [];
    for j = 1:Segn
        if ~isempty(Mipix_hemo_p{j,i})
            Mi_idx1_hemo_p = [Mi_idx1_hemo_p, j];
            Mi_idx1_hemo_n = [Mi_idx1_hemo_n, j];
            Mi_idx1_hemo_n_ff_p = [Mi_idx1_hemo_n_ff_p, j];
            Mi_idx1_hemo_n_ff_n = [Mi_idx1_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix_epi_hemo_p{j,i})
            Mi_idx1_epi_hemo_p = [Mi_idx1_epi_hemo_p, j];
            Mi_idx1_epi_hemo_n = [Mi_idx1_epi_hemo_n, j];
            Mi_idx1_epi_hemo_n_ff_p = [Mi_idx1_epi_hemo_n_ff_p, j];
            Mi_idx1_epi_hemo_n_ff_n = [Mi_idx1_epi_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix_endo_hemo_p{j,i})
            Mi_idx1_endo_hemo_p = [Mi_idx1_endo_hemo_p, j];
            Mi_idx1_endo_hemo_n = [Mi_idx1_endo_hemo_n, j];
            Mi_idx1_endo_hemo_n_ff_p = [Mi_idx1_endo_hemo_n_ff_p, j];
            Mi_idx1_endo_hemo_n_ff_n = [Mi_idx1_endo_hemo_n_ff_n, j];
        end
    end
    
    Mi_idx2_hemo_p = [];
    Mi_idx2_epi_hemo_p = [];
    Mi_idx2_endo_hemo_p = [];
    Mi_idx2_hemo_n = [];
    Mi_idx2_epi_hemo_n = [];
    Mi_idx2_endo_hemo_n = [];
    Mi_idx2_hemo_n_ff_p = [];
    Mi_idx2_epi_hemo_n_ff_p = [];
    Mi_idx2_endo_hemo_n_ff_p = [];
    Mi_idx2_hemo_n_ff_n = [];
    Mi_idx2_epi_hemo_n_ff_n = [];
    Mi_idx2_endo_hemo_n_ff_n = [];
    for j = 1:Segn
        if ~isempty(Mipix2_hemo_p{j,i})
            Mi_idx2_hemo_p = [Mi_idx2_hemo_p, j];
            Mi_idx2_hemo_n = [Mi_idx2_hemo_n, j];
            Mi_idx2_hemo_n_ff_p = [Mi_idx2_hemo_n_ff_p, j];
            Mi_idx2_hemo_n_ff_n = [Mi_idx2_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix2_epi_hemo_p{j,i})
            Mi_idx2_epi_hemo_p = [Mi_idx2_epi_hemo_p, j];
            Mi_idx2_epi_hemo_n = [Mi_idx2_epi_hemo_n, j];
            Mi_idx2_epi_hemo_n_ff_p = [Mi_idx2_epi_hemo_n_ff_p, j];
            Mi_idx2_epi_hemo_n_ff_n = [Mi_idx2_epi_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix2_endo_hemo_p{j,i})
            Mi_idx2_endo_hemo_p = [Mi_idx2_endo_hemo_p, j];
            Mi_idx2_endo_hemo_n = [Mi_idx2_endo_hemo_n, j];
            Mi_idx2_endo_hemo_n_ff_p = [Mi_idx2_endo_hemo_n_ff_p, j];
            Mi_idx2_endo_hemo_n_ff_n = [Mi_idx2_endo_hemo_n_ff_n, j];
        end
    end
    
    Mi_idx3_hemo_p = [];
    Mi_idx3_epi_hemo_p = [];
    Mi_idx3_endo_hemo_p = [];
    Mi_idx3_hemo_n = [];
    Mi_idx3_epi_hemo_n = [];
    Mi_idx3_endo_hemo_n = [];
    Mi_idx3_hemo_n_ff_p = [];
    Mi_idx3_epi_hemo_n_ff_p = [];
    Mi_idx3_endo_hemo_n_ff_p = [];
    Mi_idx3_hemo_n_ff_n = [];
    Mi_idx3_epi_hemo_n_ff_n = [];
    Mi_idx3_endo_hemo_n_ff_n = [];
    for j = 1:Segn
        if ~isempty(Mipix3_hemo_p{j,i})
            Mi_idx3_hemo_p = [Mi_idx3_hemo_p, j];
            Mi_idx3_hemo_n = [Mi_idx3_hemo_n, j];
            Mi_idx3_hemo_n_ff_p = [Mi_idx3_hemo_n_ff_p, j];
            Mi_idx3_hemo_n_ff_n = [Mi_idx3_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix3_epi_hemo_p{j,i})
            Mi_idx3_epi_hemo_p = [Mi_idx3_epi_hemo_p, j];
            Mi_idx3_epi_hemo_n = [Mi_idx3_epi_hemo_n, j];
            Mi_idx3_epi_hemo_n_ff_p = [Mi_idx3_epi_hemo_n_ff_p, j];
            Mi_idx3_epi_hemo_n_ff_n = [Mi_idx3_epi_hemo_n_ff_n, j];
        end
        if ~isempty(Mipix3_endo_hemo_p{j,i})
            Mi_idx3_endo_hemo_p = [Mi_idx3_endo_hemo_p, j];
            Mi_idx3_endo_hemo_n = [Mi_idx3_endo_hemo_n, j];
            Mi_idx3_endo_hemo_n_ff_p = [Mi_idx3_endo_hemo_n_ff_p, j];
            Mi_idx3_endo_hemo_n_ff_n = [Mi_idx3_endo_hemo_n_ff_n, j];
        end
    end
    
    % flatten and get the mean 
    Mi_idx_intersect_hemo_p = intersect(Mi_idx1_hemo_p, Mi_idx2_hemo_p);
    Mi_idx_intersect_epi_hemo_p = intersect(Mi_idx1_epi_hemo_p, Mi_idx2_epi_hemo_p);
    Mi_idx_intersect_endo_hemo_p = intersect(Mi_idx1_endo_hemo_p, Mi_idx2_endo_hemo_p);
    Mi_idx_intersect_hemo_n = intersect(Mi_idx1_hemo_n, Mi_idx2_hemo_n);
    Mi_idx_intersect_epi_hemo_n = intersect(Mi_idx1_epi_hemo_n, Mi_idx2_epi_hemo_n);
    Mi_idx_intersect_endo_hemo_n = intersect(Mi_idx1_endo_hemo_n, Mi_idx2_endo_hemo_n);
    Mi_idx_intersect_hemo_n_ff_p = intersect(Mi_idx1_hemo_n_ff_p, Mi_idx2_hemo_n_ff_p);
    Mi_idx_intersect_epi_hemo_n_ff_p = intersect(Mi_idx1_epi_hemo_n_ff_p, Mi_idx2_epi_hemo_n_ff_p);
    Mi_idx_intersect_endo_hemo_n_ff_p = intersect(Mi_idx1_endo_hemo_n_ff_p, Mi_idx2_endo_hemo_n_ff_p);
    Mi_idx_intersect_hemo_n_ff_n = intersect(Mi_idx1_hemo_n_ff_n, Mi_idx2_hemo_n_ff_n);
    Mi_idx_intersect_epi_hemo_n_ff_n = intersect(Mi_idx1_epi_hemo_n_ff_n, Mi_idx2_epi_hemo_n_ff_n);
    Mi_idx_intersect_endo_hemo_n_ff_n = intersect(Mi_idx1_endo_hemo_n_ff_n, Mi_idx2_endo_hemo_n_ff_n);
    
    % Total
    [Mipix_flat_hemo_p, Mipix_flat2_hemo_p, Mipix_flat3_hemo_p, ...
        Mipix_mean_hemo_p, Mipix_mean2_hemo_p, Mipix_mean3_hemo_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_hemo_p, Mipix_hemo_p, Mipix2_hemo_p, Mipix3_hemo_p,...
    Mipix_mean_hemo_p, Mipix_mean2_hemo_p, Mipix_mean3_hemo_p);

    [Mipix_flat_hemo_n, Mipix_flat2_hemo_n, Mipix_flat3_hemo_n, ...
        Mipix_mean_hemo_n, Mipix_mean2_hemo_n, Mipix_mean3_hemo_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_hemo_n, Mipix_hemo_n, Mipix2_hemo_n, Mipix3_hemo_n,...
    Mipix_mean_hemo_n, Mipix_mean2_hemo_n, Mipix_mean3_hemo_n);

    [Mipix_flat_hemo_n_ff_p, Mipix_flat2_hemo_n_ff_p, Mipix_flat3_hemo_n_ff_p, ...
        Mipix_mean_hemo_n_ff_p, Mipix_mean2_hemo_n_ff_p, Mipix_mean3_hemo_n_ff_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_hemo_n_ff_p, Mipix_hemo_n_ff_p, Mipix2_hemo_n_ff_p, Mipix3_hemo_n_ff_p,...
    Mipix_mean_hemo_n_ff_p, Mipix_mean2_hemo_n_ff_p, Mipix_mean3_hemo_n_ff_p);

    [Mipix_flat_hemo_n_ff_n, Mipix_flat2_hemo_n_ff_n, Mipix_flat3_hemo_n_ff_n, ...
        Mipix_mean_hemo_n_ff_n, Mipix_mean2_hemo_n_ff_n, Mipix_mean3_hemo_n_ff_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_hemo_n_ff_n, Mipix_hemo_n_ff_n, Mipix2_hemo_n_ff_n, Mipix3_hemo_n_ff_n,...
    Mipix_mean_hemo_n_ff_n, Mipix_mean2_hemo_n_ff_n, Mipix_mean3_hemo_n_ff_n);
    
    % Endo
    [Mipix_flat_endo_hemo_p, Mipix_flat2_endo_hemo_p, Mipix_flat3_endo_hemo_p, ...
        Mipix_mean_endo_hemo_p, Mipix_mean2_endo_hemo_p, Mipix_mean3_endo_hemo_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_endo_hemo_p, Mipix_endo_hemo_p, Mipix2_endo_hemo_p, Mipix3_endo_hemo_p,...
    Mipix_mean_endo_hemo_p, Mipix_mean2_endo_hemo_p, Mipix_mean3_endo_hemo_p);

    [Mipix_flat_endo_hemo_n, Mipix_flat2_endo_hemo_n, Mipix_flat3_endo_hemo_n, ...
        Mipix_mean_endo_hemo_n, Mipix_mean2_endo_hemo_n, Mipix_mean3_endo_hemo_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_endo_hemo_n, Mipix_endo_hemo_n, Mipix2_endo_hemo_n, Mipix3_endo_hemo_n,...
    Mipix_mean_endo_hemo_n, Mipix_mean2_endo_hemo_n, Mipix_mean3_endo_hemo_n);

    [Mipix_flat_endo_hemo_n_ff_p, Mipix_flat2_endo_hemo_n_ff_p, Mipix_flat3_endo_hemo_n_ff_p, ...
        Mipix_mean_endo_hemo_n_ff_p, Mipix_mean2_endo_hemo_n_ff_p, Mipix_mean3_endo_hemo_n_ff_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_endo_hemo_n_ff_p, Mipix_endo_hemo_n_ff_p, Mipix2_endo_hemo_n_ff_p, Mipix3_endo_hemo_n_ff_p,...
    Mipix_mean_endo_hemo_n_ff_p, Mipix_mean2_endo_hemo_n_ff_p, Mipix_mean3_endo_hemo_n_ff_p);

    [Mipix_flat_endo_hemo_n_ff_n, Mipix_flat2_endo_hemo_n_ff_n, Mipix_flat3_endo_hemo_n_ff_n, ...
        Mipix_mean_endo_hemo_n_ff_n, Mipix_mean2_endo_hemo_n_ff_n, Mipix_mean3_endo_hemo_n_ff_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_endo_hemo_n_ff_n, Mipix_endo_hemo_n_ff_n, Mipix2_endo_hemo_n_ff_n, Mipix3_endo_hemo_n_ff_n,...
    Mipix_mean_endo_hemo_n_ff_n, Mipix_mean2_endo_hemo_n_ff_n, Mipix_mean3_endo_hemo_n_ff_n);
    

    % Epi
    [Mipix_flat_epi_hemo_p, Mipix_flat2_epi_hemo_p, Mipix_flat3_epi_hemo_p, ...
        Mipix_mean_epi_hemo_p, Mipix_mean2_epi_hemo_p, Mipix_mean3_epi_hemo_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_epi_hemo_p, Mipix_epi_hemo_p, Mipix2_epi_hemo_p, Mipix3_epi_hemo_p,...
    Mipix_mean_epi_hemo_p, Mipix_mean2_epi_hemo_p, Mipix_mean3_epi_hemo_p);

    [Mipix_flat_epi_hemo_n, Mipix_flat2_epi_hemo_n, Mipix_flat3_epi_hemo_n, ...
        Mipix_mean_epi_hemo_n, Mipix_mean2_epi_hemo_n, Mipix_mean3_epi_hemo_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_epi_hemo_n, Mipix_epi_hemo_n, Mipix2_epi_hemo_n, Mipix3_epi_hemo_n,...
    Mipix_mean_epi_hemo_n, Mipix_mean2_epi_hemo_n, Mipix_mean3_epi_hemo_n);

    [Mipix_flat_epi_hemo_n_ff_p, Mipix_flat2_epi_hemo_n_ff_p, Mipix_flat3_epi_hemo_n_ff_p, ...
        Mipix_mean_epi_hemo_n_ff_p, Mipix_mean2_epi_hemo_n_ff_p, Mipix_mean3_epi_hemo_n_ff_p] = ...
    Func_flattenNmean(i, Mi_idx_intersect_epi_hemo_n_ff_p, Mipix_epi_hemo_n_ff_p, Mipix2_epi_hemo_n_ff_p, Mipix3_epi_hemo_n_ff_p,...
    Mipix_mean_epi_hemo_n_ff_p, Mipix_mean2_epi_hemo_n_ff_p, Mipix_mean3_epi_hemo_n_ff_p);

    [Mipix_flat_epi_hemo_n_ff_n, Mipix_flat2_epi_hemo_n_ff_n, Mipix_flat3_epi_hemo_n_ff_n, ...
        Mipix_mean_epi_hemo_n_ff_n, Mipix_mean2_epi_hemo_n_ff_n, Mipix_mean3_epi_hemo_n_ff_n] = ...
    Func_flattenNmean(i, Mi_idx_intersect_epi_hemo_n_ff_n, Mipix_epi_hemo_n_ff_n, Mipix2_epi_hemo_n_ff_n, Mipix3_epi_hemo_n_ff_n,...
    Mipix_mean_epi_hemo_n_ff_n, Mipix_mean2_epi_hemo_n_ff_n, Mipix_mean3_epi_hemo_n_ff_n);
    
    
    % Save MI_Chord_Analysis2
    MI_Chord_Analysis2(i).Mask_Segn = Mask_Segn; 
    MI_Chord_Analysis2(i).Mask_Segn2 = Mask_Segn2; 
    MI_Chord_Analysis2(i).Mask_Segn3 = Mask_Segn3;
    
    MI_Chord_Analysis2(i).Mipix_hemo_p = Mipix_hemo_p; 
    MI_Chord_Analysis2(i).Mipix2_hemo_p = Mipix2_hemo_p; 
    MI_Chord_Analysis2(i).Mipix3_hemo_p = Mipix3_hemo_p;
    MI_Chord_Analysis2(i).Mipix_mean_hemo_p = Mipix_mean_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_hemo_p = Mipix_mean2_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_hemo_p = Mipix_mean3_hemo_p;

    MI_Chord_Analysis2(i).Mipix_hemo_n = Mipix_hemo_n;
    MI_Chord_Analysis2(i).Mipix2_hemo_n = Mipix2_hemo_n;
    MI_Chord_Analysis2(i).Mipix3_hemo_n = Mipix3_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean_hemo_n = Mipix_mean_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean2_hemo_n = Mipix_mean2_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean3_hemo_n = Mipix_mean3_hemo_n;
    
    MI_Chord_Analysis2(i).Mipix_hemo_n_ff_p = Mipix_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix2_hemo_n_ff_p = Mipix2_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix3_hemo_n_ff_p = Mipix3_hemo_n_ff_p;
    MI_Chord_Analysis2(i).Mipix_mean_hemo_n_ff_p = Mipix_mean_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_hemo_n_ff_p = Mipix_mean2_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_hemo_n_ff_p = Mipix_mean3_hemo_n_ff_p;
    
    MI_Chord_Analysis2(i).Mipix_hemo_n_ff_n = Mipix_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix2_hemo_n_ff_n = Mipix2_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix3_hemo_n_ff_n = Mipix3_hemo_n_ff_n;
    MI_Chord_Analysis2(i).Mipix_mean_hemo_n_ff_n = Mipix_mean_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean2_hemo_n_ff_n = Mipix_mean2_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean3_hemo_n_ff_n = Mipix_mean3_hemo_n_ff_n;
    
    % Epi
    MI_Chord_Analysis2(i).Mask_Segn_epi = Mask_Segn_epi; 
    MI_Chord_Analysis2(i).Mask_Segn2_epi = Mask_Segn2_epi; 
    MI_Chord_Analysis2(i).Mask_Segn3_epi = Mask_Segn3_epi;
    
    MI_Chord_Analysis2(i).Mipix_epi_hemo_p = Mipix_epi_hemo_p; 
    MI_Chord_Analysis2(i).Mipix2_epi_hemo_p = Mipix2_epi_hemo_p; 
    MI_Chord_Analysis2(i).Mipix3_epi_hemo_p = Mipix3_epi_hemo_p;
    MI_Chord_Analysis2(i).Mipix_mean_epi_hemo_p = Mipix_mean_epi_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_epi_hemo_p = Mipix_mean2_epi_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_epi_hemo_p = Mipix_mean3_epi_hemo_p;

    MI_Chord_Analysis2(i).Mipix_epi_hemo_n = Mipix_epi_hemo_n;
    MI_Chord_Analysis2(i).Mipix2_epi_hemo_n = Mipix2_epi_hemo_n;
    MI_Chord_Analysis2(i).Mipix3_epi_hemo_n = Mipix3_epi_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean_epi_hemo_n = Mipix_mean_epi_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean2_epi_hemo_n = Mipix_mean2_epi_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean3_epi_hemo_n = Mipix_mean3_epi_hemo_n;
    
    MI_Chord_Analysis2(i).Mipix_epi_hemo_n_ff_p = Mipix_epi_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix2_epi_hemo_n_ff_p = Mipix2_epi_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix3_epi_hemo_n_ff_p = Mipix3_epi_hemo_n_ff_p;
    MI_Chord_Analysis2(i).Mipix_mean_epi_hemo_n_ff_p = Mipix_mean_epi_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_epi_hemo_n_ff_p = Mipix_mean2_epi_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_epi_hemo_n_ff_p = Mipix_mean3_epi_hemo_n_ff_p;
    
    MI_Chord_Analysis2(i).Mipix_epi_hemo_n_ff_n = Mipix_epi_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix2_epi_hemo_n_ff_n = Mipix2_epi_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix3_epi_hemo_n_ff_n = Mipix3_epi_hemo_n_ff_n;
    MI_Chord_Analysis2(i).Mipix_mean_epi_hemo_n_ff_n = Mipix_mean_epi_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean2_epi_hemo_n_ff_n = Mipix_mean2_epi_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean3_epi_hemo_n_ff_n = Mipix_mean3_epi_hemo_n_ff_n;
    
    
    % Endo
    MI_Chord_Analysis2(i).Mask_Segn_endo = Mask_Segn_endo; 
    MI_Chord_Analysis2(i).Mask_Segn2_endo = Mask_Segn2_endo; 
    MI_Chord_Analysis2(i).Mask_Segn3_endo = Mask_Segn3_endo;
    
    MI_Chord_Analysis2(i).Mipix_endo_hemo_p = Mipix_endo_hemo_p; 
    MI_Chord_Analysis2(i).Mipix2_endo_hemo_p = Mipix2_endo_hemo_p; 
    MI_Chord_Analysis2(i).Mipix3_endo_hemo_p = Mipix3_endo_hemo_p;
    MI_Chord_Analysis2(i).Mipix_mean_endo_hemo_p = Mipix_mean_endo_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_endo_hemo_p = Mipix_mean2_endo_hemo_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_endo_hemo_p = Mipix_mean3_endo_hemo_p;

    MI_Chord_Analysis2(i).Mipix_endo_hemo_n = Mipix_endo_hemo_n;
    MI_Chord_Analysis2(i).Mipix2_endo_hemo_n = Mipix2_endo_hemo_n;
    MI_Chord_Analysis2(i).Mipix3_endo_hemo_n = Mipix3_endo_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean_endo_hemo_n = Mipix_mean_endo_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean2_endo_hemo_n = Mipix_mean2_endo_hemo_n;
    MI_Chord_Analysis2(i).Mipix_mean3_endo_hemo_n = Mipix_mean3_endo_hemo_n;
    
    MI_Chord_Analysis2(i).Mipix_endo_hemo_n_ff_p = Mipix_endo_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix2_endo_hemo_n_ff_p = Mipix2_endo_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix3_endo_hemo_n_ff_p = Mipix3_endo_hemo_n_ff_p;
    MI_Chord_Analysis2(i).Mipix_mean_endo_hemo_n_ff_p = Mipix_mean_endo_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean2_endo_hemo_n_ff_p = Mipix_mean2_endo_hemo_n_ff_p; 
    MI_Chord_Analysis2(i).Mipix_mean3_endo_hemo_n_ff_p = Mipix_mean3_endo_hemo_n_ff_p;
    
    MI_Chord_Analysis2(i).Mipix_endo_hemo_n_ff_n = Mipix_endo_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix2_endo_hemo_n_ff_n = Mipix2_endo_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix3_endo_hemo_n_ff_n = Mipix3_endo_hemo_n_ff_n;
    MI_Chord_Analysis2(i).Mipix_mean_endo_hemo_n_ff_n = Mipix_mean_endo_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean2_endo_hemo_n_ff_n = Mipix_mean2_endo_hemo_n_ff_n; 
    MI_Chord_Analysis2(i).Mipix_mean3_endo_hemo_n_ff_n = Mipix_mean3_endo_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).fixed_eroded = fixed_eroded;
    MI_Chord_Analysis2(i).movingRegistered_eroded = movingRegistered_eroded;
    
    % 
    t1_mean_hemo_p = mean(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_p));
    t1_mean_hemo_n = mean(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n));
    t1_mean_hemo_n_ff_p = mean(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_p));    
    t1_mean_hemo_n_ff_n = mean(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_n));
    t1_sd_hemo_p = std(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_p));
    t1_sd_hemo_n = std(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n));
    t1_sd_hemo_n_ff_p = std(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_p));    
    t1_sd_hemo_n_ff_n = std(nonzeros(img .* roi_in_myo_t1(:,:,i) .* fixed_eroded .* hemo_n_ff_n));
    
    img2(img2 > 100) = 100;
    img2(img2 < 0) = 0;
    ff_mean_hemo_p = mean(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_p));
    ff_mean_hemo_n = mean(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n));
    ff_mean_hemo_n_ff_p = mean(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_p));
    ff_mean_hemo_n_ff_n = mean(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_n));
    ff_sd_hemo_p = std(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_p));
    ff_sd_hemo_n = std(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n));
    ff_sd_hemo_n_ff_p = std(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_p));
    ff_sd_hemo_n_ff_n = std(nonzeros(img2 .* movingRegistered_roi_ff .* movingRegistered_eroded .* hemo_n_ff_n));
    
    img3(img3 < 0) = 0;
    r2star_mean_hemo_p = mean(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_p));
    r2star_mean_hemo_n = mean(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n));
    r2star_mean_hemo_n_ff_p = mean(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_p));
    r2star_mean_hemo_n_ff_n = mean(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_n));
    r2star_sd_hemo_p = std(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_p));
    r2star_sd_hemo_n = std(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n));
    r2star_sd_hemo_n_ff_p = std(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_p));
    r2star_sd_hemo_n_ff_n = std(nonzeros(img3 .* movingRegistered_roi_r2star .* movingRegistered_eroded .* hemo_n_ff_n));
    
    MI_Chord_Analysis2(i).t1_mean_hemo_p = t1_mean_hemo_p;
    MI_Chord_Analysis2(i).t1_mean_hemo_n = t1_mean_hemo_n;
    MI_Chord_Analysis2(i).t1_mean_hemo_n_ff_p = t1_mean_hemo_n_ff_p;
    MI_Chord_Analysis2(i).t1_mean_hemo_n_ff_n = t1_mean_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).ff_mean_hemo_p = ff_mean_hemo_p;
    MI_Chord_Analysis2(i).ff_mean_hemo_n = ff_mean_hemo_n;
    MI_Chord_Analysis2(i).ff_mean_hemo_n_ff_p = ff_mean_hemo_n_ff_p;
    MI_Chord_Analysis2(i).ff_mean_hemo_n_ff_n = ff_mean_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).r2star_mean_hemo_p = r2star_mean_hemo_p;
    MI_Chord_Analysis2(i).r2star_mean_hemo_n = r2star_mean_hemo_n;
    MI_Chord_Analysis2(i).r2star_mean_hemo_n_ff_p = r2star_mean_hemo_n_ff_p;
    MI_Chord_Analysis2(i).r2star_mean_hemo_n_ff_n = r2star_mean_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).t1_sd_hemo_p = t1_sd_hemo_p;
    MI_Chord_Analysis2(i).t1_sd_hemo_n = t1_sd_hemo_n;
    MI_Chord_Analysis2(i).t1_sd_hemo_n_ff_p = t1_sd_hemo_n_ff_p;
    MI_Chord_Analysis2(i).t1_sd_hemo_n_ff_n = t1_sd_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).ff_sd_hemo_p = ff_sd_hemo_p;
    MI_Chord_Analysis2(i).ff_sd_hemo_n = ff_sd_hemo_n;
    MI_Chord_Analysis2(i).ff_sd_hemo_n_ff_p = ff_sd_hemo_n_ff_p;
    MI_Chord_Analysis2(i).ff_sd_hemo_n_ff_n = ff_sd_hemo_n_ff_n;
    
    MI_Chord_Analysis2(i).r2star_sd_hemo_p = r2star_sd_hemo_p;
    MI_Chord_Analysis2(i).r2star_sd_hemo_n = r2star_sd_hemo_n;
    MI_Chord_Analysis2(i).r2star_sd_hemo_n_ff_p = r2star_sd_hemo_n_ff_p;
    MI_Chord_Analysis2(i).r2star_sd_hemo_n_ff_n = r2star_sd_hemo_n_ff_n;
end
    save(MI_Chord_Analysis_fname, 'MI_Chord_Analysis2');
end