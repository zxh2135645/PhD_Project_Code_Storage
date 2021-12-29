function Func_T1FP_Chord_ReAnalysis_EndoEpi(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,LR_mdl_fname,chord_values_fname,chord_values_fname2)

t1_cell = cell(50,size(roi_in_myo_t1, 3));
ff_cell = cell(50,size(roi_in_myo_t1, 3));
mean_t1_array = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array = zeros(50,size(roi_in_myo_t1, 3));
mean_t1_hemo_array = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_hemo_array = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_hemo_array = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_t1_hemo_array = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_hemo_array = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_hemo_array = 1000*ones(50,size(roi_in_myo_t1, 3));

mean_t1_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_t1_hemo_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_hemo_array_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_hemo_array_endo = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_t1_hemo_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_hemo_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_hemo_array_endo = 1000*ones(50,size(roi_in_myo_t1, 3));

mean_t1_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_t1_hemo_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_hemo_array_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_hemo_array_epi = zeros(50,size(roi_in_myo_t1, 3));

sd_t1_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_t1_hemo_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_hemo_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_hemo_array_epi = 1000*ones(50,size(roi_in_myo_t1, 3));

mean_t1_array_remote = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_remote = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_remote = zeros(50,size(roi_in_myo_t1, 3));
sd_t1_array_remote = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array_remote = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array_remote = 1000*ones(50,size(roi_in_myo_t1, 3));

mean_t1_array_remote_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_remote_endo = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_remote_endo = zeros(50,size(roi_in_myo_t1, 3));
sd_t1_array_remote_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array_remote_endo = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array_remote_endo = 1000*ones(50,size(roi_in_myo_t1, 3));

mean_t1_array_remote_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_ff_array_remote_epi = zeros(50,size(roi_in_myo_t1, 3));
mean_r2star_array_remote_epi = zeros(50,size(roi_in_myo_t1, 3));
sd_t1_array_remote_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_ff_array_remote_epi = 1000*ones(50,size(roi_in_myo_t1, 3));
sd_r2star_array_remote_epi = 1000*ones(50,size(roi_in_myo_t1, 3));

for i = 1:size(roi_in_myo_t1, 3)
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    img3 = r2star(:,:,i);
    
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    
    
    img2(img2 > 100) = 100;
    img2(img2 < 0) = 0;
    
    figure('Position', [100 0 400 800]);
    subplot(3,2,1);
    imagesc(fixed); axis image; title('T1 (fixed)'); colormap(brewermap([],'*RdYlBu'));axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,2);
    imagesc(moving); axis image; title('FF (moving)');axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    se = strel('disk', 1);
    
    I1 = moving; I2 = fixed;
    % Set static and moving image
    S=I2; M=I1;
    
    % resizepercentag
    [movingRegistered,Bx,By,Fx,Fy] = register_images(M,S);
    movingRegistered = movingRegistered > 0.5;
    
    if (size(img, 1) ~= size(img2, 1)) || (size(img, 2) ~= size(img2, 2))
        img2 = imresize(img2,size(img),'bicubic');
        img3 = imresize(img3,size(img),'bicubic');
        myo_ff_temp = imresize(myo_ff(:,:,i),size(img),'bicubic');
        roi_in_myo_ff_temp = imresize(roi_in_myo_ff(:,:,i),size(img),'bicubic');
        remote_in_myo_ff_temp = imresize(remote_in_myo_ff(:,:,i),size(img),'bicubic');
    else
        myo_ff_temp = myo_ff(:,:,i);
        roi_in_myo_ff_temp = roi_in_myo_ff(:,:,i);
        remote_in_myo_ff_temp = remote_in_myo_ff(:,:,i);
    end
    
    img2 = movepixels(img2,Bx,By);
    img3 = movepixels(img3,Bx,By);
    
    subplot(3,2,3); imagesc(img); axis image; title('T1 map'); axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,4); imagesc(img2); axis image; caxis([0 50]); title('FF map'); axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    movingRegistered_myo_ff = movepixels(myo_ff_temp,Bx,By)>0.5;
    movingRegistered_roi_ff = movepixels(roi_in_myo_ff_temp,Bx,By)>0.5;
    movingRegistered_roi_r2star = movepixels(roi_in_myo_ff_temp,Bx,By)>0.5;
    movingRegistered_remote_ff = movepixels(remote_in_myo_ff_temp,Bx,By)>0.5;
    
    univ_roi = roi_in_myo_t1(:,:,i) & movingRegistered_roi_ff;
    
    subplot(3,2,5);
    imshowpair(fixed,movingRegistered,'Scaling','joint'); title('Registered');
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,6); imagesc(double(movingRegistered_myo_ff&fixed) + double(univ_roi) + 2*double(movingRegistered_remote_ff)); axis image;
    title(cat(2, 'Slice = ', num2str(i)));
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    saveas(gcf, cat(2, tp_dir2, 'MyocardiumRegistration_demon_Slice', num2str(i), '.png'));
    
    univ_myo = movingRegistered_myo_ff&fixed;
    fixed_eroded = imerode(myo_t1(:,:,i), se);
    BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
    center_fixed = imfill(BW_skel, 'hole');
    center_fixed = imopen(center_fixed, se); % Removing spikes
    fixedRegistered_epi = fixed_eroded - center_fixed > 0;
    fixedRegistered_endo = center_fixed + fixed_eroded > 1;
    
    movingRegistered_eroded = imerode(movingRegistered_myo_ff, se);
    BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
    center_moving = imfill(BW_skel, 'hole');
    center_moving = imopen(center_moving, se); % Removing spikes
    movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
    movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
    
    univ_myo_eroded = imerode(univ_myo, se);
    
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(img, univ_myo_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, univ_myo_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, univ_myo_eroded, Segn, Groove);
    
    [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(img, fixedRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(img, fixedRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_epi] = AHASegmentation(img2, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_endo] = AHASegmentation(img2, movingRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_epi] = AHASegmentation(img3, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_endo] = AHASegmentation(img3, movingRegistered_endo, Segn, Groove);
    
    seg_mask_overlap = zeros(size(img));
    
    for j = 1:Segn
        
        Mipix{j,i} = img(Mask_Segn .* univ_roi .* fixed_eroded == j);
        Mipix_epi{j,i} = img(Mask_Segn_epi .* univ_roi .* fixed_eroded == j);
        Mipix_endo{j,i} = img(Mask_Segn_endo .* univ_roi .* fixed_eroded == j);
        
        Mipix2{j,i} = img2(Mask_Segn2 .* univ_roi .* movingRegistered_eroded == j);
        Mipix2_epi{j,i} = img2(Mask_Segn2_epi .* univ_roi .* movingRegistered_eroded == j);
        Mipix2_endo{j,i} = img2(Mask_Segn2_endo .* univ_roi .* movingRegistered_eroded == j);
        
        Mipix3{j,i} = img3(Mask_Segn3 .* univ_roi .* movingRegistered_eroded == j);
        Mipix3_epi{j,i} = img3(Mask_Segn3_epi .* univ_roi .* movingRegistered_eroded == j);
        Mipix3_endo{j,i} = img3(Mask_Segn3_endo .* univ_roi .* movingRegistered_eroded == j);
        
        remote_mean_t1 = mean(nonzeros(remote_in_myo_t1(:,:,i) .* img));
        remote_sd_t1 = std(nonzeros(remote_in_myo_t1(:,:,i) .* img));
        thresh = remote_mean_t1 - 2*remote_sd_t1;
        hemo_mask = (img<thresh).*univ_roi;
        
        % For debugging
        seg_mask_t1 = Mask_Segn .* (univ_roi-hemo_mask) .* fixed_eroded == j;
        seg_mask_ff = Mask_Segn2 .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        seg_mask_r2star = Mask_Segn3 .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        
        seg_mask_t1_endo = Mask_Segn_endo .* (univ_roi-hemo_mask) .* fixed_eroded == j;
        seg_mask_t1_epi = Mask_Segn_epi .* (univ_roi-hemo_mask) .* fixed_eroded == j;
        
        seg_mask_ff_endo = Mask_Segn2_endo .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        seg_mask_ff_epi = Mask_Segn2_epi .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        
        seg_mask_r2star_endo = Mask_Segn3_endo .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        seg_mask_r2star_epi = Mask_Segn3_epi .* (univ_roi-hemo_mask) .* movingRegistered_eroded == j;
        
        seg_mask_t1_remote = Mask_Segn .* movingRegistered_remote_ff .* fixed_eroded == j;
        seg_mask_ff_remote = Mask_Segn2 .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        seg_mask_r2star_remote = Mask_Segn3 .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        
        hemo_mask_t1 = Mask_Segn .* hemo_mask .* fixed_eroded == j;
        hemo_mask_ff = Mask_Segn2 .* hemo_mask .* movingRegistered_eroded == j;
        hemo_mask_r2star = Mask_Segn3 .* hemo_mask .* movingRegistered_eroded == j;
        
        seg_mask_t1_remote_endo = Mask_Segn_endo .* movingRegistered_remote_ff .* fixed_eroded == j;
        seg_mask_ff_remote_endo = Mask_Segn2_endo .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        seg_mask_r2star_remote_endo = Mask_Segn3_endo .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        
        seg_mask_t1_remote_epi = Mask_Segn_epi .* movingRegistered_remote_ff .* fixed_eroded == j;
        seg_mask_ff_remote_epi = Mask_Segn2_epi .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        seg_mask_r2star_remote_epi = Mask_Segn3_epi .* movingRegistered_remote_ff .* movingRegistered_eroded == j;
        
        hemo_mask_t1_endo = Mask_Segn_endo .* hemo_mask .* fixed_eroded == j;
        hemo_mask_ff_endo = Mask_Segn2_endo .* hemo_mask .* movingRegistered_eroded == j;
        hemo_mask_r2star_endo = Mask_Segn3_endo .* hemo_mask .* movingRegistered_eroded == j;
        
        hemo_mask_t1_epi = Mask_Segn_epi .* hemo_mask .* fixed_eroded == j;
        hemo_mask_ff_epi = Mask_Segn2_epi .* hemo_mask .* movingRegistered_eroded == j;
        hemo_mask_r2star_epi = Mask_Segn3_epi .* hemo_mask .* movingRegistered_eroded == j;
        
        %
        if any(seg_mask_t1(:))
            t1_cell{j,i} = nonzeros(img.*seg_mask_t1);
            ff_cell{j,i} = nonzeros(img2.*seg_mask_ff);
            mean_t1_array(j,i) = mean(nonzeros(img.*seg_mask_t1));
            mean_ff_array(j,i) = mean(nonzeros(img2.*seg_mask_ff));
            mean_r2star_array(j,i) = mean(nonzeros(img3.*seg_mask_r2star));
            sd_t1_array(j,i) = std(nonzeros(img.*seg_mask_t1));
            sd_ff_array(j,i) = std(nonzeros(img2.*seg_mask_ff));
            sd_r2star_array(j,i) = std(nonzeros(img3.*seg_mask_r2star));
            
        end
        
        if any(seg_mask_t1_endo(:))
            mean_t1_array_endo(j,i) = mean(nonzeros(img.*seg_mask_t1_endo));
            mean_ff_array_endo(j,i) = mean(nonzeros(img2.*seg_mask_ff_endo));
            mean_r2star_array_endo(j,i) = mean(nonzeros(img3.*seg_mask_r2star_endo));
            
            sd_t1_array_endo(j,i) = std(nonzeros(img.*seg_mask_t1_endo));
            sd_ff_array_endo(j,i) = std(nonzeros(img2.*seg_mask_ff_endo));
            sd_r2star_array_endo(j,i) = std(nonzeros(img3.*seg_mask_r2star_endo));
        end
        
        if any(seg_mask_t1_epi(:))
            mean_t1_array_epi(j,i) = mean(nonzeros(img.*seg_mask_t1_epi));
            mean_ff_array_epi(j,i) = mean(nonzeros(img2.*seg_mask_ff_epi));
            mean_r2star_array_epi(j,i) = mean(nonzeros(img3.*seg_mask_r2star_epi));
            
            sd_t1_array_epi(j,i) = std(nonzeros(img.*seg_mask_t1_epi));
            sd_ff_array_epi(j,i) = std(nonzeros(img2.*seg_mask_ff_epi));
            sd_r2star_array_epi(j,i) = std(nonzeros(img3.*seg_mask_r2star_epi));
        end
        
        if any(seg_mask_t1_remote(:))
            % remote
            mean_t1_array_remote(j,i) = mean(nonzeros(img.*seg_mask_t1_remote));
            mean_ff_array_remote(j,i) = mean(nonzeros(img2.*seg_mask_ff_remote));
            mean_r2star_array_remote(j,i) = mean(nonzeros(img3.*seg_mask_r2star_remote));
            sd_t1_array_remote(j,i) = std(nonzeros(img.*seg_mask_t1_remote));
            sd_ff_array_remote(j,i) = std(nonzeros(img2.*seg_mask_ff_remote));
            sd_r2star_array_remote(j,i) = std(nonzeros(img3.*seg_mask_r2star_remote));
        end
        
        if any(seg_mask_t1_remote_endo(:))
            % remote
            mean_t1_array_remote_endo(j,i) = mean(nonzeros(img.*seg_mask_t1_remote_endo));
            mean_ff_array_remote_endo(j,i) = mean(nonzeros(img2.*seg_mask_ff_remote_endo));
            mean_r2star_array_remote_endo(j,i) = mean(nonzeros(img3.*seg_mask_r2star_remote_endo));
            sd_t1_array_remote_endo(j,i) = std(nonzeros(img.*seg_mask_t1_remote_endo));
            sd_ff_array_remote_endo(j,i) = std(nonzeros(img2.*seg_mask_ff_remote_endo));
            sd_r2star_array_remote_endo(j,i) = std(nonzeros(img3.*seg_mask_r2star_remote_endo));
        end
        
        if any(seg_mask_t1_remote_epi(:))
            % remote
            mean_t1_array_remote_epi(j,i) = mean(nonzeros(img.*seg_mask_t1_remote_epi));
            mean_ff_array_remote_epi(j,i) = mean(nonzeros(img2.*seg_mask_ff_remote_epi));
            mean_r2star_array_remote_epi(j,i) = mean(nonzeros(img3.*seg_mask_r2star_remote_epi));
            sd_t1_array_remote_epi(j,i) = std(nonzeros(img.*seg_mask_t1_remote_epi));
            sd_ff_array_remote_epi(j,i) = std(nonzeros(img2.*seg_mask_ff_remote_epi));
            sd_r2star_array_remote_epi(j,i) = std(nonzeros(img3.*seg_mask_r2star_remote_epi));
        end

        if any(hemo_mask_t1(:))
            mean_t1_hemo_array(j,i) = mean(nonzeros(img .* hemo_mask_t1));
            mean_ff_hemo_array(j,i) = mean(nonzeros(img2 .* hemo_mask_ff));
            mean_r2star_hemo_array(j,i) = mean(nonzeros(img3 .* hemo_mask_r2star));
            sd_t1_hemo_array(j,i) = std(nonzeros(img .* hemo_mask_t1));
            sd_ff_hemo_array(j,i) = std(nonzeros(img2 .* hemo_mask_ff));
            sd_r2star_hemo_array(j,i) = std(nonzeros(img3 .* hemo_mask_r2star));
        end
        
        if any(hemo_mask_t1_endo(:))
            mean_t1_hemo_array_endo(j,i) = mean(nonzeros(img .* hemo_mask_t1_endo));
            mean_ff_hemo_array_endo(j,i) = mean(nonzeros(img2 .* hemo_mask_ff_endo));
            mean_r2star_hemo_array_endo(j,i) = mean(nonzeros(img3 .* hemo_mask_r2star_endo));
            sd_t1_hemo_array_endo(j,i) = std(nonzeros(img .* hemo_mask_t1_endo));
            sd_ff_hemo_array_endo(j,i) = std(nonzeros(img2 .* hemo_mask_ff_endo));
            sd_r2star_hemo_array_endo(j,i) = std(nonzeros(img3 .* hemo_mask_r2star_endo));
        end
        
        if any(hemo_mask_t1_epi(:))
            mean_t1_hemo_array_epi(j,i) = mean(nonzeros(img .* hemo_mask_t1_epi));
            mean_ff_hemo_array_epi(j,i) = mean(nonzeros(img2 .* hemo_mask_ff_epi));
            mean_r2star_hemo_array_epi(j,i) = mean(nonzeros(img3 .* hemo_mask_r2star_epi));
            sd_t1_hemo_array_epi(j,i) = std(nonzeros(img .* hemo_mask_t1_epi));
            sd_ff_hemo_array_epi(j,i) = std(nonzeros(img2 .* hemo_mask_ff_epi));
            sd_r2star_hemo_array_epi(j,i) = std(nonzeros(img3 .* hemo_mask_r2star_epi));
        end
    end
end

%% Plot one linear regression
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_ff_array_nz = nonzeros(mean_ff_array);
mean_t1_array_nz = nonzeros(mean_t1_array);

mean_t1_hemo_array_nz = nonzeros(mean_t1_hemo_array);
mean_r2star_hemo_array_nz = nonzeros(mean_r2star_hemo_array);

mean_ff_array_nz_new = mean_ff_array_nz;
mean_t1_array_nz_new = mean_t1_array_nz;

mean_t1_hemo_array_nz_new = mean_t1_hemo_array_nz;
mean_r2star_hemo_array_nz_new = mean_r2star_hemo_array_nz;

mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];

mean_t1_hemo_array_nz_new(mean_r2star_hemo_array_nz<=0) = [];
mean_r2star_hemo_array_nz_new(mean_r2star_hemo_array_nz<=0) = [];

mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote);
mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote);
mean_r2star_array_nz_remote = nonzeros(mean_r2star_array_remote);

mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
mean_r2star_array_nz_new_remote = mean_r2star_array_nz_remote;

mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
mean_r2star_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];

plot_save = cat(2, tp_dir2, 'LinearReg/');
if ~exist(plot_save, 'dir')
    mkdir(plot_save);
end

mdl = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);

if any(mean_r2star_hemo_array_nz_new)
    mdl2 = fitlm(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new);
    
    figure('Position', [100 100 600 800]);
    title(cat(2, name, ' ', time_point))
    scatter(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_r2star_hemo_array_nz_new .* mdl2.Coefficients.Estimate(2) + mdl2.Coefficients.Estimate(1);
    hold on;
    plot(mean_r2star_hemo_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_r2star_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    xlabel('R2star (s^{-1})');
    ylabel('T1 (ms)');
    yl = ylim;
    xl = xlim;
    text(0.5*xl(2), yl(1)+200, cat(2,'Y = ', num2str(mdl2.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl2.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl2.Rsquared.Ordinary,3)), 'FontSize', 12)
    saveas(gcf, cat(2, plot_save, 'r2starVST1_1line_demon_in_hemo.png'));
else
    mdl2 = struct;
    
    figure('Position', [100 100 600 800]);
    title(cat(2, name, ' ', time_point))
    scatter(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = [];
    hold on;
    plot(mean_r2star_hemo_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_r2star_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    xlabel('R2star (s^{-1})');
    ylabel('T1 (ms)');
    yl = ylim;
    xl = xlim;
    text(0.5*xl(2), 0.95*yl(2), 'NA', 'FontSize', 12)
    saveas(gcf, cat(2, plot_save, 'r2starVST1_1line_demon_in_hemo.png'));
end



figure('Position', [100 100 600 800]);
title(cat(2, name, ' ', time_point))
rows = ceil(size(roi_in_myo_t1, 3) / 2) + 1;
subplot(rows,1,rows);
scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_nz_new .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
yl = ylim;
xl = xlim;
text(0.5*xl(2), yl(1)+200, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

mdl_general_ffvst1 = mdl;
mdl_general_r2vst1_in_hemo = mdl2;

mdl_slc_cell = cell(size(roi_in_myo_t1, 3), 1);
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz = nonzeros(mean_ff_array(:,i));
    mean_t1_array_nz = nonzeros(mean_t1_array(:,i));
    mean_ff_array_nz_new = mean_ff_array_nz;
    mean_t1_array_nz_new = mean_t1_array_nz;
    mean_ff_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_t1_array_nz_new(mean_ff_array_nz<=0) = [];
    mean_ff_array_nz_remote = nonzeros(mean_ff_array_remote(:,i));
    mean_t1_array_nz_remote = nonzeros(mean_t1_array_remote(:,i));
    mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
    mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
    mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
    mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<=0) = [];
    
    mdl_slc = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);
    
    subplot(rows,2,i);
    scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_ff_array_nz_new .* mdl_slc.Coefficients.Estimate(2) + mdl_slc.Coefficients.Estimate(1);
    hold on;
    plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    text(0.2*xl(2), yl(1)+100, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)
    
    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
    
    mdl_slc_cell{i} = mdl_slc;
end

saveas(gcf, cat(2, plot_save, 'T1vsFF_1line_demon.png'));

mdl_results = struct;
mdl_results.mdl_general_ffvst1 = mdl;
mdl_results.mdl_general_r2vst1_in_hemo = mdl2;
mdl_results.mdl_slc = mdl_slc_cell;

save(LR_mdl_fname, '-struct', 'mdl_results');

chord_value_results = struct;
chord_value_results.mean_t1_array = mean_t1_array;
chord_value_results.mean_ff_array = mean_ff_array;
chord_value_results.mean_r2star_array = mean_r2star_array;
chord_value_results.mean_t1_array_remote = mean_t1_array_remote;
chord_value_results.mean_ff_array_remote = mean_ff_array_remote;
chord_value_results.mean_r2star_array_remote = mean_r2star_array_remote;
chord_value_results.mean_t1_hemo_array = mean_t1_hemo_array;
chord_value_results.mean_ff_hemo_array = mean_ff_hemo_array;
chord_value_results.mean_r2star_hemo_array = mean_r2star_hemo_array;

chord_value_results.sd_t1_array = sd_t1_array;
chord_value_results.sd_ff_array = sd_ff_array;
chord_value_results.sd_r2star_array = sd_r2star_array;
chord_value_results.sd_t1_array_remote = sd_t1_array_remote;
chord_value_results.sd_ff_array_remote = sd_ff_array_remote;
chord_value_results.sd_r2star_array_remote = sd_r2star_array_remote;
chord_value_results.sd_t1_hemo_array = sd_t1_hemo_array;
chord_value_results.sd_ff_hemo_array = sd_ff_hemo_array;
chord_value_results.sd_r2star_hemo_array = sd_r2star_hemo_array;
save(chord_values_fname, '-struct', 'chord_value_results');

chord_value_results2 = struct;
chord_value_results2.mean_t1_array_endo = mean_t1_array_endo;
chord_value_results2.mean_ff_array_endo = mean_ff_array_endo;
chord_value_results2.mean_r2star_array_endo = mean_r2star_array_endo;
chord_value_results2.mean_t1_array_remote_endo = mean_t1_array_remote_endo;
chord_value_results2.mean_ff_array_remote_endo = mean_ff_array_remote_endo;
chord_value_results2.mean_r2star_array_remote_endo = mean_r2star_array_remote_endo;
chord_value_results2.mean_t1_hemo_array_endo = mean_t1_hemo_array_endo;
chord_value_results2.mean_ff_hemo_array_endo = mean_ff_hemo_array_endo;
chord_value_results2.mean_r2star_hemo_array_endo = mean_r2star_hemo_array_endo;
chord_value_results2.mean_t1_array_epi = mean_t1_array_epi;
chord_value_results2.mean_ff_array_epi = mean_ff_array_epi;
chord_value_results2.mean_r2star_array_epi = mean_r2star_array_epi;
chord_value_results2.mean_t1_array_remote_epi = mean_t1_array_remote_epi;
chord_value_results2.mean_ff_array_remote_epi = mean_ff_array_remote_epi;
chord_value_results2.mean_r2star_array_remote_epi = mean_r2star_array_remote_epi;
chord_value_results2.mean_t1_hemo_array_epi = mean_t1_hemo_array_epi;
chord_value_results2.mean_ff_hemo_array_epi = mean_ff_hemo_array_epi;
chord_value_results2.mean_r2star_hemo_array_epi = mean_r2star_hemo_array_epi;

chord_value_results2.sd_t1_array_endo = sd_t1_array_endo;
chord_value_results2.sd_ff_array_endo = sd_ff_array_endo;
chord_value_results2.sd_r2star_array_endo = sd_r2star_array_endo;
chord_value_results2.sd_t1_array_remote_endo = sd_t1_array_remote_endo;
chord_value_results2.sd_ff_array_remote_endo = sd_ff_array_remote_endo;
chord_value_results2.sd_r2star_array_remote_endo = sd_r2star_array_remote_endo;
chord_value_results2.sd_t1_hemo_array_endo = sd_t1_hemo_array_endo;
chord_value_results2.sd_ff_hemo_array_endo = sd_ff_hemo_array_endo;
chord_value_results2.sd_r2star_hemo_array_endo = sd_r2star_hemo_array_endo;
chord_value_results2.sd_t1_array_epi = sd_t1_array_epi;
chord_value_results2.sd_ff_array_epi = sd_ff_array_epi;
chord_value_results2.sd_r2star_array_epi = sd_r2star_array_epi;
chord_value_results2.sd_t1_array_remote_epi = sd_t1_array_remote_epi;
chord_value_results2.sd_ff_array_remote_epi = sd_ff_array_remote_epi;
chord_value_results2.sd_r2star_array_remote_epi = sd_r2star_array_remote_epi;
chord_value_results2.sd_t1_hemo_array_epi = sd_t1_hemo_array_epi;
chord_value_results2.sd_ff_hemo_array_epi = sd_ff_hemo_array_epi;
chord_value_results2.sd_r2star_hemo_array_epi = sd_r2star_hemo_array_epi;
save(chord_values_fname2, '-struct', 'chord_value_results2');
end