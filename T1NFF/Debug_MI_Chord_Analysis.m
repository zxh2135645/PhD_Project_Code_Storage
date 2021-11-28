%for i = 4:4

%% This is the one to generate rotating segment video
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    img3 = r2star(:,:,i);
    %fixed = myo_t1(:,:,i);
    %moving = myo_ff(:,:,i);
    fixed = roi_in_myo_t1(:,:,i);
    moving = roi_in_myo_ff(:,:,i);
    
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    img2 = imwarp(img2,tform,'OutputView',imref2d(size(fixed)));
    img3 = imwarp(img3,tform,'OutputView',imref2d(size(fixed)));

    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
    
    movingRegistered_myo_ff = imwarp(myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    movingRegistered_roi_ff = movingRegistered;
    movingRegistered_roi_r2star = movingRegistered;
    %movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    
    fixed_eroded = imerode(myo_t1(:,:,i), se);
    BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
    center_mask_t1(:,:,i) = imfill(BW_skel, 'hole');
    center_fixed = center_mask_t1(:,:,i);
    center_fixed = imopen(center_fixed, se); % Removing spikes
    fixedRegistered_epi = fixed_eroded - center_fixed > 0;
    fixedRegistered_endo = center_fixed + fixed_eroded > 1;
    
    movingRegistered_eroded = imerode(movingRegistered_myo_ff, se);
    BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
    center_mask_ff(:,:,i) = imfill(BW_skel, 'hole');
    center_moving = center_mask_ff(:,:,i);
    center_moving = imopen(center_moving, se); % Removing spikes
    movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
    movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
    
    %movingRegistered_roi_ff = imwarp(roi_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    %movingRegistered_roi_r2star = imwarp(roi_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
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
        
    figure('Position', [100 0 1600 1600]);
    seg_mask_overlap = zeros(size(img));
    vid_save = cat(2, 'myVideoFile.mov');
    myVideo = VideoWriter(vid_save); %open video file
    myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
    
    open(myVideo)
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
        
         % For debugging
        seg_mask_t1 = Mask_Segn .* roi_in_myo_t1(:,:,i) .* fixed_eroded == j;
        seg_mask_ff = Mask_Segn2 .* movingRegistered_roi_ff .* movingRegistered_eroded == j;
        
        if ~isempty(nonzeros(seg_mask_t1)) && ~isempty(nonzeros(seg_mask_ff))
            seg_mask_overlap = seg_mask_overlap + seg_mask_t1 + 2*seg_mask_ff;
            ax1 = subplot(1,2,1);
            imagesc(ax1, seg_mask_overlap); axis image; caxis([0 3]);
            title(['Segment = ', num2str(j)]);
            colormap(brewermap([],'*RdYlBu'));
            ax2 = subplot(1,2,2);
            imagesc(ax2, seg_mask_t1 + 2*seg_mask_ff); axis image; caxis([0 3]);
            title(['Segment = ', num2str(j)]); %colorbar;
            colormap(brewermap([],'*RdYlBu'));
            pause(.5);
            frame = getframe(gcf); %get frame
            writeVideo(myVideo, frame);
        end
    end
    close(myVideo);