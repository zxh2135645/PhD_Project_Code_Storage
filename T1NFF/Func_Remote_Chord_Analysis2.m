function [Remote_Chord_Analysis2, center_mask_ff] = ...
    Func_Remote_Chord_Analysis2(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star, ...
    Remote_Chord_Analysis_fname)
% Ripped off from Func_MI_Chord_Analysis2
% Making ff10 masks and R2*80 masks
Remotepix = cell(Segn, size(t1, 3));
Remotepix_epi = cell(Segn, size(t1, 3));
Remotepix_endo = cell(Segn, size(t1, 3));

Remotepix2 = cell(Segn, size(t1, 3));
Remotepix2_epi = cell(Segn, size(t1, 3));
Remotepix2_endo = cell(Segn, size(t1, 3));

Remotepix3 = cell(Segn, size(t1, 3));
Remotepix3_epi = cell(Segn, size(t1, 3));
Remotepix3_endo = cell(Segn, size(t1, 3));

se = strel('disk', 1);
[optimizer, metric] = imregconfig('multimodal');
Remote_Chord_Analysis2 = struct;
center_mask_ff = zeros(size(remote_in_myo_t1));
center_mask_t1 = zeros(size(remote_in_myo_t1));

for i = 1:size(t1, 3)
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    img3 = r2star(:,:,i);
    fixed = myo_t1(:,:,i);
    moving = myo_ff(:,:,i);
    
    % ff is moving and t1 is fixed
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    img2 = imwarp(img2,tform,'OutputView',imref2d(size(fixed)));
    img3 = imwarp(img3,tform,'OutputView',imref2d(size(fixed)));
    
    movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
    
    fixed_eroded = imerode(fixed, se);
    %BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
    %center_mask_t1(:,:,i) = imfill(BW_skel, 'hole');
    %center_fixed = center_mask_t1(:,:,i);
    %center_fixed = imopen(center_fixed, se); % Removing spikes
    %fixedRegistered_epi = fixed_eroded - center_fixed > 0;
    %fixedRegistered_endo = center_fixed + fixed_eroded > 1;
    
    movingRegistered_eroded = imerode(movingRegistered, se);
    BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
    center_mask_ff(:,:,i) = imfill(BW_skel, 'hole');
    center_moving = center_mask_ff(:,:,i);
    center_moving = imopen(center_moving, se); % Removing spikes
    movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
    movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
    
    movingRegistered_remote_ff = imwarp(remote_in_myo_ff(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    movingRegistered_remote_r2star = imwarp(remote_in_myo_r2star(:,:,i),tform,'OutputView',imref2d(size(fixed))) > 0.5;
    %movingRegistered_epi = imwarp(myo_epi_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    %movingRegistered_endo = imwarp(myo_endo_moving,tform,'OutputView',imref2d(size(fixed))) > 0.5;
    
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(img, movingRegistered_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2] = AHASegmentation(img2, movingRegistered_eroded, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3] = AHASegmentation(img3, movingRegistered_eroded, Segn, Groove);
    
    [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(img, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(img, movingRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_epi] = AHASegmentation(img2, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn2_endo] = AHASegmentation(img2, movingRegistered_endo, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_epi] = AHASegmentation(img3, movingRegistered_epi, Segn, Groove);
    [Segmentpix, stats, Mask_Segn3_endo] = AHASegmentation(img3, movingRegistered_endo, Segn, Groove);
    
    for j = 1:Segn
        Remotepix{j,i} = img(Mask_Segn .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        Remotepix_epi{j,i} = img(Mask_Segn_epi .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        Remotepix_endo{j,i} = img(Mask_Segn_endo .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        
        Remotepix2{j,i} = img2(Mask_Segn2 .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        Remotepix2_epi{j,i} = img2(Mask_Segn2_epi .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        Remotepix2_endo{j,i} = img2(Mask_Segn2_endo .* movingRegistered_remote_ff .* movingRegistered_eroded == j);
        
        Remotepix3{j,i} = img3(Mask_Segn3 .* movingRegistered_remote_r2star .* movingRegistered_eroded == j);
        Remotepix3_epi{j,i} = img3(Mask_Segn3_epi .* movingRegistered_remote_r2star .* movingRegistered_eroded == j);
        Remotepix3_endo{j,i} = img3(Mask_Segn3_endo .* movingRegistered_remote_r2star .* movingRegistered_eroded == j);
    end
    
    Remote_idx1 = [];
    Remote_idx1_epi = [];
    Remote_idx1_endo = [];
    for j = 1:Segn
        if ~isempty(Remotepix{j,i})
            Remote_idx1 = [Remote_idx1, j];
        end
        
        if ~isempty(Remotepix_epi{j,i})
            Remote_idx1_epi = [Remote_idx1_epi, j];
        end
        
        if ~isempty(Remotepix_endo{j,i})
            Remote_idx1_endo = [Remote_idx1_endo, j];
        end
    end
    
    Remote_idx2 = [];
    Remote_idx2_epi = [];
    Remote_idx2_endo = [];
    for j = 1:Segn
        if ~isempty(Remotepix2{j,i})
            Remote_idx2 = [Remote_idx2, j];
        end
        if ~isempty(Remotepix2_epi{j,i})
            Remote_idx2_epi = [Remote_idx2_epi, j];
        end
        if ~isempty(Remotepix2_endo{j,i})
            Remote_idx2_endo = [Remote_idx2_endo, j];
        end
    end
    
     % flatten and get the mean 
    Remote_idx_intersect = intersect(Remote_idx1, Remote_idx2);
    Remote_idx_intersect_epi = intersect(Remote_idx1_epi, Remote_idx2_epi);
    Remote_idx_intersect_endo = intersect(Remote_idx1_endo, Remote_idx2_endo);
    
     % Save Remote_Chord_Analysis2
    Remote_Chord_Analysis2(i).Mask_Segn = Mask_Segn; 
    Remote_Chord_Analysis2(i).Mask_Segn2 = Mask_Segn2; 
    Remote_Chord_Analysis2(i).Mask_Segn3 = Mask_Segn3;
    
    Remote_Chord_Analysis2(i).Remotepix = Remotepix; 
    Remote_Chord_Analysis2(i).Remotepix2 = Remotepix2; 
    Remote_Chord_Analysis2(i).Remotepix3 = Remotepix3;

    % Epi
    Remote_Chord_Analysis2(i).Mask_Segn_epi = Mask_Segn_epi; 
    Remote_Chord_Analysis2(i).Mask_Segn2_epi = Mask_Segn2_epi; 
    Remote_Chord_Analysis2(i).Mask_Segn3_epi = Mask_Segn3_epi;
    
    Remote_Chord_Analysis2(i).Remotepix_epi = Remotepix_epi; 
    Remote_Chord_Analysis2(i).Remotepix2_epi = Remotepix2_epi; 
    Remote_Chord_Analysis2(i).Remotepix3_epi = Remotepix3_epi;
    
    % Endo
    Remote_Chord_Analysis2(i).Mask_Segn_endo = Mask_Segn_endo; 
    Remote_Chord_Analysis2(i).Mask_Segn2_endo = Mask_Segn2_endo; 
    Remote_Chord_Analysis2(i).Mask_Segn3_endo = Mask_Segn3_endo;
    
    Remote_Chord_Analysis2(i).Remotepix_endo = Remotepix_endo; 
    Remote_Chord_Analysis2(i).Remotepix2_endo = Remotepix2_endo; 
    Remote_Chord_Analysis2(i).Remotepix3_endo = Remotepix3_endo;
    
    Remote_Chord_Analysis2(i).fixed_eroded = fixed_eroded;
    Remote_Chord_Analysis2(i).movingRegistered_eroded = movingRegistered_eroded;
    
    % Mean values
    t1_mean = mean(nonzeros(img .* movingRegistered_remote_ff .* movingRegistered_eroded));
    t1_sd = std(nonzeros(img .* movingRegistered_remote_ff .* movingRegistered_eroded));
    
    img2(img2 > 100) = 100;
    img2(img2 < 0) = 0;
    ff_mean = mean(nonzeros(img2 .* movingRegistered_remote_ff .* movingRegistered_eroded));
    ff_sd = std(nonzeros(img2 .* movingRegistered_remote_ff .* movingRegistered_eroded));
    
    img3(img3 < 0) = 0;
    r2star_mean = mean(nonzeros(img3 .* movingRegistered_remote_r2star .* movingRegistered_eroded));
    r2star_sd = std(nonzeros(img3 .* movingRegistered_remote_r2star .* movingRegistered_eroded));
    
    Remote_Chord_Analysis2(i).t1_mean = t1_mean;
    Remote_Chord_Analysis2(i).ff_mean = ff_mean;
    Remote_Chord_Analysis2(i).r2star_mean = r2star_mean;
    
    Remote_Chord_Analysis2(i).t1_sd = t1_sd;
    Remote_Chord_Analysis2(i).ff_sd = ff_sd;
    Remote_Chord_Analysis2(i).r2star_sd = r2star_sd;

end
    save(Remote_Chord_Analysis_fname, 'Remote_Chord_Analysis2');
end