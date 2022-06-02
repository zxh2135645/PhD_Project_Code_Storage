function [sep_ratio, centroids, mask_centerline, thera_i_array, thera_s_array]= Func_septum_ratio(mask_lv_new, mask_blood_new, xs, ys, xi, yi)

se = strel('disk', 1);
mask_myocardium = mask_lv_new - mask_blood_new;
mask_myocardium(isnan(mask_myocardium)) = 0;

BW_skel = zeros(size(mask_myocardium));
mask_center = zeros(size(mask_myocardium));
mask_centerline = zeros(size(mask_myocardium));
centroids = cell(1, size(mask_myocardium, 3));

for i = 1:size(mask_myocardium, 3)
    if any(any(mask_myocardium(:,:,i)))
        BW_skel(:,:,i) = bwmorph(mask_myocardium(:,:,i), 'skel', Inf);
        mask_center(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
        mask_center(:,:,i) = imopen(mask_center(:,:,i), se);
        mask_center_erode = imerode(mask_center(:,:,i), se);
        mask_centerline(:,:,i) = mask_center(:,:,i) - mask_center_erode;
        
        stats = regionprops(mask_blood_new(:,:,i));
        centroids{i} = stats.Centroid;
    end
end

num_slc = size(mask_centerline, 3);
idx_s = zeros(1, num_slc);
idx_i = zeros(1, num_slc);
thera_i_array = zeros(1, num_slc);
thera_s_array = zeros(1, num_slc);

row_col_sqrt = ceil(sqrt(num_slc));

figure(); 
for slc = 1:num_slc
    if any(any(mask_myocardium(:,:,slc)))
        [y,x] = find(mask_centerline(:,:,slc));
        
        for j = 1:num_slc
            dist_s = sqrt(abs(y - ys(j)).^2 + abs(x - xs(j)).^2);
            [~, idx_s(j)] = min(dist_s);
            dist_i = sqrt(abs(y - yi(j)).^2 + abs(x - xi(j)).^2);
            [~, idx_i(j)] = min(dist_i);
        end
        
        centroid = centroids{slc};
        
        %mask_centerline(y(idx_s(slc)),x(idx_s(slc)),1) = 2;
        %mask_centerline(y(idx_i(slc)),x(idx_i(slc)),1) = 3;
        %mask_centerline(round(centroid(2)),round(centroid(1)),1) = 3;
        
        theta_s = atan2(-(y(idx_s(slc)) - centroid(2)), x(idx_s(slc)) - centroid(1));
        theta_i = atan2(-(y(idx_i(slc)) - centroid(2)), x(idx_i(slc)) - centroid(1));
        
        theta_array = atan2(-(y - centroid(2)), x - centroid(1));
        
        thera_i_array(slc) = atan2(-(yi(slc) - centroid(2)), xi(slc) - centroid(1));
        thera_s_array(slc) = atan2(-(ys(slc) - centroid(2)), xs(slc) - centroid(1));
        
        n = length(theta_array);
        if theta_s < 0 theta_s = theta_s + 2*pi; end
        if theta_i < 0 theta_i = theta_i + 2*pi; end
        
        for i = 1:n
            if theta_array(i) < 0 theta_array(i) = theta_array(i) + 2*pi; end
        end
        

        idx_theta = find(theta_array <= theta_i & theta_array >= theta_s);
        idx_theta_2 = find(theta_array >= theta_i & theta_array <= theta_s);
        
        if length(idx_theta) > 1 && length(idx_theta_2) <= 1
            for i = 1:length(idx_theta)
                mask_centerline(y(idx_theta(i)), x(idx_theta(i)), slc) = 2;
            end
        elseif length(idx_theta) <= 1 && length(idx_theta_2) > 1
            for i = 1:length(idx_theta_2)
                mask_centerline(y(idx_theta_2(i)), x(idx_theta_2(i)), slc) = 2;
            end
        end
        
        subplot(row_col_sqrt,row_col_sqrt,slc);
        imagesc(mask_centerline(:,:,slc)); axis image; axis off;
    end
end

septum = zeros(1, num_slc);
rest = zeros(1, num_slc);
for slc = 1:num_slc
   rest(slc) = sum(sum(mask_centerline(:,:,slc) == 1));
   septum(slc) = sum(sum(mask_centerline(:,:,slc) == 2));
   
   if rest(slc) == 0
      sep_ratio(slc) = nan;
   else
      sep_ratio(slc) = septum(slc) ./ rest(slc);
   end  
end

if any(sep_ratio>2)
    septum = zeros(1, num_slc);
    rest = zeros(1, num_slc);
    for slc = 1:num_slc
        rest(slc) = sum(sum(mask_centerline(:,:,slc) == 2));
        septum(slc) = sum(sum(mask_centerline(:,:,slc) == 1));
        if rest(slc) && septum(slc) == 0
            sep_ratio(slc) = nan;
        elseif septum(slc) ~= 0 && rest(slc) == 0
            sep_ratio(slc) = 0;
        else
            sep_ratio(slc) = septum(slc) ./ rest(slc);
        end
    end
end


