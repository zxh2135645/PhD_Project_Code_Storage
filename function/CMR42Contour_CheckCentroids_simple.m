function [excludeContour] = ...
    CMR42Contour_CheckCentroids_simple(num_slice, excludeContour, volume_image);


exclude_centroids = zeros(num_slice, 2);


for i = 1:num_slice

    % Get centroid of excludeContour
    if ~isempty(excludeContour{i})
        coords = excludeContour{i};
        exclude_centroids(i, :) = mean(coords, 1);
    end
end

% Check centroids for each contour type
centroid_types = {'exclude_centroids'};
centroid_vars = {exclude_centroids};
coords_vars = {excludeContour};
img_height = size(volume_image, 1);
img_width = size(volume_image, 2);

for k = 1:numel(centroid_types)
    centroids = centroid_vars{k};
    for i = 1:num_slice
        if any(centroids(i,:))
            x = centroids(i,1);
            y = centroids(i,2);
            if x < img_height/8 || x > 7*img_height/8 || y < img_width/8 || y > 7*img_width/8
                warning('%s centroid for slice %d is out of expected range: (%.2f, %.2f)', centroid_types{k}, i, x, y);

                if x < img_height/8
                    coords_vars{k}{i}(:,1) =  coords_vars{k}{i}(:,1) * (8);
                elseif x > 7*img_height/8
                    coords_vars{k}{i}(:,1) =  coords_vars{k}{i}(:,1) * (1/8);
                end
                if y < img_width/8
                    coords_vars{k}{i}(:,2) =  coords_vars{k}{i}(:,2) * (8);
                elseif y > 7*img_width/8
                    coords_vars{k}{i}(:,2) =  coords_vars{k}{i}(:,2) * (1/8);
                end
                centroid_vars{k}(i,:) = mean(coords_vars{k}{i}, 1);
            end
        end
    end
end


excludeContour = coords_vars{1};


end
