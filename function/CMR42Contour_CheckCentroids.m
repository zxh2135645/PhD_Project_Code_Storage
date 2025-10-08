   function [endo_flow, epi_flow, excludeContour, myoRef, NoReFlow, freeROI] = ...
        CMR42Contour_CheckCentroids(num_slice, endo_flow, epi_flow, excludeContour, myoRef, NoReFlow, freeROI, volume_image);
    
        endo_centroids = zeros(num_slice, 2);
        epi_centroids = zeros(num_slice, 2);
        exclude_centroids = zeros(num_slice, 2);
        myoRef_centroids = zeros(num_slice, 2);
        noReFlow_centroids = zeros(num_slice, 2);
        freeROI_centroids = zeros(num_slice, 2);
        
    for i = 1:num_slice
        % Get centroid of endo_flow
        if ~isempty(endo_flow{i})
            coords = endo_flow{i}; % [N,2]
            endo_centroids(i, :) = mean(coords, 1); % [1,2] centroid (X,Y)
        end
        % Get centroid of epi_flow
        if ~isempty(epi_flow{i})
            coords = epi_flow{i};
            epi_centroids(i, :) = mean(coords, 1);
        end
        % Get centroid of excludeContour
        if ~isempty(excludeContour{i})
            coords = excludeContour{i};
            exclude_centroids(i, :) = mean(coords, 1);
        end
        % Get centroid of myoRef
        if ~isempty(myoRef{i})
            coords = myoRef{i};
            myoRef_centroids(i, :) = mean(coords, 1);
        end
        % Get centroid of NoReFlow
        if ~isempty(NoReFlow{i})
            coords = NoReFlow{i};
            noReFlow_centroids(i, :) = mean(coords, 1);
        end
        % Get centroid of freeROI
        if ~isempty(freeROI{i})
            coords = freeROI{i};
            freeROI_centroids(i, :) = mean(coords, 1);
        end
    end

        % Check centroids for each contour type
    centroid_types = {'endo_centroids', 'epi_centroids', 'exclude_centroids', 'myoRef_centroids', 'noReFlow_centroids', 'freeROI_centroids'};
    centroid_vars = {endo_centroids, epi_centroids, exclude_centroids, myoRef_centroids, noReFlow_centroids, freeROI_centroids};
    coords_vars = {endo_flow, epi_flow, excludeContour, myoRef, NoReFlow, freeROI};
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

endo_flow = coords_vars{1};
epi_flow = coords_vars{2};
excludeContour = coords_vars{3};
myoRef = coords_vars{4};   
NoReFlow = coords_vars{5};
freeROI = coords_vars{6};

    end
