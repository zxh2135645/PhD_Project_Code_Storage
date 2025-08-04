function coords = GroovePickGeneric(vol_image, mask_heart, lb, ub)


p = 1; % number of landmarks to pick
coords = struct;
x_array = zeros(size(mask_heart, 3), p);
y_array = zeros(size(mask_heart, 3), p);
x_centroid_array = zeros(size(mask_heart, 3), 1);
y_centroid_array = zeros(size(mask_heart, 3), 1);


for slice_num = 1:size(vol_image,3)

    [x_heart, y_heart] = find(mask_heart(:,:,slice_num) ~= 0);


    x_centroid_array(slice_num) = round(mean(x_heart),1);
    y_centroid_array(slice_num) = round(mean(y_heart),1);

    if isempty(x_heart) || isempty(y_heart)
        y_array(slice_num) = nan;
        x_array(slice_num) = nan;
    else
        figure();
        imagesc(vol_image(:,:,slice_num));
        caxis([lb, ub]);
        truesize([3*size(vol_image,1), 3*size(vol_image,2)]);
        axis off;

        for i = 1:p
            [y_array(slice_num,i), x_array(slice_num,i)] = ginput(1);
            y_array(slice_num,i) = round(y_array(slice_num,i), 1);
            x_array(slice_num,i) = round(x_array(slice_num,i), 1);
            h1 = text(y_array(slice_num,i), x_array(slice_num,i), '*', 'HorizontalAlignment', 'center', 'Color', [0 0 0], 'FontSize', 16);
            h2 = text(y_array(slice_num,i), x_array(slice_num,i), num2str(i), 'HorizontalAlignment', 'center', 'Color', [1 0 0], 'FontSize', 14);
        end
    end
end

close all;
coords.x_array = x_array;
coords.y_array = y_array;
coords.x_centroid_array = x_centroid_array;
coords.y_centroid_array = y_centroid_array;


end