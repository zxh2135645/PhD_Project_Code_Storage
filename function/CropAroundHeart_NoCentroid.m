function [img_cropped] = CropAroundHeart_NoCentroid(img)
% No centroid
img_size = size(img);
if size(img, 3) > 1
    % img_cropped = cell(img_size(3), 1);
    for i = 1:img_size(3)
        img_2D = img(:,:,i);
        if img_size(1) <= img_size(2)
            x_dis = round(img_size(2)/6);
            y_dis = round(img_size(1)/6);
            img_cropped_2d = imcrop(img_2D, [round(0.5*x_dis), round(0.5*y_dis), 5*x_dis, 5*y_dis]);
        else
            x_dis = round(img_size(2)/6);
            y_dis = round(img_size(1)/6);
            img_cropped_2d = imcrop(img_2D, [round(0.5*x_dis), round(0.5*y_dis), 5*x_dis, 5*y_dis]);
        end
        img_cropped(:,:,i) = img_cropped_2d;
    end
    
else
    if img_size(1) == img_size(2)
        x_dis = round(img_size(2)/6);
        y_dis = round(img_size(1)/6);
        img_cropped = imcrop(img, [round(1.5*x_dis), round(1.5*y_dis), 3*x_dis, 3*y_dis]);
    elseif img_size(1) < img_size(2)
        x_dis = round(img_size(2)/6);
        y_dis = round(img_size(1)/6);
        img_cropped = imcrop(img, [round(1.5*x_dis), round(1.5*y_dis), 3*x_dis, 3*y_dis]);
    else
        x_dis = round(img_size(2)/6);
        y_dis = round(img_size(1)/6);
        img_cropped = imcrop(img, [x_dis, y_dis, 3*x_dis, 3*y_dis]);
    end
end
end