function [img_cropped] = CropAroundHeart(centroid, img)
% centroid is 1x2
x = centroid(1);
y = centroid(2);
x_dis = round(size(img,2)/5);
y_dis = round(size(img,1)/5);
img_cropped = imcrop(img, [x-x_dis, y-y_dis, 2*x_dis, 2*y_dis]);

end