function mask = genMask2D(img,Nbins)

if nargin < 2
    Nbins = 3;
end

img = abs(img);
if ndims(img) > 3
    img = img(:,:,:,1);
end

mask = zeros(size(img));
for slice = 1:size(img,3)
    tempimg = img(:,:,slice);
    tempimg = imgaussfilt(tempimg,2);
    tempimg = imfill(tempimg,8);
    tempimg(tempimg>prctile(tempimg(:),80)) = prctile(tempimg(:),80);
    [im_mask,centroids] = kmeans(tempimg(:),Nbins);
    [~,idx] = min(centroids);
    tempimg(im_mask==idx) = 0;
    tempimg(tempimg>0) = 1;
    tempimg = imfill(tempimg,8);
    tempimg = imopen(tempimg,strel('disk',2));
    mask(:,:,slice) = tempimg;
end
