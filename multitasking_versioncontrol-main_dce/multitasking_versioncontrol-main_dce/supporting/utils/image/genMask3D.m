function mask = genMask3D(img,Nbins)

if nargin < 2
    Nbins = 4;
end

img = abs(img);
if ndims(img) > 3
    img = img(:,:,:,1);
end

tempimg = imgaussfilt3(img,2);
tempimg = imfill(tempimg,8);
tempimg(tempimg>prctile(tempimg(:),80)) = prctile(tempimg(:),80);
[im_mask,centroids] = kmeans(tempimg(:),Nbins);
[~,idx] = min(centroids);
tempimg(im_mask==idx) = 0;
tempimg(tempimg>0) = 1;
tempimg = imfill(tempimg,8);
tempimg = imfill(tempimg,26);
tempimg = imopen(tempimg,strel('disk',2));

mask = tempimg;

