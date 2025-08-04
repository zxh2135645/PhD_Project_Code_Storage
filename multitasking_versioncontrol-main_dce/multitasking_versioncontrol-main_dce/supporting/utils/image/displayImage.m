function displayImage(img,ylim,cmap,slice)

img = squeeze(abs(img));
img = img(:,:,:);

if nargin < 4
    slice = floor(size(img,3)/2) + 1;
elseif slice < 1 || slice > size(img,3)
    slice = floor(size(img,3)/2) + 1;
end

if nargin < 3
    cmap = 'gray';
end
if nargin < 2
    ylim = [-Inf Inf];
elseif isempty(ylim)
    ylim = [-Inf Inf];
end

imagesc(img(:,:,slice),ylim);axis equal tight;colormap(cmap);
