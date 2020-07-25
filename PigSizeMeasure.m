clear all;
close all;

[volume_image, slice_data, image_meta_data] = dicom23D('D:\Data\Cedars_2020_pig\DHARMAKUMAR_5257_20P10_WEEK1\BIRI_RESEARCH_DIANE_20200427_085459_137000\TFL_LOC_MULTI_IPAT_TRA_0004\');
res = slice_data(1).PixelSpacing(2);
%%
figure();
for i = 1:size(volume_image,3)
    subplot(3,5,i)
    imagesc(volume_image(:,:,i)); axis image;
    
end

%%
img = volume_image(:,:,end);
figure();
imagesc(img); axis image
[x, y] = ginput(2);
%%
res*(x(2) - x(1))

% 217 mm

%%
img = volume_image(:,:,end);
figure();
imagesc(img); axis image
[x2, y2] = ginput(2);

%%
res*(y2(2) - y2(1))

% is about 300 mm