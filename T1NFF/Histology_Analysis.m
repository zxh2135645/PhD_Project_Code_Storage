clear all;
close all;

addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
dir_name = cat(2, base_dir, '/Vol_Histology/');
fname = 'Merry_Trichrome-1_magnified.PNG';

f_to_read = imread(cat(2, dir_name, fname));

%%
im = im2double(f_to_read);
figure();
imagesc(im);

numColors = 3;
L = imsegkmeans(f_to_read,numColors);
B = labeloverlay(f_to_read,L);
imshow(B)
title("Labeled Image RGB");

% Step 3: Convert Image from RGB Color Space to L*a*b* Color Space
lab_he = rgb2lab(f_to_read);

%%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
pixel_labels = imsegkmeans(ab,numColors,NumAttempts=3);

figure();
B2 = labeloverlay(f_to_read,pixel_labels);
imshow(B2)
title("Labeled Image a*b*")

%% Step 5: Create Images that Segment H&E Image by Color
figure();
mask1 = pixel_labels == 1;
cluster1 = f_to_read.*uint8(mask1);
imshow(cluster1)
title("Objects in Cluster 1");

figure();
mask2 = pixel_labels == 2;
cluster2 = f_to_read.*uint8(mask2);
imshow(cluster2)
title("Objects in Cluster 2");

figure();
mask3 = pixel_labels == 3;
cluster3 = f_to_read.*uint8(mask3);
imshow(cluster3)
title("Objects in Cluster 3");

% figure();
% mask4 = pixel_labels == 4;
% cluster4 = f_to_read.*uint8(mask4);
% imshow(cluster4)
% title("Objects in Cluster 4");

%% Another k-means?
% pixel_labels_2 = imsegkmeans(ab.*mask1,numColors,NumAttempts=3);
% 
% figure();
% mask1_2 = pixel_labels_2 == 1;
% cluster1 = f_to_read.*uint8(mask1_2);
% imshow(cluster1)
% title("Objects in Cluster 1");
% 
% figure();
% mask2_2 = pixel_labels_2 == 2;
% cluster2 = f_to_read.*uint8(mask2_2);
% imshow(cluster2)
% title("Objects in Cluster 2");
% 
% figure();
% mask3_2 = pixel_labels == 3;
% cluster3 = f_to_read.*uint8(mask3_2);
% imshow(cluster3)
% title("Objects in Cluster 3");

%%

figure(); imshow(mask1, []);

se = strel('disk',2);
afterOpening = imopen(mask1,se);
figure();
imshow(uint8(afterOpening).*f_to_read,[]);

%%
BW = ~mask1;
BW2 = imfill(BW,"holes");
figure
imshow(~BW2)
title('Filled Image');

se = strel('disk',5);
closeBW = imclose(BW2,se);
figure, imshow(closeBW)

[X Y] = size(closeBW);
C1 = regiongrowing(closeBW, 1, 10);
figure();
imshow(C1);

C2 = regiongrowing(closeBW, X, Y);
figure();
imshow(C2);

%%
ff_1 = uint8(~(C1+C2)) .* B2 .* uint8(~(~afterOpening|BW2));

figure();
imshow(ff_1);

a = X*Y - numel(nonzeros(C1+C2))
b = numel(nonzeros(~afterOpening|BW2))

ff = (a - b) / a


E = entropyfilt(ff_1);
S = stdfilt(ff_1,ones(9));
R = rangefilt(ff_1,ones(9));

Eim = rescale(E);
Sim = rescale(S);

figure();
montage({Eim,Sim,R},'Size',[1 3],'BackgroundColor','w',"BorderSize",20)
title('Texture Images Showing Local Entropy, Local Standard Deviation, and Local Range')


BW11 = imbinarize(rgb2gray(Eim),0.5);
figure();
imshow(BW11)
title('Thresholded Texture Image')

figure();
imshow(uint8(BW11).*B2)

%%
se = strel('disk',5);
BW11_afteropening = imopen(BW11,se);
figure();
imshow(uint8(BW11_afteropening).*B2)
%%
b2 = numel(nonzeros(BW11))
ff2 = b2 / a

b2 = numel(nonzeros(BW11_afteropening))
ff2 = b2 / a

BWao = bwareaopen(BW11_afteropening,300);
figure();
imshow(BWao)
title('Area-Opened Texture Image')

BWao_fill = imfill(BWao,"holes");

b2 = numel(nonzeros(BWao_fill))
ff2 = b2 / a

figure();
imshow(uint8(BWao_fill) .* f_to_read);

figure();
imshow(uint8(mask3) .* f_to_read);

%%
ff_1_mask = ~(C1+C2) .* ~(~afterOpening|BW2);

figure();
imshow(ff_1_mask);

ff_1_mask_BWao = bwareaopen(ff_1_mask,20);
figure();
imshow(ff_1_mask_BWao)
title('Area-Opened Texture Image')

b3 = numel(nonzeros(ff_1_mask_BWao))
ff3 = b3 / a

ff_1_mask_BWao_fill = imfill(ff_1_mask_BWao,"holes");
figure();
imshow(uint8(ff_1_mask_BWao_fill) .* f_to_read);

b4 = numel(nonzeros(ff_1_mask_BWao_fill))
ff4 = b4 / a

% P = 500;
% BW3 = bwareaopen(~BW2,P);
% figure
% imshow(BW3);
% title('Removed Image');

% [C,h] = imcontour(~BW2,10);
% 
% figure(); 
% imcontour(I,[128 128]); 

%% Step 6: Segment Nuclei
% L = lab_he(:,:,1);
% L_blue = L.*double(mask1);
% L_blue = rescale(L_blue);
% idx_light_blue = imbinarize(nonzeros(L_blue));
% 


% blue_idx = find(mask1);
% mask_dark_blue = mask1;
% mask_dark_blue(blue_idx(idx_light_blue)) = 0;
% 
% figure();
% blue_nuclei = f_to_read.*uint8(mask_dark_blue);
% imshow(blue_nuclei)
% title("Blue Nuclei")
%%
% figure();
% imagesc(im_gray>0.75);


% K = wiener2(im_gray,[3 3]);
% figure;
% imagesc(K);
% title('Portion of the Image with Noise Removed by Wiener Filter');