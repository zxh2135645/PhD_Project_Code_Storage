rootFolder = '/Users/jameszhang/Documents/Data/Anzhen_FID111140/YERJNBSN';
outFolder  = '/Users/jameszhang/Documents/Data/Anzhen_FID111140/GIFs';
rootFolder = '/Users/jameszhang/Documents/Data/Anzhen_FID176224/R21NZM5R';
outFolder  = '/Users/jameszhang/Documents/Data/Anzhen_FID176224/GIFs';
rootFolder = '/Users/jameszhang/Documents/Data/YANG 25881 3534 W STRESS/DICOM/'
outFolder  = '/Users/jameszhang/Documents/Data/YANG 25881 3534 W STRESS/GIFs';
rootFolder = '/Volumes/Extreme SSD/Anzhen/DICOM/New0123456amidicom64/0060349562/YPYKGM05/'
outFolder  = '/Volumes/Extreme SSD/Anzhen/DICOM/New0123456amidicom64/0060349562/GIFs';
rootFolder = '/Users/jameszhang/Documents/Data/26P02_E_AC2/BIRI_Research_Randy_20260216_081058.400000/'
outFolder  = '/Users/jameszhang/Documents/Data/26P02_E_AC2/GIFs';
rootFolder = '/Users/jameszhang/Documents/Data/26P03_D_AC1/BIRI_Research_Randy_20260220_080242.300000/'
outFolder  = '/Users/jameszhang/Documents/Data/26P03_D_AC1/GIFs';
% Only process series whose description contains this (leave empty to process all)
seriesNameFilter = 'cine_sax'; 
seriesNameFilter = 'tf2d12_retro'; 
%seriesNameFilter = 'TF2D12_RETRO'; 


export_cine_dicoms_to_gifs(rootFolder, outFolder, true, seriesNameFilter);
%%
cineData = read_all_cine_dicoms(rootFolder, true, 'cine_sax');
%%
for Nt = 1:size(cineData(1).vol, 4)
    contourData = draw_freehand_contours_on_cine(cineData, 1, Nt);
end
%%
z = 1;
I = cineData(1).vol(:,:,z,contourData.frameIdx);

figure;
imagesc(I); axis image off; colormap gray; hold on;
if contourData.sliceDrawn(z)
    visboundaries(contourData.masks(:,:,z), 'Color', 'r', 'LineWidth', 1);
end
title(sprintf('Slice %d', z));