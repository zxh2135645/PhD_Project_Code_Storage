rootFolder = '/Users/jameszhang/Documents/Data/26P03_Heart_Exvivo_Phantom_040726P03_19900407/DICOM';
outFolder  = '/Users/jameszhang/Documents/Data/26P03_Heart_Exvivo_Phantom_040726P03_19900407/GIFs';


% Only process series whose description contains this (leave empty to process all)
seriesNameFilter = '3DmGRE'; 

for i = 1:8
    export_cine_dicoms_to_gifs(rootFolder, outFolder, true, seriesNameFilter, i);
end
%%
cineData = read_all_cine_dicoms(rootFolder, true, seriesNameFilter);
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