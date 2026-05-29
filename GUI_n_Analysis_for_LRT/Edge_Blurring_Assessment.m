
im_LRT = imread('/Users/jameszhang/Documents/MATLAB/LRT/Results/Lisbon_D6/Resp_motion_V2Tiffs/FID17210_Lisbon_D6__LRT_CINE_SA5_RespPhase2_WholeHeart_DB.tif');
imgDist = double(im_LRT);

im_CMR = dicomread('/Users/jameszhang/Documents/MATLAB/LRT/DICOM/LISBON_D6/TFL25_T1_PSIR_SAX6_MAG_0146/LISBON_DAY6.MR.HEART_CEDARS.0146.0001.2021.11.15.12.46.40.753818.16159629.IMA');
imgRef = double(im_CMR);

% remove first and last 31 columns (y-direction)
cols = size(imgDist,2);
if cols > 62
    imgDist = imgDist(:, 32:cols-31);
else
    error('Cannot crop: image has %d columns, need >62.', cols);
end
figure(); imagesc(imgDist); axis image; title('Cropped imgDist');

% resize imgDist to match imgRef
[tr,tc,~] = size(imgRef);
[dr,dc,~] = size(imgDist);
if dr ~= tr || dc ~= tc
    imgDist = imresize(imgDist, [tr tc], 'bicubic');
end

% show the resized result
figure; imshow(imgDist, []); title('Interpolated imgDist');


% rescale imgDist and imgRef to the range [0,255]
mn = min(imgDist(:)); mx = max(imgDist(:));
if mx > mn
    imgDist = (imgDist - mn) / (mx - mn) * 255;
else
    imgDist = zeros(size(imgDist));
end

mn = min(imgRef(:)); mx = max(imgRef(:));
if mx > mn
    imgRef = (imgRef - mn) / (mx - mn) * 255;
else
    imgRef = zeros(size(imgRef));
end

% convert to uint8
imgDist = uint8(round(imgDist));
imgRef  = uint8(round(imgRef));
%%
figure(); imagesc(imgDist); axis image;
figure(); imagesc(imgRef); axis image;

%% draw binary mask around heart on the reference image
figure; imshow(imgRef, []); title('Draw heart ROI (double-click to finish)');
hROI = drawfreehand('Color','r','LineWidth',1.5);
mask = createMask(hROI);

% show results
figure; imshow(mask); title('Binary mask (heart region)');
figure; imshowpair(mat2gray(imgRef), mask, 'blend'); title('Heart mask overlay');

%% optional: apply mask to reference image
imgRefMasked = imgRef;
imgRefMasked(~mask) = 0;
figure; imagesc(imgRefMasked); title('Masked reference image');

imgDistMasked = imgDist;
imgDistMasked(~mask) = 0;
figure; imagesc(imgDistMasked); title('Masked LRT image');
%%
[similarity,similarityMaps,weightMaps] = HaarPSI(imgRefMasked,imgDistMasked,1);

