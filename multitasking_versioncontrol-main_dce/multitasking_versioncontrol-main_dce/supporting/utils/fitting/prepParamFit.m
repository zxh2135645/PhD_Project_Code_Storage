function prepParamFit(prefix, img, tIdx, echoSpacing, flipAngle)
% Create image binary and header file forT1 fitting
%
% Input:
%
%   prefix      - file name prefix for image and header
%   img         - 4D image matrix
%   tIdx        - time index for parametric fitting, should be the same
%                 length as the last dimension of img
%
%   -- Hsu-Lei Lee, last edited on 11/16/2021


%% Extract header info

if nargin < 3
    tIdx = 1:size(img, length(size(img)));
end

img = reshape(img, [], numel(tIdx)).';
iNvoxel = size(img, 2);

tIdxstring = [];
for n = 1:numel(tIdx)
    tIdxstring = [tIdxstring num2str(tIdx(n)) ','];
end
tIdxstring = tIdxstring(1:end-1);

% write to header file
filenameHdr = [prefix '_params'];
fileID      = fopen(filenameHdr,'w');

fprintf(fileID,'iLength=%d\n' , size(img,1));
fprintf(fileID,'iNvoxel=%d\n' , iNvoxel);
fprintf(fileID,'timeIdx=%s\n' , tIdxstring);
fprintf(fileID,'echoSpacing=%f\n' , echoSpacing*1000);
fprintf(fileID,'flipAngle=%f\n' , flipAngle);
fprintf(fileID,'##END');
                    
fclose(fileID);


%% Write image binary

filenameImg = [prefix '_fullrecon'];
writeImgBinary(filenameImg, single(img),'single');



