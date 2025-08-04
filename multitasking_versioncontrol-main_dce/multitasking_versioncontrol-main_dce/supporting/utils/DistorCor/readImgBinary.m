function img = readImgBinary(filename,imgSize,precision)

if nargin < 3
    precision = 'single';
end

fileID = fopen(filename,'rb');
img = fread(fileID,imgSize, precision);
fclose(fileID);
