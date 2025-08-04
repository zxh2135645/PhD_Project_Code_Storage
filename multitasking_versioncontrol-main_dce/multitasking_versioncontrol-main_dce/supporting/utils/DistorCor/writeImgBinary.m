function writeImgBinary(filename,img,precision)

if nargin < 3
    precision = 'single';
end

fileID = fopen(filename,'wb');
fwrite(fileID,img,precision);
fclose(fileID);
