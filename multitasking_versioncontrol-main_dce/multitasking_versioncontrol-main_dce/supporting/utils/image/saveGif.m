function saveGif(img,filePath,fileName,fps,cmap)

if nargin < 5
    cmap = gray(256);
end

if nargin < 4
    fps = 20;
end

if nargin < 3
    fileName = '_temp.gif';
end

if nargin < 2 || isempty(filePath)
    filePath = '.';
end

if ~strcmp(fileName(end-3:end),'.gif')
    fileName = [fileName '.gif'];
end

imwrite(uint8(255*abs(reshape(img,size(img,1),size(img,2),1,[]))),cmap,[filePath '/' fileName],'gif','loopcount',inf','delaytime',1/fps)

