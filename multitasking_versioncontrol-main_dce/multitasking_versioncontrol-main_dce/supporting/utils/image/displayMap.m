function map_color = displayMap(map,cmap,range,fps,fileName,filePath)

if nargin == 5
    filePath = [];
    if isempty(fileName)
        fileName = 'map.gif';
    end
end

if nargin < 4
    fps = 4;
end

if nargin < 3
    range = [min(map(:)) max(map(:))];
end
range = range - range(1);
map   = (map - range(1))/range(2);

if nargin < 2
    cmap = parula(256);
end

 %% play cardiac cycle
map_color = zeros(size(repmat(map(:,:,1,:),[1 1 3 1])));

for n = 1:size(map_color,4)
    map_color(:,:,:,n) = ind2rgb(round(255*(map(:,:,1,n))),cmap);
end
if size(map_color,4) == 1
    map_color(:,:,:,2) = map_color;
end
implayZoom(map_color,fps);


%% save colored gif
if nargin > 3
    try
        filename1 = [filePath '_' fileName];
        temp = uint8(round(255*(map(:,:,1,:))));
        imwrite(temp(:,:,:,1),cmap,filename1,'gif','LoopCount',inf','DelayTime',1/fps);
        for n = 2:size(temp,4)
            imwrite(temp(:,:,:,n),cmap,filename1,'gif','WriteMode','append','delaytime',1/fps);
        end
    catch
        fprintf("Saving images to %s failed. Saving them in current folder.", filePath);
        temp = uint8(round(255*(map(:,:,1,:))));
        imwrite(temp(:,:,:,1),cmap,filename,'gif','LoopCount',inf','DelayTime',1/fps);
        for n = 2:size(temp,4)
            imwrite(temp(:,:,:,n),cmap,filename,'gif','WriteMode','append','delaytime',1/fps);
        end
    end
end
