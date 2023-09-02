function OutImage = ImPostProc(InputImage)
    % InputImage is a 3D binary image
    % Remove 1 pixel clusters
    OutImage = InputImage;
    for j = 1:size(InputImage, 3)
        [L, n] = bwlabel(InputImage(:,:,j));
        NumOfPixel = 1;
        Temp = OutImage(:,:,j);
        for i = 1:n
            if sum(L(:) == i) <= NumOfPixel
                Temp(L == i) = 0;
            end
        end
        
        % Fill holes of one infarct
        % Temp = imfill(Temp, 'holes');
        OutImage(:,:,j) = Temp;
    end
end