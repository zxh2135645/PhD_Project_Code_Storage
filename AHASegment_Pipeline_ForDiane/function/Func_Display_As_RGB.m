function [RGB_out] = Func_Display_As_RGB(img, R, G, B)
        
    img = img - min(img(:));
    img = img ./ max(img(:));
    img = im2uint8(img);
    map = gray(256);
    Gray = ind2rgb(img, map);
    RGB_out = Gray;
    
if nargin < 2 || isempty(R)
    R = zeros(size(Gray, 1), size(Gray, 2));
end
if nargin < 3 || isempty(G)
    G = zeros(size(Gray, 1), size(Gray, 2));
end
if nargin < 4 || isempty(B)
    B = zeros(size(Gray, 1), size(Gray, 2));
end

    RGB_out(:,:,1) = Gray(:,:,1) + R;
    RGB_out(:,:,2) = Gray(:,:,2) + G;
    RGB_out(:,:,3) = Gray(:,:,3) + B;

end