function [PixelMask] = DrawPixelMask_GUI2020(loc, img_size)
% PixelMask is 2D
x = loc(1,:);
y = loc(2,:);

PixelMask = zeros(img_size(1), img_size(2));

for i = 1:numel(x)
    PixelMask(round(y(i)), round(x(i))) = 1;
end
end