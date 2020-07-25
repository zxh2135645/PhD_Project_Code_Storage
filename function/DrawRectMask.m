function [RectMask] = DrawRectMask(img_size)
% Since drawrectangle only works for 2018b or newer version
% It might not hurt to just create one 

    [x, y] = ginput(2);
    widt = abs(round(x(2) - x(1)));
    leng = abs(round(y(2) - y(1)));
    RectPts = rectangle('Position', [round(x(1)) round(y(1)) widt, leng], 'FaceColor',[0 .5 .5]);
    
    RectMask = zeros(img_size);
    RectMask(RectPts.Position(2):(RectPts.Position(2)+RectPts.Position(4)), ...
        RectPts.Position(1):(RectPts.Position(1)+RectPts.Position(3))) = 1;
end