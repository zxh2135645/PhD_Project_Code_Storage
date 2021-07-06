function [RectMask] = DrawRectMask_V2_GUI2020(loc, img_size, ax)
% Since drawrectangle only works for 2018b or newer version
% It might not hurt to just create one
% RectMask is 2D
    x = loc(1,:);
    y = loc(2,:);
    
    widt = abs(round(x(2) - x(1)));
    leng = abs(round(y(2) - y(1)));
    if x(1) <= x(2) && y(1) <= y(2)
        x_start = x(1);
        y_start = y(1);
    elseif x(1) > x(2) && y(1) > y(2)
        x_start = x(2);
        y_start = y(2);
    elseif x(1) <= x(2) && y(1) > y(2)
        x_start = x(1);
        y_start = y(2);
    elseif x(1) > x(2) && y(1) <= y(2)
        x_start = x(2);
        y_start = y(1);
    end
    
    RectPts = rectangle('Position', [round(x_start) round(y_start) widt, leng], 'EdgeColor',[0.8500, 0.3250, 0.0980], ...
        'LineWidth', 1.5, 'Parent', ax);
    RectMask = zeros(img_size(1), img_size(2));
    % te = [RectPts.Position(2):(RectPts.Position(2)+RectPts.Position(4)), ...
    %  RectPts.Position(1):(RectPts.Position(1)+RectPts.Position(3))];
    
    if x(1) > x(2) && y(1) > y(2)
        % reverse
        RectMask(round(y(2)):round(y(1)), round(x(2)):round(x(1))) = 1;
    elseif x(1) <= x(2) && y(1) <= y(2)
        % This is where default goes
        RectMask(round(y(1)):round(y(2)), round(x(1)):round(x(2))) = 1;
    elseif x(1) > x(2) && y(1) <= y(2)
        RectMask(round(y(1)):round(y(2)), round(x(2)):round(x(1))) = 1;
    elseif x(1) <= x(2) && y(1) > y(2)
        RectMask(round(y(2)):round(y(1)), round(x(1)):round(x(2))) = 1;
    end
    % RectMask(RectPts.Position(2):(RectPts.Position(2)+RectPts.Position(4)), ...
    %   RectPts.Position(1):(RectPts.Position(1)+RectPts.Position(3))) = 1;
end