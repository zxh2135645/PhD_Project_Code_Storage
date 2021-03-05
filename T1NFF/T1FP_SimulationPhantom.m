clear all;
close all;
%%

% Self gridding
res_array = [0.4, 0.6, 0.8, 1.2, 1.6, 2];
half_height = 5;
width_max = 1;
for width = 1:width_max
    figure();
    for i = 1:length(res_array)
        subplot(3,2,i);
        rectangle('Position', [-width/2, -half_height, width, half_height*2], 'FaceColor', [0.8500, 0.3250, 0.0980]);
        set(gca,'Color', [0.25, 0.25, 0.25]);
        hold on;
        
        res = res_array(i);
        Plot_Grids(res, half_height);
        
        axis equal;
        title(sprintf('Resolution: %1.1f x %1.1f mm^2', res, res))
        axis([-10 10 -5 5]);
    end
end


% hold on; 
% rectangle('Position', [-1, -5, 2, 10]);
% axis([-10 10 -5 5]);
% 
% rectangle('Position', [-1.5, -5, 3, 10]);
% axis([-10 10 -5 5]);
% 
% rectangle('Position', [-2, -5, 4, 10]);
% axis([-10 10 -5 5]);
% 
% rectangle('Position', [-2.5, -5, 5, 10]);
% axis([-10 10 -5 5]);

% 
T2_myo = 50;
T2_hemo = 50;
T2star_hemo = 20;
T2star_myo = 35;
%%
function Plot_Grids(res, half_height)
lim = floor(half_height / res);
for i = 1:lim
    xline(res*(2*i-1)/2, 'Color', [1,1,1]);
    yline(res*(2*i-1)/2, 'Color', [1,1,1]);
    xline(-res*(2*i-1)/2, 'Color', [1,1,1]);
    yline(-res*(2*i-1)/2, 'Color', [1,1,1]);
end
end

% imageExample
function imageExample()
img = imread('peppers.png');
sz = size(img);
hFig = figure();
hAx = axes();
image([1 sz(2)]+100, [1 sz(1)]+200, img)
set(hFig, 'WindowButtonDownFcn', @mouseDown)

    function mouseDown(o, e)
        p = get(hAx, 'CurrentPoint');
        p = p(1, 1:2);
        
        x = round(axes2pix(sz(2), [1 sz(2)], p(1)));
        y = round(axes2pix(sz(1), [1 sz(1)], p(2)));
        title(sprintf('image pixel = (%d, %d)', x, y))
    end
end

