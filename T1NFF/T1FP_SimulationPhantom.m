clear all;
close all;
%%

% Self gridding
res_array = [0.4, 0.6, 0.8, 1.2, 1.6, 2];
half_height = 5;
width_max = [0.3, 0.6, 1, 2];
for w = 1:length(width_max)
    width = width_max(w);
    figure();
    for i = 1:length(res_array)
        subplot(3,2,i);
        rectangle('Position', [-width/2, -half_height, width, half_height*2], 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', [0.8500, 0.3250, 0.0980]);
        set(gca,'Color', [0.25, 0.25, 0.25]);
        hold on;
        
        res = res_array(i);
        Func_Plot_Grids(res, half_height);
        
        axis equal;
        title(sprintf('%1.1f x %1.1f mm^2', res, res))
        %xlabel('mm'); ylabel('mm');
        set(gca, 'FontSize', 12);
        axis([-5 5 -5 5]);
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

%%
res_array = [0.4, 0.6, 0.8, 1.2, 1.6, 2];

T2_myo = 50;
T2_hemo = 50;

fs = 3;
gamma = 42.57;

main_inhom = 0.05;
gxy = 0.3;
T2star_hemo = zeros(1, length(res_array));
T2star_myo = zeros(1, length(res_array));
width_max = [0.3, 0.6, 1, 2];

for i = 1:length(res_array)
    %alpha = 0.3 / res_array(i); 
    T2star_hemo(i) = 1/(1/(T2_hemo/1000) + fs*gamma*(main_inhom +gxy*res_array(i)/res_array(1))) * 1000;
    T2star_myo(i) = 1/(1/(T2_myo/1000) + fs*gamma*(main_inhom)) * 1000;
end

T2star_hybrid = zeros(length(width_max), length(res_array));
for i = 1:length(res_array)
    for w = 1:length(width_max)
        width = width_max(w);
        alpha = width / res_array(i);
        if alpha > 1
            if mod(floor(alpha), 2) == 0
                alpha_real = (width - (floor(alpha)-1)*res_array(i)) / 2;
            else
                alpha_real = (alpha - floor(alpha))/2;
            end
        elseif alpha == 1
            alpha_real = 0;
        else
            alpha_real = alpha;
        end
        
        T2star_hybrid(w, i) = alpha_real * T2star_hemo(i) + (1-alpha_real) * T2star_myo(i);
    end
end

%% TODO Plot
for w = 1:length(width_max)
%w = 4;
width = width_max(w);


% rectangle('Position', [-width/2, -half_height, width, half_height*2], 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', [0.8500, 0.3250, 0.0980]);
T2star_lb = T2star_hemo(end);
T2star_ub = T2star_myo(end);
colorbar_ub = [0.8500, 0.3250, 0.0980];
colorbar_lb = [0.25, 0.25, 0.25];

k = (colorbar_lb - colorbar_ub) ./ (T2star_ub - T2star_lb);

figure();
for i = 1:length(res_array)
    
    grid_width = res_array(i);
    f = floor(width / grid_width);
    color_map = k * (T2star_hybrid(w, i) - T2star_lb) + colorbar_ub;
    inner_color = k * (T2star_hemo(i) - T2star_lb) + colorbar_ub;
    subplot(3,2,i);
    if f == 0
        rectangle('Position', [-grid_width/2 - f*grid_width, -half_height, grid_width, half_height*2], 'FaceColor', color_map, 'EdgeColor', color_map);
    else
        if mod(f,2) == 0
            f = f - 1;
        end
        for ii = 1:f
            iii = floor(ii/2);
            rectangle('Position', [-grid_width/2 + (-1)^ii * iii * grid_width, -half_height, grid_width, half_height*2], 'FaceColor', inner_color, 'EdgeColor', inner_color);
            % rectangle('Position', [-grid_width/2 + (ii-1)*grid_width, -half_height, grid_width, half_height*2], 'FaceColor', inner_color, 'EdgeColor', inner_color);
        end
        iif = (f + 1)/2;
        rectangle('Position', [-grid_width/2 - iif*grid_width, -half_height, grid_width, half_height*2], 'FaceColor', color_map, 'EdgeColor', color_map);
        rectangle('Position', [-grid_width/2 + iif*grid_width, -half_height, grid_width, half_height*2], 'FaceColor', color_map, 'EdgeColor', color_map);
    end
    set(gca,'Color', [0.25, 0.25, 0.25]);
    hold on;
    
    res = res_array(i);
    Func_Plot_Grids(res, half_height);
    
    axis equal;
    title(sprintf('%1.1f x %1.1f mm^2', res, res))
    set(gca, 'FontSize', 12);
    axis([-5 5 -5 5]);
end
end
%%
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

