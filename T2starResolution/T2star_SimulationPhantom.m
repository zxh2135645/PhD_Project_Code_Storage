clear all;
close all;
%%
addpath('../function/')
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
    %T2star_hemo(i) = 1/(1/(T2_hemo/1000) + fs*gamma*(main_inhom +gxy*res_array(i)/res_array(1))) * 1000;
    T2star_myo(i) = 1/(1/(T2_myo/1000) + fs*gamma*(main_inhom)) * 1000;
end

T2star_hemo = repmat(15, [1,6]); % From iron phantoms, resolution is not playing a major role



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

% Output is T2star_hybrid;
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
%% Matrixization
dx = 0.1; % mm
dy = 0.1; % mm
X = 10;
Y = 10;

Nx = X / dx;
Ny = Y / dy;
figure();
sigma = 2;

for w = 1:length(res_array)
    width = res_array(w);
    t = randn(Ny, Nx) * sigma^2 + T2star_myo(w);
    hemo_w = res_array(w)/dx;
    hemo_w = fix(hemo_w/2);
    if mod(hemo_w, 2) == 0
        t(:, (Nx/2-hemo_w):(Nx/2+hemo_w-1)) = randn(Ny, length((Nx/2-hemo_w):(Nx/2+hemo_w-1))) * sigma + T2star_hybrid(1,w);
    else
        t(:, (Nx/2-hemo_w):(Nx/2+hemo_w)) = randn(Ny, length((Nx/2-hemo_w):(Nx/2+hemo_w))) * sigma + T2star_hybrid(1,w);
    end
    
    
    subplot(3,2,w);
    imagesc(t);
    axis image; axis off;
    colormap(brewermap([],'RdBu'));
    caxis([0 50]);
    
    title(sprintf('%1.1f x %1.1f mm^2', width, width))
    set(gca, 'FontSize', 12);
    %axis([-5 5 -5 5]);
end


%% Is it necessary to simulate SNR effect?
dx = 0.1; % mm
dy = 0.1; % mm
X = 10;
Y = 10;

Nx = X / dx;
Ny = Y / dy;
gre = @(rho, alpha, TR, T1, TE, T2star) rho*sin(alpha)*(1-exp(-TR./T1))./(1-cos(alpha)*exp(-TR/T1))*exp(-TE./T2star);
d2r = @(x) x/180*pi;
vec = @(x) x(:);

% Width of the tissue
width_max = [0.3, 0.6, 1, 2];
t_width = width_max(1);
hemo_w = round(t_width / dx);
hemo_ww = fix(hemo_w/2);
t_hemo = zeros(Ny, Nx);
t_remote = zeros(Ny, Nx);

if mod(hemo_w, 2) == 0
    t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1)) = 1;
    t_remote = ~t_hemo;
else
    t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww)) = 1;
    t_remote = ~t_hemo;
end

%gre(1, d2r(60), 790, 1200, 40, 15)
TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';
signal_gre = zeros(Ny, Nx, length(TE_array));
remote_gre = zeros(Ny, Nx, length(TE_array),1);
sigma = 0.02;
hemo_t2star = ones(Ny, Nx) * 15;
remote_t2star = ones(Ny, Nx) * 38;

for i = 1:length(TE_array)
    TE = TE_array(i);
    noise = randn(Ny, Nx) * sigma; 
    signal_gre(:,:,i) = gre(1, d2r(18), 121, 1200, TE, hemo_t2star) + noise;
    remote_gre(:,:,i) = gre(1, d2r(18), 121, 1200, TE, remote_t2star) + noise;
end

t = signal_gre .* t_hemo + remote_gre .* t_remote;
res_array = [0.4, 0.6, 0.8, 1.2, 1.6, 2];
C = zeros([Ny, Nx, length(res_array), size(t,3)]);

for i = 1:length(res_array)
    res = res_array(i);
    for te = 1:length(TE_array)
        C(:,:,i,te) = Func_map_to_bloc(dx, Nx, res, t(:,:,te));
    end
end
%

C_gre = reshape(C, [], length(TE_array));
C_t2star_fit = zeros(size(C_gre, 1), 1);
options = fitoptions('Method', 'NonlinearLeastSquares');
options.Lower = [0 -10];
options.Upper = [1, -0.01];
for i = 1:Nx*Ny*length(res_array)
    f_t = fit(TE_array, C_gre(i,:)', 'exp1', options);
    C_t2star_fit(i) = -1/f_t.b;
    %f_remote = fit(TE_array, remote_gre, 'exp1');
end

C_t2star_fit_reshape = reshape(C_t2star_fit, Ny, Nx, length(res_array));

%% Plot

%% Directly partial voluming on T2* values
t_gre = reshape(t, [], length(TE_array));
t_t2star_fit = zeros(size(t_gre, 1), 1);
options = fitoptions('Method', 'NonlinearLeastSquares');
options.Lower = [0 -10];
options.Upper = [1, -0.01];
for i = 1:Nx*Ny
    f_t = fit(TE_array, t_gre(i,:)', 'exp1', options);
    t_t2star_fit(i) = -1/f_t.b;
    %f_remote = fit(TE_array, remote_gre, 'exp1');
end

t_t2star_fit_reshape = reshape(t_t2star_fit, Ny, Nx);

% Matrices
figure();
for i = 1:length(res_array)
    
res = res_array(i);
C = Func_map_to_bloc(dx, Nx, res, t_t2star_fit_reshape);
subplot(3,2,i); imagesc(C); caxis([0 50]);
colormap(brewermap([],'RdYlBu')); axis image; axis off;

end



%% DEPRECATED
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

