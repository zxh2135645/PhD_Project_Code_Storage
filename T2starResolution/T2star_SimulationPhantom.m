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
%% Matrixization (Should be deprecated)
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


%% Is it necessary to simulate SNR effect? (Yes I'm doing that)
% To do that in work stastion; 
% TODO need to iterate tissue width: [0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
% This has become the main body

res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2];

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
width_max = [0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0];
TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';
t2star_hemo = 15;
t2star_remote = 38;
tissue_canvas = struct;
tissue_canvas_cell = cell(10, 1);


for cc = 1:length(tissue_canvas_cell)
    for w = 1:length(width_max)
        t_width = width_max(w);
        disp(['width is: ', num2str(t_width)])
        tic;
        tissue_canvas(w).t_width = t_width;
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
        signal_gre = zeros(Ny, Nx, length(TE_array));
        remote_gre = zeros(Ny, Nx, length(TE_array),1);
        hemo_t2star = ones(Ny, Nx) * t2star_hemo;
        remote_t2star = ones(Ny, Nx) * t2star_remote;
        % sigma = 0.02;

        sigma_array = 0.01:0.01:0.1;
        C_t2star_fit_reshape = zeros(Ny, Nx, length(res_array), length(sigma_array));

        for s = 1:length(sigma_array)
            sigma = sigma_array(s);
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
            %options = fitoptions('Method', 'NonlinearLeastSquares');
            %options.Lower = [0 -10];
            %options.Upper = [1, -0.01];
            for i = 1:Nx*Ny*length(res_array)
                f_t = fit(TE_array, C_gre(i,:)', 'exp1', 'Lower', [0 -10], 'Upper', [1 -0.01]);
                C_t2star_fit(i) = -1/f_t.b;
                %f_remote = fit(TE_array, remote_gre, 'exp1');
            end

            C_t2star_fit_reshape(:,:,:,s) = reshape(C_t2star_fit, Ny, Nx, length(res_array));
        end
        toc;
        tissue_canvas(w).C_t2star_fit_reshape = C_t2star_fit_reshape;
    end
    tissue_canvas_cell{cc} = tissue_canvas;
end

%% Plot
C_t2star_fit_reshape = tissue_canvas(1).C_t2star_fit_reshape;
for s = 1:length(sigma_array)
figure();
for i = 1:length(res_array)
   subplot(3,2,i);
   imagesc(C_t2star_fit_reshape(:,:,i,s));
   caxis([0 50]);
end
end
%% Save as mat
% SimPhantom_04042021.res_array = res_array;
% SimPhantom_04042021.sigma_array = sigma_array;
% SimPhantom_04042021.C_t2star_fit_reshape = C_t2star_fit_reshape;
% save_dir = uigetdir;
% fname = 'SimPhantom_04042021';
% save(cat(2, save_dir, '/', fname), 'SimPhantom_04042021');

% Need to send from remote desktop

% SimPhantom_04202021.res_array = res_array;
% SimPhantom_04202021.sigma_array = sigma_array;
% SimPhantom_04202021.width_max = width_max;
% SimPhantom_04202021.tissue_canvas = tissue_canvas;

SimPhantom_10132021.res_array = res_array;
SimPhantom_10132021.sigma_array = sigma_array;
SimPhantom_10132021.width_max = width_max;
SimPhantom_10132021.tissue_canvas = tissue_canvas;

save_dir = uigetdir;
fname = 'SimPhantom_10132021';
save(cat(2, save_dir, '/', fname), 'SimPhantom_10132021');
%% Load SimPhantom_04202021.mat
clear all;
close all;
addpath('../function/');

base_dir = uigetdir;
% f_to_read = cat(2, base_dir, '/Simulation_Results/Phantom/SimPhantom_04202021.mat');
f_to_read = cat(2, base_dir, '/SimPhantom_04202021.mat');
load(f_to_read);
%% Save images 04202021
res_array = SimPhantom_04202021.res_array;
sigma_array = SimPhantom_04202021.sigma_array;
width_max = SimPhantom_04202021.width_max;
save_dir = cat(2, base_dir, '/img/');
if ~exist(save_dir)
    mkdir(save_dir);
end

for v = 1:length(width_max)
    w = width_max(v);
    C_t2star_fit_reshape = SimPhantom_04202021.tissue_canvas(v).C_t2star_fit_reshape;
for s = 1:length(sigma_array)
    figure();
    for i = 1:length(res_array)
        subplot(3,2,i);
        imagesc(C_t2star_fit_reshape(:,:,i,s));
        caxis([0 50]); axis image; axis off;
        colormap(brewermap([],'RdBu'));
    end
    sigma = sigma_array(s);
    sigma_str = num2str(sigma,'%0.2f');
    fname = cat(2, 'SimPhantom_04202021_Sigma_', sigma_str([1 3 4]), '_', num2str(w), 'mm', '.tif');
    saveas(gcf, cat(2, save_dir, fname));
end
close all;
end
%% SimPhantom_04202021 analysis
dx = 0.1;
v = 1;
SimPhantom_analysis = struct;

for v = 1:length(width_max)
w = width_max(v);
C_t2star_fit_reshape = SimPhantom_04202021.tissue_canvas(v).C_t2star_fit_reshape;
Nx = size(C_t2star_fit_reshape, 2);
Ny = size(C_t2star_fit_reshape, 1);
hemo_mask = zeros(size(squeeze(C_t2star_fit_reshape(:,:,:,1))));
myo_mask = zeros(size(squeeze(C_t2star_fit_reshape(:,:,:,1))));
for i = 1:length(res_array)
    res = res_array(i);
    
    % Nx_hemo = res / dx;
    hemo_w = res/dx;
    hemo_ww = fix(hemo_w/2);
    
    if mod(hemo_ww, 2) == 0
        hemo_mask(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1), i) = ones(Ny, length((Nx/2-hemo_ww):(Nx/2+hemo_ww-1)));
        myo_mask(:,:,i) = ~hemo_mask(:,:,i);
    else
        hemo_mask(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww), i) = ones(Ny, length((Nx/2-hemo_ww):(Nx/2+hemo_ww)));
        myo_mask(:,:,i) = ~hemo_mask(:,:,i);
    end
end

mean_hemo_array = zeros(length(sigma_array), length(res_array));
mean_myo_array = zeros(length(sigma_array), length(res_array));
std_myo_array = zeros(length(sigma_array), length(res_array));
CNR_array = zeros(length(sigma_array), length(res_array));

for s = 1:length(sigma_array)
    hemo_temp = reshape(hemo_mask(:,:,:) .* C_t2star_fit_reshape(:,:,:,s), [], length(res_array));
    myo_temp = reshape(myo_mask(:,:,:) .* C_t2star_fit_reshape(:,:,:,s), [], length(res_array));
    for i = 1:length(res_array)
        mean_hemo_array(s,i) = mean(nonzeros(hemo_temp(:,i)));
        mean_myo_array(s,i) = mean(nonzeros(myo_temp(:,i)));
        std_myo_array(s,i) = std(nonzeros(myo_temp(:,i)));
    end
end

CNR_array = abs(mean_hemo_array - mean_myo_array) ./ std_myo_array;
SNR_array = mean_myo_array ./ std_myo_array;
SimPhantom_analysis(v).t_width = w;
SimPhantom_analysis(v).CNR_array = CNR_array;
SimPhantom_analysis(v).SNR_array = SNR_array;
SimPhantom_analysis(v).res_array = res_array;
SimPhantom_analysis(v).sigma_array = sigma_array;
end
save_dir = cat(2, base_dir, '/Simulation_Results/Phantom/');
fname = 'SimPhantom04202021_analysis.mat';
save(cat(2, save_dir, fname), 'SimPhantom_analysis');

%% Heatmap of CNR (04/20/2021)

save_dir = cat(2, base_dir, '/img/');
for i = 1:length(SimPhantom_analysis)
    CNR_array = SimPhantom_analysis(i).CNR_array;
    SNR_array = SimPhantom_analysis(i).SNR_array;
    w = width_max(i);
    figure();
    subplot(1,2,1);
    imagesc(CNR_array); axis image; axis off;
    colorbar;
    subplot(1,2,2);
    imagesc(SNR_array); axis image; axis off;
    caxis([7 30]);
    %colormap(brewermap([],'*RdYlBu'));
    colormap(brewermap([],'*YlGnBu'));
    colorbar;
    fname = cat(2, 'SimPhantom_04202021_CNRSNR_', num2str(w), 'mm', '.tif');
    saveas(gcf, cat(2, save_dir, fname));
end
%% Save images
% res_array = SimPhantom_04202021.res_array;
% sigma_array = SimPhantom_04202021.sigma_array;
% save_dir = cat(2, base_dir, '/img/');
% for s = 1:length(sigma_array)
%     figure();
%     for i = 1:length(res_array)
%         subplot(3,2,i);
%         imagesc(C_t2star_fit_reshape(:,:,i,s));
%         caxis([0 50]); axis image; axis off;
%         colormap(brewermap([],'RdBu'));
%     end
%     sigma = sigma_array(s);
%     sigma_str = num2str(sigma,'%0.2f');
%     fname = cat(2, 'SimPhantom_04202021_Sigma_', sigma_str([1 3 4]) ,'.tif');
%     %saveas(gcf, cat(2, save_dir, fname));
% end
%% Load SimPhantom_04042021.mat
clear all;
close all;
addpath('../function/');

base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/Simulation_Results/Phantom/SimPhantom_04042021.mat');
load(f_to_read);
%% Save images
res_array = SimPhantom_04042021.res_array;
sigma_array = SimPhantom_04042021.sigma_array;
save_dir = cat(2, base_dir, '/img/Simulation_Phantom/');
for s = 1:length(sigma_array)
    figure();
    for i = 1:length(res_array)
        subplot(3,2,i);
        imagesc(SimPhantom_04042021.C_t2star_fit_reshape(:,:,i,s));
        caxis([0 50]);; axis image; axis off;
        colormap(brewermap([],'RdBu'));
    end
    sigma = sigma_array(s);
    sigma_str = num2str(sigma,'%0.2f');
    fname = cat(2, 'SimPhantom_04222021_Sigma_', sigma_str([1 3 4]) ,'.tif');
    %saveas(gcf, cat(2, save_dir, fname));
end

%% Save images 2023
clear all;
close all;
addpath('../function/');

base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/Simulation_Results/Phantom/SimPhantom_03162023.mat');
load(f_to_read);
res_array = SimPhantom_03162023.res_array;
sigma_array = SimPhantom_03162023.sigma_array;
width_max = SimPhantom_03162023.width_max;
save_dir = cat(2, base_dir, '/img/Simulation_Phantom/');
%%
for s = 1:length(sigma_array)
    figure();
    for i = 1:4
        img = SimPhantom_03162023.tissue_canvas{1}(1).C_cell{1,s};
        subplot(2,2,i);
        imagesc(img(:,:,i));
        caxis([0 50]);; axis image; axis off;
        colormap(brewermap([],'RdBu'));
    end
    sigma = sigma_array(s);
    sigma_str = num2str(sigma,'%0.2f');
    fname = cat(2, 'SimPhantom_03162023_w6_Sigma_', sigma_str([1 3 4]) ,'.tif');
    saveas(gcf, cat(2, save_dir, fname));
end
%% SimPhantom_04042021 analysis
dx = 0.1;
Nx = size(SimPhantom_04042021.C_t2star_fit_reshape, 2);
Ny = size(SimPhantom_04042021.C_t2star_fit_reshape, 1);
hemo_mask = zeros(size(squeeze(SimPhantom_04042021.C_t2star_fit_reshape(:,:,:,1))));
myo_mask = zeros(size(squeeze(SimPhantom_04042021.C_t2star_fit_reshape(:,:,:,1))));
for i = 1:length(res_array)
    res = res_array(i);
    
    % Nx_hemo = res / dx;
    hemo_w = res/dx;
    hemo_ww = fix(hemo_w/2);
    
    if mod(hemo_ww, 2) == 0
        hemo_mask(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1), i) = ones(Ny, length((Nx/2-hemo_ww):(Nx/2+hemo_ww-1)));
        myo_mask(:,:,i) = ~hemo_mask(:,:,i);
    else
        hemo_mask(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww), i) = ones(Ny, length((Nx/2-hemo_ww):(Nx/2+hemo_ww)));
        myo_mask(:,:,i) = ~hemo_mask(:,:,i);
    end
end

mean_hemo_array = zeros(length(sigma_array), length(res_array));
mean_myo_array = zeros(length(sigma_array), length(res_array));
std_myo_array = zeros(length(sigma_array), length(res_array));
CNR_array = zeros(length(sigma_array), length(res_array));

for s = 1:length(sigma_array)
    hemo_temp = reshape(hemo_mask(:,:,:) .* SimPhantom_04042021.C_t2star_fit_reshape(:,:,:,s), [], length(res_array));
    myo_temp = reshape(myo_mask(:,:,:) .* SimPhantom_04042021.C_t2star_fit_reshape(:,:,:,s), [], length(res_array));
    for i = 1:length(res_array)
        mean_hemo_array(s,i) = mean(nonzeros(hemo_temp(:,i)));
        mean_myo_array(s,i) = mean(nonzeros(myo_temp(:,i)));
        std_myo_array(s,i) = std(nonzeros(myo_temp(:,i)));
    end
end

CNR_array = abs(mean_hemo_array - mean_myo_array) ./ std_myo_array;
SNR_array = mean_myo_array ./ std_myo_array;
SimPhantom_analysis.CNR_array = CNR_array;
SimPhantom_analysis.SNR_array = SNR_array;
SimPhantom_analysis.res_array = res_array;
SimPhantom_analysis.sigma_array = sigma_array;
save_dir = cat(2, base_dir, '/Simulation_Results/Phantom/');
fname = 'SimPhantom_analysis.mat';
%save(cat(2, save_dir, fname), 'SimPhantom_analysis');
%% Heatmap of CNR
figure();
subplot(1,2,1);
imagesc(CNR_array); axis image; axis off;
colorbar;
subplot(1,2,2);
imagesc(SNR_array); axis image; axis off;
caxis([7 30]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

%% Pure noise map
addpath('../function/');

Ny = 10;
Nx = 10;
sigma = 0.05;
noise_map = randn(Ny, Nx) * sigma;

figure(); imagesc(noise_map);
axis image; axis off;
colormap(brewermap([],'*RdYlBu'));

%% hemo/remote + noise
Ny = 11;
Nx = 11;
gray1 = [233 234 235]/255;
gray2 = [214 216 217]/255;
gray3 = [195 198 199]/255;
gray4 = [176 178 181]/255;
gray5 = [157 159 162]/255;

hemo1 = [214 216 217]/255;
hemo2 = [195 198 199]/255;
hemo3 = [176 178 181]/255;
hemo4 = [139 141 144]/255;
hemo5 = [99 100 102]/255;

% 1
hemo_map = zeros(Ny, Nx, 3);
hemo_map(:,:,1) = ones(Ny,Nx).*gray1(1);
hemo_map(:,:,2) = ones(Ny,Nx).*gray1(2);
hemo_map(:,:,3) = ones(Ny,Nx).*gray1(3);
hemo_map(:,6,1) = ones(Ny,1).*hemo1(1);
hemo_map(:,6,2) = ones(Ny,1).*hemo1(2);
hemo_map(:,6,3) = ones(Ny,1).*hemo1(3);

noise_map = randn(Ny, Nx) * sigma;
hemo_graymap = rgb2gray(hemo_map);
figure(); imagesc(hemo_graymap+noise_map);
axis image; axis off; caxis([0 1]);
colormap gray;

% 2
hemo_map = zeros(Ny, Nx, 3);
hemo_map(:,:,1) = ones(Ny,Nx).*gray2(1);
hemo_map(:,:,2) = ones(Ny,Nx).*gray2(2);
hemo_map(:,:,3) = ones(Ny,Nx).*gray2(3);
hemo_map(:,6,1) = ones(Ny,1).*hemo2(1);
hemo_map(:,6,2) = ones(Ny,1).*hemo2(2);
hemo_map(:,6,3) = ones(Ny,1).*hemo2(3);
noise_map = randn(Ny, Nx) * sigma;
hemo_graymap = rgb2gray(hemo_map);
figure(); imagesc(hemo_graymap+noise_map);
axis image; axis off; caxis([0 1]);
colormap gray;

% 3
hemo_map = zeros(Ny, Nx, 3);
hemo_map(:,:,1) = ones(Ny,Nx).*gray3(1);
hemo_map(:,:,2) = ones(Ny,Nx).*gray3(2);
hemo_map(:,:,3) = ones(Ny,Nx).*gray3(3);
hemo_map(:,6,1) = ones(Ny,1).*hemo3(1);
hemo_map(:,6,2) = ones(Ny,1).*hemo3(2);
hemo_map(:,6,3) = ones(Ny,1).*hemo3(3);
noise_map = randn(Ny, Nx) * sigma;
hemo_graymap = rgb2gray(hemo_map);
figure(); imagesc(hemo_graymap+noise_map);
axis image; axis off; caxis([0 1]);
colormap gray;

% 4
hemo_map = zeros(Ny, Nx, 3);
hemo_map(:,:,1) = ones(Ny,Nx).*gray4(1);
hemo_map(:,:,2) = ones(Ny,Nx).*gray4(2);
hemo_map(:,:,3) = ones(Ny,Nx).*gray4(3);
hemo_map(:,6,1) = ones(Ny,1).*hemo4(1);
hemo_map(:,6,2) = ones(Ny,1).*hemo4(2);
hemo_map(:,6,3) = ones(Ny,1).*hemo4(3);
noise_map = randn(Ny, Nx) * sigma;
hemo_graymap = rgb2gray(hemo_map);
figure(); imagesc(hemo_graymap+noise_map);
axis image; axis off; caxis([0 1]);
colormap gray;

% 5
hemo_map = zeros(Ny, Nx, 3);
hemo_map(:,:,1) = ones(Ny,Nx).*gray5(1);
hemo_map(:,:,2) = ones(Ny,Nx).*gray5(2);
hemo_map(:,:,3) = ones(Ny,Nx).*gray5(3);
hemo_map(:,6,1) = ones(Ny,1).*hemo5(1);
hemo_map(:,6,2) = ones(Ny,1).*hemo5(2);
hemo_map(:,6,3) = ones(Ny,1).*hemo5(3);
noise_map = randn(Ny, Nx) * sigma;
hemo_graymap = rgb2gray(hemo_map);
figure(); imagesc(hemo_graymap+noise_map);
axis image; axis off; caxis([0 1]);
colormap gray;
%% Directly partial voluming on T2* values (needs to be deprecated too)
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

