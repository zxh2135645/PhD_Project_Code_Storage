clear all;
close all;

save_dir = uigetdir;
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
% width_max = [0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0];
width_max = [0.6, 2.0];
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
            res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2];
            C_cell = cell(length(res_array), 1);
            C = zeros([Ny, Nx, size(t,3)]);

            for i = 1:length(res_array)
                res = res_array(i);
                for te = 1:length(TE_array)
                    C(:,:,te) = Func_map_to_bloc_Exaustive(dx, Nx, res, t(:,:,te));
                end
            end
            %

            C_gre = reshape(C, [], length(TE_array));
            C_t2star_fit = zeros(size(C_gre, 1), 1);
            %options = fitoptions('Method', 'NonlinearLeastSquares');
            %options.Lower = [0 -10];
            %options.Upper = [1, -0.01];
            parfor i = 1:Nx*Ny*length(res_array)
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

SimPhantom_02282023.res_array = res_array;
SimPhantom_02282023.sigma_array = sigma_array;
SimPhantom_02282023.width_max = width_max;
SimPhantom_02282023.tissue_canvas = tissue_canvas_cell;

fname = 'SimPhantom_02282023';
save(cat(2, save_dir, '/', fname), 'SimPhantom_02282023');


%% Extaustive Sliding window 03/10/2023
% Renew the hemo_mask and and myo_mask
%% Load SimPhantom_04042021.mat
addpath('../function/');

base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_02282023.mat');
load(f_to_read);

res_array = SimPhantom_02282023.res_array;
sigma_array = SimPhantom_02282023.sigma_array;
width_max = SimPhantom_02282023.width_max;

C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(1).C_t2star_fit_reshape;
C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(3).C_t2star_fit_reshape;
% Width x Height x Res x Sigma
figure();
imagesc(C_t2star_fit_reshape(:,:,1,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,7,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,1,10))
figure();
imagesc(C_t2star_fit_reshape(:,:,7,10))
%%
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(1).C_t2star_fit_reshape;
Nx = size(C_t2star_fit_reshape, 2);
Ny = size(C_t2star_fit_reshape, 1);
myo_mask = zeros(size(squeeze(C_t2star_fit_reshape(:,:,:,1))));

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww), :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape,3));
myo_mask(:, (Nx/2+3*hemo_ww)+1:end, :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape, 3));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask_cell = cell(length(res_array),1);
%thresh_mat = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
hemo_mask_frame_cell = cell(length(res_array),1);

for iter = 1:length(SimPhantom_02282023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{iter}(w).C_t2star_fit_reshape;

        % Nx_hemo = res / dx;
        hemo_w = round(t_width/dx);
        hemo_ww = fix(hemo_w/2);

        if mod(hemo_w, 2) == 0
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1),w) = 1;
        else
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww),w) = 1;
        end

        for i = 1:length(res_array)
            res = res_array(i);
            for j = 1:length(sigma_array)
                blocs = Func_map_to_bloc_encoded2(dx, Nx, res, C_t2star_fit_reshape(:,:,i,j));
                hemo_mask_frame = zeros(Ny, Nx, length(res_array), length(width_max), size(blocs,3));
                hemo_mask = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas), size(blocs,3));
                for xx = 1:size(blocs, 3)
                    blocs_overlap = t_hemo(:,:,w) .* blocs(:,:,xx);
                    idx_blocs = nonzeros(unique(blocs_overlap));
                    for k = 1:length(idx_blocs)
                        [row, col] = ind2sub([Ny Nx], find(blocs == idx_blocs(k)));
                        hemo_mask_frame(row, col, i, w, xx) = 1;
                    end
                    thresh = mean(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j))) - 2*std(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j)));
                    hemo_mask(:,:,i,j,w,iter,xx) = (C_t2star_fit_reshape(:,:,i,j) < thresh); %.* hemo_mask_frame(:,:,i,w);
                    %thresh_mat(:,:,i,j,w) = thresh;
                end
                hemo_mask_cell{xx} = hemo_mask;
                hemo_mask_frame_cell{xx} = hemo_mask_frame;
            end
        end
    end
end