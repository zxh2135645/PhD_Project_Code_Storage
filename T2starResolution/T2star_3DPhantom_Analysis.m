clear all; close all;
%%
addpath('../function/create-a-heart-phantom-master/')
addpath('../function/')
load('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/Cylinder_Phantom.mat');
load('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/T2star_Parameters_1024x1024x24_SNR_divided_by_100.mat');

%%
dx = 0.04; % in mm
dz = 1; % mm
res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4]; % in mm
res_through_array = [2, 4, 6, 8]; % in mm
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];

% C = zeros(Nx, Ny, sz(3), length(res_array), length(res_through_array), length(Phantom_shape_cell));
% label_mat = zeros(Nx, Ny, sz(3), length(res_array), length(res_through_array), length(Phantom_shape_cell));
MxyTE_remote_echo = T2star_Parameters.MxyTE_remote_echo;
MxyTE_hemo_echo = T2star_Parameters.MxyTE_hemo_echo;
snr_array = T2star_Parameters.snr_array;

Nx = 1024;
Ny = 1024;
% noise_map = T2star_Parameters.noise_map(:,:,:,1);
% clear T2star_Parameters


for k = 1:1
%for k = 1:length(Phantom_shape_cell)

    transmural = transmural_array(k);
    VObj = Phantom_shape_cell{k};
    sz = size(VObj.t2star);
    Nz = sz(3);
    
    [X, Y] = meshgrid(1:sz(1));
    [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));

    t2star_crop_interp = interp2(X, Y, VObj.t2star(:,:,1), Xq, Yq, 'spline');


    % figure();
    % imagesc(t2star_crop_interp); colormap gray;

    t2star_4d_interp = zeros([size(t2star_crop_interp), sz(3), length(MxyTE_remote_echo)]);      %  
    t2star_4d_interp_mask = zeros([size(t2star_crop_interp), sz(3), length(MxyTE_remote_echo)]); % In consideration for ellipsoid
    for nte = 1:length(MxyTE_remote_echo)
        t2star_3d_interp = repmat(t2star_crop_interp, [1, 1, sz(3)]);
        t2star_3d_interp(t2star_3d_interp >= 0.01375) = 0; % Blood Pool
        t2star_3d_interp(t2star_3d_interp < 0.01375 & t2star_3d_interp >= 0.0075) = MxyTE_remote_echo(nte); % Remote
        t2star_3d_interp(t2star_3d_interp < 0.0075 & t2star_3d_interp >= 0.0025) = MxyTE_hemo_echo(nte);    % Hemo
        t2star_3d_interp(t2star_3d_interp < 0.0025 & t2star_3d_interp > 0) = 0;                             % Air
        t2star_4d_interp(:,:,:,nte) = t2star_3d_interp;     
        t2star_4d_interp_mask(:,:,:,nte) = abs(t2star_3d_interp)>0;
    end
    
    
    for n = length(snr_array):length(snr_array)
    % for n = 1:length(snr_array)
    % Add noise map
        noise_map = T2star_Parameters.noise_map(:,:,:,n);
        t2star_4d_interp_noised = t2star_4d_interp + noise_map;
        
        C = zeros(Nx, Ny, sz(3), length(MxyTE_remote_echo), length(res_array), length(res_through_array));
        % label_mat = zeros(Nx, Ny, sz(3), length(MxyTE_remote_echo), length(res_array), length(res_through_array));
        
   
   
        tic
        figure('Position', [100 100 1600 200]);
        %tiledlayout(1, 8, 'Padding', 'none', 'TileSpacing', 'compact');
        for i = 1:length(res_array)
        %for i = 8:8
            res = res_array(i);
            for j = 1:length(res_through_array)
            %for j = 4:4
                res_through = res_through_array(j);
                for nte = 1:length(MxyTE_remote_echo)
                    
                    [C(:,:,:,nte,i,j), ~] = Func_map_to_bloc_3D(dx, dz, Nx, Nz, res, res_through,  t2star_4d_interp_noised(:,:,:,nte));
                    
                    subplot(1,5,nte);
                    imagesc(abs(C(:,:,1,nte,i,j))); axis off; axis equal; colormap gray; caxis([0 0.2])
                    
                end
            end
        end
        toc
        
        %% See if the Noise level is enough
        figure();
        imagesc(abs(C(:,:,1,nte,i,j))); axis off; axis equal; colormap gray; caxis([0 0.2])
        roi = drawcircle('StripeColor','y');
        bw = createMask(roi);
        
        snr = mean(nonzeros(bw .* abs(C(:,:,1,nte,i,j)))) / std(nonzeros(bw .* abs(C(:,:,1,nte,i,j))));
        
        
        %% T2* fitting
        %C(:,:,:,i,j) .* 
    end
    
end

%%



    
        save(cat(2, 'C:\Users\MRI\Documents\James_Data\T2star_SimulationPhantom\Ellipsoid_Blocked_Transmural', num2str(transmural), '.mat'), 'Phantom_shape_cell', '-v7.3');

%end