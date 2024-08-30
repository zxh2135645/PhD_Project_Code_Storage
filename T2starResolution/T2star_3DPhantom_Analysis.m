clear all; close all;

addpath('.\function\')
%load('C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Cylinder_Phantom.mat');
load('C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellipsoid_Phantom.mat');
load('C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\T2star_Parameters_1024x1024x24_SNR_divided_by_100.mat');

vec = @(x) x(:);
%%
% assume the heart shape from equatorial plane to the top of apex is 3.6cm
dx = 0.04; % in mm
dz = 0.1; % mm
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

TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';
A = TE_array .^ [0, 1];

% parfor k = 1:1
for k = 1:length(Phantom_shape_cell)

    transmural = transmural_array(k);
    VObj = Phantom_shape_cell{k};
    sz = size(VObj.t2star);
    Nz = sz(3);

    
    [X, Y] = meshgrid(1:sz(1));
    [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));

    t2star_crop_interp = interp2(X, Y, VObj.t2star(:,:,1), Xq, Yq, 'spline');


    % figure();
    % imagesc(t2star_crop_interp); colormap gray;

    %%
    t2star_4d_interp = zeros([size(t2star_crop_interp), sz(3), length(MxyTE_remote_echo)]);      %  
    t2star_4d_interp_mask = zeros([size(t2star_crop_interp), sz(3), length(MxyTE_remote_echo)]); % In consideration for ellipsoid
    

    for nte = 1:length(MxyTE_remote_echo)
        t2star_3d_interp = repmat(t2star_crop_interp, [1, 1, sz(3)]);
        t2star_3d_interp(abs(t2star_3d_interp) >= 0.01375) = 0; % Blood Pool
        t2star_3d_interp(abs(t2star_3d_interp) < 0.01375 & abs(t2star_3d_interp) >= 0.0075) = MxyTE_remote_echo(nte); % Remote
        t2star_3d_interp(abs(t2star_3d_interp) < 0.0075 & abs(t2star_3d_interp) >= 0.0025) = MxyTE_hemo_echo(nte);    % Hemo
        t2star_3d_interp(abs(t2star_3d_interp) < 0.0025 & abs(t2star_3d_interp) > 0) = 0;                             % Air
        t2star_4d_interp(:,:,:,nte) = t2star_3d_interp;     
        t2star_4d_interp_mask(:,:,:,nte) = abs(t2star_3d_interp)>0;
    end
    
    % hemo_mask_gt = repmat(t2star_crop_interp, [1, 1, sz(3)]);
    % %t2star_3d_interp = repmat(t2star_crop_interp, [1, 1, sz(3)]);
    % hemo_mask_gt(abs(hemo_mask_gt) >= 0.01375) = 0; % Blood Pool
    % hemo_mask_gt(abs(hemo_mask_gt) < 0.0025 & abs(hemo_mask_gt) > 0) = 0;
    % hemo_mask_gt(abs(hemo_mask_gt) < 0.01375 & abs(hemo_mask_gt) >= 0.0075) = 2; % Remote
    % hemo_mask_gt(abs(hemo_mask_gt) < 0.0075 & abs(hemo_mask_gt) >= 0.0025) = 3;    % Hemo
    % 
    % se = strel('disk', 2);
    % hemo_mask_gt = imerode(hemo_mask_gt,se);
    
        
    %for n = length(snr_array):length(snr_array)
    for n = 1:length(snr_array)
    % Add noise map
        
        noise_map = T2star_Parameters.noise_map(:,:,:,n);
        t2star_4d_interp_noised = t2star_4d_interp + noise_map;
        
        C = zeros(Nx, Ny, sz(3), length(MxyTE_remote_echo), length(res_array), length(res_through_array));
        label_mat = zeros(Nx, Ny, sz(3), length(res_array), length(res_through_array));
        
        t2star_map = zeros(Nx*Ny*sz(3), length(res_array), length(res_through_array));
   
        %tic
        %figure('Position', [100 100 1600 200]);
        %tiledlayout(1, 8, 'Padding', 'none', 'TileSpacing', 'compact');
        for i = 1:length(res_array)
        %for i = 1:1
            res = res_array(i);
            for j = 1:length(res_through_array)
            %for j = 1:1
                res_through = res_through_array(j);
                for nte = 1:length(MxyTE_remote_echo)

                    %tic;
                    [C(:,:,:,nte,i,j), label_mat(:,:,:,i,j)] = Func_map_to_bloc_3D(dx, dz, Nx, Nz, res, res_through,  t2star_4d_interp_noised(:,:,:,nte));
                    %toc;

                %     if nte == 1
                %         label_mat_masked = label_mat(:,:,:,i,j) .* t2star_4d_interp_mask(:,:,:,nte);
                %         label_array = vec(label_mat_masked);
                %         [unique_blocs, ia, ic] = unique(label_array);
                %         C_gre = zeros(length(unique_blocs), length(MxyTE_remote_echo));
                %         % C_t2star_fit = zeros(length(unique_blocs), 1);
                %     end
                % 
                %     %subplot(1,5,nte);
                %     %imagesc(abs(C(:,:,1,nte,i,j))); axis off; axis equal; colormap gray; caxis([0 0.2])
                % 
                %     C_blocked_array = vec(C(:,:,:,nte,i,j) .* t2star_4d_interp_mask(:,:,:,nte));
                %     C_gre(:, nte) = C_blocked_array(ia);
                % 
                %     %for bloc = 1:(length(unique_blocs) - 1)
                %         % label = unique_blocs(bloc+1); % exclude 0;
                %         % idx = find(label_array == label); % All the weighting in this idx are the same
                %         % C_gre(bloc, nte) = C_blocked_array(idx(1));
                % 
                %         %if nte == length(MxyTE_remote_echo)
                %             %f_t = fit(TE_array, abs(C_gre(bloc,:))', 'exp1', 'Lower', [0 -10], 'Upper', [1 -0.01]);
                %             %C_t2star_fit(bloc) = -1/f_t.b;
                % 
                % 
                %         %end
                %     %end
                end
                % 
                % f_t = A\log(abs(C_gre).');
                % C_t2star_fit = abs(1./f_t(2,:));
                % t2star_map(:,i,j) = C_t2star_fit(ic);               
            end
        end


        % t2star_map = reshape(t2star_map, Nx, Ny, sz(3), length(res_array), length(res_through_array));

        % save(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\T2starMap_Blocked_LinReg_Transmural', num2str(transmural), '_NoiseLevel', num2str(n), '.mat'), 't2star_map', '-v7.3');
        save(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\T2starW_Blocked_Transmural', num2str(transmural), '_NoiseLevel', num2str(n), '.mat'), 'C', '-v7.3');
    end
end
        
        %% See if the Noise level is enough
        % nte = 1;
        % figure('Position', [100 100 1600 200]);
        % for nte = 1:length(MxyTE_remote_echo)
        %     subplot(1,5,nte);
        %     imagesc(abs(C(:,:,1,nte,i,j))); axis off; axis equal; colormap gray; caxis([0 0.2])
        % end

        %roi = drawcircle('StripeColor','y');
        %bw = createMask(roi);
        
        %snr = mean(nonzeros(bw .* abs(C(:,:,1,nte,1,1)))) / std(nonzeros(bw .* abs(C(:,:,1,nte,1,1))))
        
        % figure();
        % imagesc(abs(t2star_map(:,:,1,1,1))); axis off; axis equal; colormap gray; caxis([0 100]);
        % roi = drawcircle('StripeColor','y');
        % bw = createMask(roi);
        % 
        % nte = 1;
        % snr_array = zeros(length(res_array), length(res_through_array));
        % for i = 1:length(res_array)
        % %for i = 8:8
        %     for j = 1:length(res_through_array)
        %         snr_array(i,j) = mean(nonzeros(bw .* abs(C(:,:,1,nte,i,j)))) / std(nonzeros(bw .* abs(C(:,:,1,nte,i,j))));
        %     end
        % end
        % 
        % t2star_map(t2star_map <= 0.198) = 0.198;
        % t2star_map(t2star_map >= 100) = 100;
        % 
        % t2starnr_array = zeros(length(res_array), length(res_through_array));
        % for i = 1:length(res_array)
        % %for i = 8:8
        %     for j = 1:length(res_through_array)
        %         t2starnr_array(i,j) = mean(nonzeros(bw .* abs(t2star_map(:,:,1,i,j)))) / std(nonzeros(bw .* abs(t2star_map(:,:,1,i,j))));
        %     end
        % end

