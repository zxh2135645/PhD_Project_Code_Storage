
[fid_file, fid_path] = uigetfile('*.mat');

load(strcat(fid_path, fid_file));

%load(strcat(fid_path, fid_file), 'Params', 'SpatialCoeff', 'TemporalBasis');

%% What's this struct for?
vec = @(x) x(:);
dispim = @(x) x(:,:,3,:);
resp_phase = 1;
card_phase = 20;
t1_recov = 81;
Gr = TemporalBasis.Gr;
Phi = TemporalBasis.Phi;

L = size(Gr,1);
temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,end,1), L, []);
temp = reshape(reshape(dispim(reshape(SpatialCoeff.U_TV, Params.Ny, Params.Nx, Params.Nz, Params.Necho, [])), [], L)*temp, Params.Ny, Params.Nx, Params.Necho, []);
cw = max(vec(abs(temp(:,:,1,:))));

recon = abs(temp(:,:,1,:)/cw);
ax1 = implay(recon);

%% Old code version 
vec = @(x) x(:);
dispim = @(x) x(:,:,1,:);
resp_phase = 2;
card_phase = 20;
t1_recov = 41;

eco = 1;
L = size(Gr,1);
temp = Gr\reshape(Phi(:,t1_recov,:,resp_phase,end), L, []);
temp = reshape(reshape(dispim(reshape(U_tv, Ny, Nx, Nz, params.NEco, [])), [], L)*temp, Ny, Nx, params.NEco, []);
cw = 0.8*max(vec(abs(temp(:,:,eco,:))));

recon = squeeze(abs(temp(:,:,eco,:)/cw));
recon = fftshift(recon,1);
ax1 = implay(recon);

%% Stack-of-stars
vec = @(x) x(:);
dispim = @(x) x(:,:,13,:);
resp_phase = 1;
card_phase = 20;
t1_recov = 81;
Gr = FitParams.Gr;
Phi = FitParams.Phi;
eco = 5;
L = size(Gr,1);
temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,end), L, []);
temp = reshape(reshape(dispim(reshape(FitParams.U_TV, FitParams.Ny, FitParams.Nx, FitParams.Nz, FitParams.Necho, [])), [], L)*temp, FitParams.Ny, FitParams.Nx, FitParams.Necho, []);
cw = 0.4*max(vec(abs(temp(:,:,eco,:))));

recon = squeeze(abs(temp(:,:,eco,:)/cw));
ax1 = implay(recon);
%%
%exportReconGif(recon, cat(2, FitParams.filePath, '/T1recovery_Late.gif'), 20);
exportReconGif(recon, cat(2, '/home/zhangx1/Documents/Data/Anzhen_AMI_VTR_Radial_FID166925/meas_MID01137_FID166925_1111BEAT_MT_Unified_20250818_20260124T115107/T1recovery_Late_Eco4.gif'), 20);

%%
% --- Reconstruct ALL slices across cardiac phases (replace dispim() slice-14 shortcut) ---

vec = @(x) x(:);

resp_phase = 3;
t1_recov   = 41;
eco        = 1;

L  = size(Gr,1);

% Pull the Phi block for this T1 recovery + resp phase
% Expected dims here: [L, card_phase, (maybe other dims...)] -> we reshape to [L, Nt]
Phi_blk = squeeze(Phi(:, t1_recov, :, resp_phase, end));     % -> [L x Nt] (typically Nt = cardiac phases)
temp_c  = Gr \ reshape(Phi_blk, L, []);                      % -> [L x Nt]
Nt      = size(temp_c, 2);

% Basis images for this echo: U_tv is [Ny x Nx x Nz x NEco x L] (as in your code)
U5   = reshape(U_tv, Ny, Nx, Nz, params.NEco, []);           % -> [Ny x Nx x Nz x NEco x L]
Ueco = squeeze(U5(:,:,:,eco,:));                             % -> [Ny x Nx x Nz x L]

% Multiply basis by coefficients to get images: [Ny x Nx x Nz x Nt]
temp_img = reshape( reshape(Ueco, [], L) * temp_c, Ny, Nx, Nz, Nt );

% Normalize + optional shift (keep what you had)
cw    = 0.5*max(vec(abs(temp_img)));                             % global scale
recon = abs(temp_img) / max(cw, eps);                        % [0..1]
recon = fftshift(recon, 1);
recon = fftshift(recon, 3);

% --- Export as GIF: each frame = montage of ALL slices at one cardiac phase ---
gif_name  = cat(2, fid_path, fid_file(15:23), '_cardiac_motion_all_slices_L16.gif');
delay_sec = 0.07;                 % adjust playback speed
loopcount = inf;

% Choose a grid for tiling slices
ncol = ceil(sqrt(Nz));
nrow = ceil(Nz / ncol);

for t = 1:Nt
    vol = recon(:,:,:,t);         % [Ny x Nx x Nz]

    % Tile slices into one 2D frame (requires Image Processing Toolbox: imtile)
    frame2d = imtile(vol, 'GridSize', [nrow ncol], 'BorderSize', [2 2]);

    % Convert to indexed image for GIF
    frame2d = im2uint8(mat2gray(frame2d));   % uint8 grayscale
    [A,map] = gray2ind(frame2d, 256);

    if t == 1
        imwrite(A, map, gif_name, 'gif', 'LoopCount', loopcount, 'DelayTime', delay_sec);
    else
        imwrite(A, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
    end
end

fprintf('Saved GIF: %s (Nt=%d frames, Nz=%d slices tiled %dx%d)\n', gif_name, Nt, Nz, nrow, ncol);

% --- Optional: also save ONE GIF PER SLICE (true "through all slices", separated) ---
%{
outdir = 'gif_per_slice';
if ~exist(outdir, 'dir'); mkdir(outdir); end
for z = 1:Nz
    gifz = fullfile(outdir, sprintf('cardiac_motion_slice_%02d.gif', z));
    for t = 1:Nt
        img = recon(:,:,z,t);                 % [Ny x Nx]
        img = im2uint8(mat2gray(img));
        [A,map] = gray2ind(img, 256);

        if t == 1
            imwrite(A, map, gifz, 'gif', 'LoopCount', inf, 'DelayTime', delay_sec);
        else
            imwrite(A, map, gifz, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
        end
    end
end
fprintf('Saved per-slice GIFs in: %s\n', outdir);
%}

%% ===== Export Late Enhancement (t1_recov = 200) for ALL slices =====
% Produces:
%   1) One GIF where each frame is a montage of all slices (Nz) at one cardiac phase
%   2) (Optional) One GIF per slice across cardiac phases

vec = @(x) x(:);

resp_phase = 1;
t1_recov_LE = [31,41,51,61,81,101,121,161,201,221];   % late enhancement recovery index
card_phase = 1;
eco = 1;

L  = size(Gr,1);

% --- Coefficients for this recovery + resp phase ---
Phi_blk_LE = squeeze(Phi(:, t1_recov_LE, card_phase, resp_phase, end));   % -> [L x Nt]
coef_LE    = Gr \ reshape(Phi_blk_LE, L, []);                    % -> [L x Nt]
Nt_LE      = size(coef_LE, 2);

% --- Basis images for this echo ---
U5   = reshape(U_tv, Ny, Nx, Nz, params.NEco, []);               % -> [Ny x Nx x Nz x NEco x L]
Ueco = squeeze(U5(:,:,:,eco,:));                                 % -> [Ny x Nx x Nz x L]

% --- Reconstruct volume over time: [Ny x Nx x Nz x Nt] ---
img_LE = reshape( reshape(Ueco, [], L) * coef_LE, Ny, Nx, Nz, Nt_LE );

% --- Normalize (global) + optional shift (match your pipeline) ---
cw_LE   = max(vec(abs(img_LE)));
reconLE = abs(img_LE) / max(cw_LE, eps);                         % [0..1]
reconLE = fftshift(reconLE, 1);
reconLE = fftshift(reconLE, 3);

% --- Export montage GIF (all slices tiled per cardiac phase) ---
gif_LE_name  = cat(2, fid_path, 'late_enhancement_all_slices.gif');
delay_sec   = 0.50;     % slower is often nicer for LGE
loopcount   = inf;

ncol = ceil(sqrt(Nz));
nrow = ceil(Nz / ncol);

for t = 1:Nt_LE
    vol = reconLE(:,:,:,t);                                      % [Ny x Nx x Nz]
    frame2d = imtile(vol, 'GridSize', [nrow ncol], 'BorderSize', [2 2]);

    frame2d = im2uint8(mat2gray(frame2d));
    [A,map] = gray2ind(frame2d, 256);

    if t == 1
        imwrite(A, map, gif_LE_name, 'gif', 'LoopCount', loopcount, 'DelayTime', delay_sec);
    else
        imwrite(A, map, gif_LE_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
    end
end

fprintf('Saved LGE montage GIF: %s (Nt=%d frames, Nz=%d slices tiled %dx%d)\n', ...
    gif_LE_name, Nt_LE, Nz, nrow, ncol);

% --- Optional: One GIF per slice across cardiac phases ---
%{
outdir = 'late_enhancement_gif_per_slice';
if ~exist(outdir, 'dir'); mkdir(outdir); end

for z = 1:Nz
    gifz = fullfile(outdir, sprintf('late_enhancement_slice_%02d.gif', z));
    for t = 1:Nt_LE
        img = reconLE(:,:,z,t);
        img = im2uint8(mat2gray(img));
        [A,map] = gray2ind(img, 256);

        if t == 1
            imwrite(A, map, gifz, 'gif', 'LoopCount', inf, 'DelayTime', delay_sec);
        else
            imwrite(A, map, gifz, 'gif', 'WriteMode', 'append', 'DelayTime', delay_sec);
        end
    end
end

fprintf('Saved LGE per-slice GIFs in: %s\n', outdir);
%}