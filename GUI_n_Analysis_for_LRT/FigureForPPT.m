%% Load Data
clear all;
close all;

[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');

Nseg = size(Phi,2); Nsect = size(Phi,5);
%% Montage of T1 evolution
%slc = 3;
%dispim = @(x)fftshift(x(:,:,slc,:),1);


img_save = cat(2, fid_path, 'img_montage/');
if ~exist(img_save, 'dir')
   mkdir(img_save); 
end
%%  Draw Heart masks and compute centroid
soi = [6, 5, 4, 3, 2, 1, 14, 13, 12, 11];
% soi = [5, 4, 3, 2, 1, 16, 15, 14, 13, 12];
t1oi = [21, 41, 61, 81, 101, 121, 141, 161, 181];
wx = 64;
wy = 64;
mask = zeros(Ny, Nx, Nz);
figure('Position',[100, 100, 720, 800]);
% tight_subplot(Nh, Nw, gap, marg_h, marg_w)
for i = 1:length(soi)
    slc = soi(i);
    dispim = @(x) fftshift(x(:,:,slc,:), 1);
    for j = 2:2
        % Seg#41
        t1_idx = t1oi(j);
        temp = Gr\reshape(Phi(:,t1_idx,1,end,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
        cw = max(vec(abs(temp)));
        imagesc(abs(temp(:,:,1)/cw)); axis image; 
        colormap gray; axis off;
        
        roi = drawpolygon;
        mask(:,:,slc) = createMask(roi);
    end
end

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    slc = soi(i);
    subplot(4,4,i)
    imagesc(mask(:,:,slc))
    s = regionprops(mask(:,:,slc),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end
%% Plot Montage of T1 evolution
odd_array = (1:Nseg/2) * 2 - 1;
odd_array = [11, 21, 31, 41, 61, 81, 101, 121, 141, 161, 181];
depth = length(odd_array);
caxis_cell = {[0.3 0.9],[0.2 0.75]};
Nsect_array = [1, 15];
%Nsect_array = 1:15;
figure('Position',[0, 100, 1800, 800]);
[ha, pos] = tight_subplot(2,1,[.01 -.1],[.01 .01],[.01 .01]);

temp_crop_4d = zeros(wy, wx, depth, Nsect);
for n = 1:length(Nsect_array)
    nsec = Nsect_array(n);
    for slc = 3:3
            % slc = soi(i);
            dispim = @(x) fftshift(x(:,:,slc,:), 1);
            temp = Gr\reshape(Phi(:,:,1,end,:), L, []);
            temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
            
            %phase_temp = angle(temp);
            %phase_diff = phase_temp - phase_temp(:,:,end);
            %cos_phase_diff = cos(phase_diff);
            %cos_phase_diff(cos_phase_diff >= 0) = 1;
            %cos_phase_diff(cos_phase_diff <  0) = -1;
            
            idx = find(soi == slc);
            centroid = centroids{idx};
            
            %temp_crop = imcrop3(abs(temp(:,:,odd_array,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            temp_crop = imcrop3(abs(temp(:,:,odd_array,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (depth-1)]);
            %temp_psir = abs(temp(:,:,odd_array,nsec)) .* cos_phase_diff(:,:,odd_array,nsec);
            %psir = (temp_psir - min(temp_psir(:))) ./ (max(temp_psir(:)) - min(temp_psir(:)));
            %psir_crop = imcrop3(psir, [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            % cos_phase_diff_crop = imcrop3(cos_phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            %phase_diff_crop = imcrop3(phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            cw1 = 0.8*max(vec(abs(temp(:,:,odd_array,nsec))));
            %cw2 = 0.5*max(vec(phase_diff(:,:,odd_array,nsec)));
            %cw3 = 0.5*max(vec(psir));
            
            
            axes(ha(n));
            temp_crop = imrotate(flip(temp_crop, 1), 90);
            montage(temp_crop./cw1, 'Size', [1 11]);


            %axes(ha(2));
            %montage(phase_diff_crop./cw2, 'Size', [12 8]);
            %caxis([-pi pi]);
            %axes(ha(3));
            %montage(psir_crop, 'Size', [12 8]);
            %caxis(caxis_cell{n});
            temp_crop_4d(:,:,:,n) = temp_crop;
    end
end

% cw = max(vec(abs(temp_crop_4d(:,:,2:end,n))));
% axes(ha(1));
% %axes(ha(n));
% temp_crop_4d = imrotate(flip(temp_crop_4d, 1), 90);
% temp_crop_reshape = reshape(temp_crop_4d, wy, wx, []);
% montage(temp_crop_reshape./cw, 'Size', [length(Nsect_array), 11]);

fname = cat(2, fid_file(1:17), fid_file(end-25:end-20), '_Nseg1-15.png');
saveas(gcf,cat(2, img_save, fname));
close all;

%% Make gifs for cardiac, respiratory, T1, T2* evolution
% T1 evolution
wx = 80;
wy = 80;
slc = 3;

idx = find(soi == 3);
dispim = @(x)fftshift(x(:,:,slc,:),1);
odd_array = 11:5:141;

temp = Gr\reshape(Phi(:,odd_array,1,4,15), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));
temp = imrotate(permute(temp, [2,1,3]), 180);

% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
ax2 = implay(abs(temp(:,:,:)/cw));


save_path = cat(2, fid_path, 'cine/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

fh = figure('Position', [100 100 300 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% RR = slice_data(end).TriggerTime + slice_data(end).RepetitionTime; % ms
% time_resolution_cine = RR / length(slice_data) / 1000;
temproal_resolution = 0.05;
for n = 1:size(temp,3)
    centroid = centroids{idx};
    % Draw plot
    temp_crop = imcrop(abs(temp(:,:,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);

    imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_CMR_T1Evo', '_WholeHeart2.gif');
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'WriteMode','append');
    end
end


%% Cardiac Motion
idx = find(soi == 3);
dispim = @(x)fftshift(x(:,:,slc,:),1);
t1_idx = 41;

temp = Gr\reshape(Phi(:,t1_idx,:,4,15), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));
temp = imrotate(permute(temp, [2,1,3]), 180);

% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
ax2 = implay(abs(temp(:,:,:)/cw));


save_path = cat(2, fid_path, 'cine/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

fh = figure('Position', [100 100 300 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% RR = slice_data(end).TriggerTime + slice_data(end).RepetitionTime; % ms
% time_resolution_cine = RR / length(slice_data) / 1000;
temproal_resolution = 0.05;
for n = 1:size(temp,3)
    centroid = centroids{idx};
    % Draw plot
    temp_crop = imcrop(abs(temp(:,:,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);

    imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_CMR_CardiacMotion', '_WholeHeart.gif');
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'WriteMode','append');
    end
end

%% Respiratory Motion
idx = find(soi == 3);
dispim = @(x)fftshift(x(:,:,slc,:),1);
t1_idx = 61;

temp = Gr\reshape(Phi(:,t1_idx,2,:,15), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));
temp = imrotate(permute(temp, [2,1,3]), 180);

% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
ax2 = implay(abs(temp(:,:,:)/cw));


save_path = cat(2, fid_path, 'cine/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

fh = figure('Position', [100 100 300 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% RR = slice_data(end).TriggerTime + slice_data(end).RepetitionTime; % ms
% time_resolution_cine = RR / length(slice_data) / 1000;
temproal_resolution = 0.75;
for n = 1:size(temp,3)
    centroid = centroids{idx};
    % Draw plot
    temp_crop = imcrop(abs(temp(:,:,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);

    imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_CMR_RespiratoryMotion', '_WholeHeart.gif');
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'WriteMode','append');
    end
end

%% Show mixture of multiple dimensions
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes', 'Phi_rt_full');
%%
slc = 3;
idx = find(soi == slc);
dispim = @(x)fftshift(x(:,:,slc,:),1);
ds = round(1/(params.lEchoSpacing)); %downsampling to get to 20fps
recon=reshape(U,Ny,Nx,Nz,L);
recon=dispim(recon); %dispim uses 1st partition
recon=reshape(recon,[],L);

recon=(reshape(recon*Phi_rt_full(:,13:6:1000),Ny,Nx,[])); % XZ

cw=max(recon(:));
implay(abs(recon)/abs(cw));
%%
fh = figure('Position', [100 100 300 400]);
axis tight manual
wx = 80;
wy = 80;

mask = zeros(Ny, Nx, Nz);
temproal_resolution = params.lEchoSpacing;
save_path = cat(2, fid_path, 'cine/');

for n = 1:size(recon,3)
    centroid = centroids{idx};
    % Draw plot
    recon_rotate = imrotate(flip(recon, 1), 90);
    recon_crop = imcrop(abs(recon_rotate(:,:,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);

    imagesc(abs(recon_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_CMR_MixtureMotion', '_Cropped.gif');

    
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', temproal_resolution, 'WriteMode','append');
    end
end