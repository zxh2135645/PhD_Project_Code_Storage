clear all;
close all;
%% Load BinningResult
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi_rt_full_init', 'L', 'U_init', 'Ny', 'Nx', 'Nz', 'vec','params', 'Phicard', 'Phiresp', 'ccL', 'Norig');
%% 1) Before binning - realtime play
dispim = @(x,st)fftshift(x(:,:,3,:),1);
L_init = 32;
recon = reshape(U_init, Ny, Nx, Nz, L_init);
recon = dispim(recon);
recon = reshape(recon, [], L_init);

recon = reshape(recon*Phi_rt_full_init(:, 3:1:400), Ny, Nx, []);

cw = max(recon(:));
implay(abs(recon)/abs(cw));

%% 2) After binning
Norig = 192;
ccL = 1;
resp = 4;
temp = reshape(reshape(dispim(reshape(U_init, Ny, Nx, Nz, [])),[],L_init)*reshape(reshape(Phiresp, [], 1) ,[], L_init).', Ny, Norig, []);

cw = max(temp(:));
implay(abs(temp)/abs(cw));

% cardiac phase
temp2 = reshape(reshape(dispim(reshape(U_init, Ny, Nx, Nz, [])),[],L_init)*reshape(reshape(permute(Phicard(resp,:,:,:), [2,1,3,4]), [], ccL) ,[], L_init).', Ny, Nx, []);

cw2 = max(temp2(:));
implay(abs(temp2)/abs(cw2));

%% (optional) plot each cardiac phase
% figure();
% for i = 1:size(Phicard, 2)
%     subplot(5,5,i);
%     imagesc(abs(temp2(:,:,i))/abs(cw2)); axis image; axis off; colormap gray;
% end
%
wx = 80;
wy = 80;
X = size(temp2,1);
Y = size(temp2,2);
centroid = [X/2, Y/2];
delay_time = 1;

save_path = cat(2, fid_path, 'CardiacBinning_Debug/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

fh = figure();
for i = 1:size(Phicard, 2)
    
    temp_crop = imcrop(abs(temp2(:,:,i)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    temp_crop = flip(imrotate(temp_crop, 90),2);
    
    imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    title(cat(2, 'cbin = ', num2str(i)));
    % pause(1);
    
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_LRT_CINE_cbin', num2str(size(Phicard, 2)), '.gif');
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', delay_time, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', delay_time, 'WriteMode','append');
    end
end

%% Load Reconstructed Image
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');

%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,3,:),1);

temp = Gr\reshape(Phi(:,41,:,resp,end), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));


ax1 = implay(abs(temp/cw));
% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;

%% Optional
fh = figure();
for i = 1:size(Phicard, 2)
    
    temp_crop = imcrop(abs(temp(:,:,i)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    temp_crop = flip(imrotate(temp_crop, 90),2);
    
    imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    title(cat(2, 'cbin = ', num2str(i)));
    % pause(1);
    
    drawnow
    % Capture the plot as an image
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), '_LRT_CINE_cbin', num2str(size(Phicard, 2)), '.gif');
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', delay_time, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', delay_time, 'WriteMode','append');
    end
end