clear all; 
close all;
%% 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,3,:),1);

temp = Gr\reshape(Phi(:,:,1,1), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
%cw = 0.5*max(vec(abs(temp)));


% ax1 = implay(abs(temp/cw));
%figure();
%ax2 = imagesc(abs(temp(:,:,180)/cw)); axis image; colormap gray;axis off;
seg = 110;
card = 3;
temp_3D = zeros(Ny, Nx, Nz);
for i = 1:Nz
    dispim = @(x) fftshift(x(:,:,i,:), 1);
    temp = Gr\reshape(Phi(:,:,card,end), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    temp_3D(:,:,i) = abs(temp(:,:,seg));
end

%% Save temp_3D as tif  (Optional)
addpath('../function/RY/')
save_dir = cat(2, fid_path, 'Card', num2str(card), '_Seg', num2str(seg), '/');
mat2tif(temp_3D, save_dir);

%% crop image
x_array = zeros(Nz, 1);
y_array = zeros(Nz, 1);
figure();
wx = 96;
wy = 96;
temp_3D_crop = zeros(wy, wx, Nz);
for i = 1:size(temp_3D, 3)
    if i == 1
        imagesc(temp_3D(:,:,i)); axis off; axis image;
        [x, y] = getpts;
    end
    temp_3D_crop(:,:,i) = imcrop(temp_3D(:,:,i), [x-wx/2, y-wy/2, wx-1, wy-1]);
end

%% montage
figure();
cw = 0.5*max(vec(temp_3D_crop));
montage(fftshift(temp_3D_crop/cw,3)); axis image; colormap gray;

%% volume show
fullViewPnl = uipanel(figure,'Title','Original Volume');
volshow(squeeze(fftshift(temp_3D_crop/cw,3)),'Parent',fullViewPnl);

%% Cardiac motion
dispim = @(x) fftshift(x(:,:,i,:), 1);
temp = Gr\reshape(Phi(:,:,1,end), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = 0.5*max(vec(abs(temp)));
ax1 = implay(abs(temp/cw));



