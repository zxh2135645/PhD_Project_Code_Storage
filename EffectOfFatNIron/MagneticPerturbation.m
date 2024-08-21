%%
clear all; close all;

%%
addpath('../function/')

B0_dir = [0 0 1]';

matrix_size = [2048 2048 24];
voxel_size = [2, 2, 1];

domain = 'kspace';
D=dipole_kernel(matrix_size, voxel_size, B0_dir, domain);


%% 

figure();
imshow3D(abs(fftshift(fftshift(fftshift(D,1),2),3)))

figure();
imshow3D(abs(D))

%% Susceptibility map
load('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/Cylinder_Phantom.mat');
suscept_hemo = 110.2;
suscept_remote = 1.36;
suscept_blood_lv = 1.36;


t2star = Phantom_shape_cell{2}.t2star;
%VObj = Phantom_shape_cell{k};
%sz = size(VObj.t2star);
sz = size(t2star);

[X, Y] = meshgrid(1:sz(1));
[Xq, Yq] = meshgrid(linspace(1, sz(1), 2048));

t2star_crop_interp = interp2(X, Y, t2star(:,:,1), Xq, Yq, 'nearest');
t2star_3d_interp = repmat(t2star_crop_interp, [1, 1, sz(3)]);

%%
susc_map = zeros(size(t2star_3d_interp));

%unique_values = unique(t2star);
susc_map(t2star_3d_interp >= 0.0175) = suscept_remote;
susc_map(t2star_3d_interp < 0.01375 & t2star_3d_interp >= 0.0075) = suscept_remote;
susc_map(t2star_3d_interp < 0.0075 & t2star_3d_interp >= 0.0025) = suscept_hemo;
susc_map(t2star_3d_interp < 0.0025) = 0;

susc_map = susc_map * 10^-9;
% susc_map(t2star == 0.005) = suscept_hemo;
% susc_map(t2star == 0.01) = suscept_remote;
% susc_map(t2star == 0.0175) = suscept_remote;

susc_map_air = ones(size(t2star_3d_interp)) * 9 * 10^-9;;

figure(); imshow3D(susc_map);

%%
perturb_map = fftn(fftshift(fftn(susc_map)) .* D) + fftn(fftshift(fftn(susc_map_air)) .* D);
%%
figure(); imagesc(abs(perturb_map(:,:,1)));

