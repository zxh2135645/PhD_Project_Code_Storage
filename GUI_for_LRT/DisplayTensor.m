clear all; 
close all;
%% 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
% dispim = @(x)fftshift(x(:,:,SliceSlider.Value,:),1);
temp = Gr\reshape(Phi(:,end,:,1,:), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = 0.5*max(vec(abs(temp)));

ax1 = implay(abs(temp(:,:,:,1)/cw));
