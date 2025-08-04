%% 
% [fid_file, fid_path] = uigetfile('*.mat');
% load(strcat(fid_path, fid_file), 'TemporalBasis', 'ReconOptions', 'SpatialCoeff', 'Params');
%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,1,:),1);
vec = @(x) x(:);
Phi = TemporalBasis.Phi;
Ny = Params.Ny;
Nx = Params.Nx;
Nz = Params.Nz;
L = ReconOptions.L;
U = SpatialCoeff.U;
Gr = TemporalBasis.Gr;

temp = Gr\reshape(Phi(:,:,1,1,1), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Params.Necho);
cw = max(vec(abs(temp)))*0.5;

%temp = temp(192-48+1:192+48,:,:);
ax1 = implay(abs(temp/cw));

%% Try T1 mapping

% figure(); imagesc(abs(temp(:,:,300))); axis image;
% %%
% figure(); plot(20:size(Phi,2), squeeze(abs(temp(207, 102,20:end))));
% %%
% figure(); plot(1:size(Phi,5), squeeze(abs(temp(207, 102, :))));