clear all;
close all;

%% T1 recovery of the current LRT parameters settings
% And This is set up for VTR version
% Updated in July 11, 2022
% Dictionary fitting looks more promising to me

%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file));
%%
addpath('../function/')
params = FitParams.params;
reconOptions = FitParams.reconOptions;
[curves, T1s] = genCurveSubspace_for_ParamFit(params, FitParams, reconOptions);

cutoff = 10;
Mz_dict_norm_abs = abs(squeeze(curves(:,:,1,1,end,end))) ./ max(max(abs(squeeze(curves(:,:,1,1,end,end)))));
Mz_dict_norm_abs_truc = Mz_dict_norm_abs((cutoff+1):end,:);
%Mz_dict_norm_abs_truc = Mz_dict_norm_abs(:,(cutoff+1):2:end);

%% This is trying to fit from reconstructed images
%% Should be fitting to FitParams
Ny = FitParams.Ny;
Nx = FitParams.Nx;
Nz = FitParams.Nz;
Phi = FitParams.Phi;
Gr = FitParams.Gr;
L = size(Phi, 1);
U = FitParams.U_TV;
Neco = params.Necho;

Nydisp = params.Nydisp;
Nxdisp = params.Nxdisp;

%%
% Load mask
vec = @(x) x(:);
mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Nydisp, Nxdisp, Nz);
if ~exist(mask_f)
    for i = 1:Nz
        dispim = @(x,st) x(:,:,i,:,:);
        for j = 1:1
            for k = 1:1
                temp = Gr\reshape(Phi(:,:,j,k,end), L, []);
                temp = reshape(reshape(dispim(U),[],L) * temp, Nydisp, Nxdisp, Neco, []);
                cw = 0.5*max(vec(abs(temp(:,:,1,:))));
                figure();
                imagesc(abs(temp(:,:,1,end))/cw); axis image;
                roi = drawpolygon;
                mask(:,:,i) = createMask(roi);
            end
        end
    end
    save(mask_f, 'mask');
else
    load(mask_f);
end

 close all;

 %%
 dispim = @(x) x(:,:,11,:,:);
 resp_phase = 1;
 card_phase = 1;
 temp = Gr\reshape(Phi(:,41,:,resp_phase,end), L, []);
 temp = reshape(reshape(dispim(U),[],L)*temp, Nydisp, Nxdisp, Neco, []);
 cw = 0.8*max(vec(abs(temp(:,:,1,:))));


 ax1 = implay(abs(temp(:,:,1,:)/cw));
%%
% cardiac phase and resp phase needs to be encoded
% Dave_D8      [12, 4]



%% IR
N_nt = 8;
N_nt = size(Phi, 5);
slc = 3;
t1_map_3d_nt = zeros(Nydisp, Nxdisp, Nz, N_nt);
card_phase_array = [1];
resp_phase = 1;
t1_map_4d_nt = zeros(Nydisp, Nxdisp, Nz, N_nt, length(card_phase_array));
Nseg = size(Phi,2);

for nt = 1:N_nt
    %for nt = 1:2
    for i = 1:Nz
        %for i = slc:slc
        dispim = @(x,st) x(:,:,i,:,:);
        for j = 1:length(card_phase_array)
            card_phase = card_phase_array(j);

            temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,nt), L, []);
            temp = reshape(reshape(dispim(U),[],L)*temp, Nydisp, Nxdisp, Neco, []);
            cw = 0.5*max(vec(abs(temp(:,:,1,:))));

            % ipt_2d = abs(reshape(temp(:,:,11:end), [], (Nseg-20)/2));
            ipt_2d = abs(reshape(temp(:,:,1,(cutoff+1):end), [], (Nseg-cutoff)));
            % ipt_2d = abs(reshape(temp(:,:,:), [], Nseg));
            % mask = roipoly(abs(temp(:,:,41)) / cw); axis image;
            mask_1d = vec(mask(:,:,i));


            tic;
            % T1Mapping_DictFit_Func2_Anzhen;
            t1_map_2d = T1Mapping_DictFit_Func2_Anzhen(ipt_2d, mask_1d, Mz_dict_norm_abs_truc, Nydisp, Nxdisp, T1s);
            toc;

            %t1_map_3d_nt(:,:,i,nt) = t1_map_2d;
            t1_map_4d_nt(:,:,i,nt,j) = t1_map_2d;
        end
    end
end

%% Display T1 maps
%t1_map_3d_nt_shifted = fftshift(t1_map_3d_nt, 3);
t1_map_3d_nt_shifted = t1_map_4d_nt;

figure();
for slc = 1:size(t1_map_3d_nt_shifted, 3)
    subplot(4,5,slc);
    imagesc(sum(t1_map_3d_nt_shifted(:,:,slc,end),5));
    clim([0 800]);
end

%%
figure();
slc = 11

for i = 1:10
    imagesc(sum(t1_map_3d_nt_shifted(:,:,slc,i),5));
    clim([0 800]);
    colormap gray
    pause(.5)
end

%%
figure();
for i = 1:10
    subplot(3,4,i)
    imagesc(sum(t1_map_3d_nt_shifted(:,:,slc,i),5));
    clim([0 650]);
    colormap gray
    title(cat(2, num2str(i), ' min'))
    colorbar;
    axis image;
    axis off;
    %pause(.5)
end

%%
figure(); plot(squeeze(t1_map_3d_nt_shifted(77,105,slc,:)));
hold on; plot(squeeze(t1_map_3d_nt_shifted(92,95,slc,:)));


