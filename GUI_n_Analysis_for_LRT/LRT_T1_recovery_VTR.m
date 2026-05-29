%% T1 recovery of the current LRT parameters settings
% Updated in July 11, 2022
% Dictionary fitting looks more promising to me
TR1 = 0.0035;
TR2 = 0.0114;
Nseg = 192;
TR = TR1 + TR2;
alpha_deg = 5;
alpha0_deg = 180; % For IR
TI = 0.0105;


alpha = 5 * pi / 180;
e1 = @(R1) exp(-TR1*R1);
e2 = @(R1) exp(-TR2*R1);

Mss = @(e,alpha) (1-e) / (1-cos(alpha)*e); % Steady State
n = 1:Nseg;
Sint = @(A,e,alpha,B) A * Mss(e1,alpha) * (1 - (B-1)*(e1*cos(alpha)).^(n-1));
S = @(A,R1,alpha,B) Sint(A, e(R1), alpha, B);
R1 = 1 ./ logspace(log10(1),log10(1),1);
B0 = -1;
%%
[params, FitParams, dictOptions] = defaultT1DictionaryVTRInputs();
dict = genT1Dictionary_VTR(params, FitParams, dictOptions);

T1s = dict.T1s;
Mz_dict_norm_abs_truc = dict.Mz_dict_norm_abs_truc.';

figure(); plot(Mz_dict_norm_abs_truc(41,:))

%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'U_tv');
resp = 2;
card = 1;
% Load mask
mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Ny, Nx, Nz);
if ~exist(mask_f)
    for i = 1:Nz
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
        for j = card:card
            for k = resp:resp
                temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                cw = 0.5*max(vec(abs(temp)));
                figure();
                imagesc(abs(temp(:,:,end,1))/cw); axis image;
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
 %% IR
% N_nt = 15;
% N_nt = 8;
N_nt = size(Phi, 5);
% slc = 5;
t1_map_3d_nt = zeros(Ny, Nx, Nz, N_nt);
card_phase_array = [24];
resp_phase = 2;
t1_map_4d_nt = zeros(Ny, Nx, Nz, N_nt, length(card_phase_array));
cutoff = dictOptions.cutoff;
Nseg = params.linesPerShot
params.NEco = 1;
R1s = 1 ./ dictOptions.T1s;
for nt = 1:N_nt
    %for nt = 1:2
    for i = 1:Nz
        %for i = slc:slc
        dispim = @(x,st) fftshift(x(:,:,i,:), 1);
        for j = 1:length(card_phase_array)
            card_phase = card_phase_array(j);

            temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,nt), L, []);
            temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
            cw = 0.5*max(vec(abs(temp)));

            % ipt_2d = abs(reshape(temp(:,:,11:end), [], (Nseg-20)/2));
            ipt_2d = abs(reshape(temp(:,:,(cutoff+1):end), [], (Nseg-cutoff)));
            % ipt_2d = abs(reshape(temp(:,:,:), [], Nseg));
            % mask = roipoly(abs(temp(:,:,41)) / cw); axis image;
            mask_1d = vec(mask);

            tic;
            T1Mapping_DictFit_Func2;
            toc;

            %t1_map_3d_nt(:,:,i,nt) = t1_map_2d;
            t1_map_4d_nt(:,:,i,nt,j) = t1_map_2d;
        end
    end
end

%% Display T1 maps
t1_map_3d_nt_shifted = fftshift(t1_map_4d_nt, 3);

figure();
% for slc = 1:size(t1_map_3d_nt_shifted, 3)
for slc = 5:14
    subplot(2,5,slc-4);
    imagesc(t1_map_3d_nt_shifted(:,:,slc,end,1));
    clim([0 1000]);
    axis image;
end

%% Save As MAT
save_dir = cat(2, fid_path, 'dDCE_T1_Dict_Diastole/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%save(cat(2, save_dir, 'T1Map_PostCon_Seg15_June26th_', num2str(card_phase), '.mat'), 't1_map_3d_nt_shifted');
save(cat(2, save_dir, 'T1Map_PostCon_Seg10', '.mat'), 't1_map_3d_nt_shifted');
