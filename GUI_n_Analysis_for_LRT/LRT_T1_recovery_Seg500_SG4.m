clear all;
close all;

%% T1 recovery of the current LRT parameters settings
% Updated in July 11, 2022
% Dictionary fitting looks more promising to me


% Nseg = 680;
% TR = 0.00357;
Nseg = 500;
TR = 0.00357;
alpha_deg = 5;
alpha0_deg = 180; % For IR
% alpha0_deg = 90; % For SR
TI = 0.0105;

alpha = 5 * pi / 180;
e = @(R1) exp(-TR*R1);
Mss = @(e,alpha) (1-e) / (1-cos(alpha)*e); % Steady State
n = 1:Nseg;
Sint = @(A,e,alpha,B) A * Mss(e,alpha) * (1 - (B-1)*(e*cos(alpha)).^(n-1));
S = @(A,R1,alpha,B) Sint(A, e(R1), alpha, B);

R1 = 1 ./ logspace(log10(1),log10(1),1);

B0 = -1;
% curve = Sint(1, e(R1), alpha, B0);

%% New fitting
cutoff = 41;


E1 = @(t, R1) exp(-t*R1);
M0= 1;
M00 = -M0;
reps = 20; % Why is this 20?? to reach a pseudo-steady-state
M01 = M0 * (1 - 2*E1(TI, R1));
Mz = zeros(Nseg*reps, 1);
R1s = 1 ./ logspace(log10(.1),log10(2),401); % For IR
%R1s = 1 ./ logspace(log10(.1),log10(4),401); % For SR

Mz_dict = zeros(length(R1s), Nseg);
Mz_dict_norm = zeros(length(R1s), Nseg);

for k = 1:length(R1s)
    R1 = R1s(k);
    for r = 1:reps
        M01 = M00 * E1(TI,R1) + M0 * (1 - E1(TI,R1));
        for i = 1:Nseg
            Mz(i + (r-1)*Nseg) = M01 * cos(alpha) * E1(TR, R1)  + M0 * (1 - E1(TR, R1)); % This is it
            M01 = Mz(i+(r-1)*Nseg);
        end
        M00 = -M01;
    end
    Mz_dict(k,:) = Mz(end-Nseg+1:end);
    Mz_dict_norm(k,:) = Mz_dict(k,:) ./ max(Mz_dict(k,:));
end

%figure();
%plot(Mz)

Mz_dict_norm_abs = abs(Mz_dict_norm);
Mz_dict_norm_abs_truc = Mz_dict_norm_abs(:,(cutoff+1):end);
%Mz_dict_norm_abs_truc = Mz_dict_norm_abs(:,(cutoff+1):2:end);
%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx');

%%
% Load mask
mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Ny, Nx, Nz);
if ~exist(mask_f)
    for i = 1:Nz
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
        for j = 1:1
            for k = 1:1
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
%%
% cardiac phase and resp phase needs to be encoded
% Dave_D8      [12, 4]



%% IR
N_nt = 15;
N_nt = 8;
N_nt = size(Phi, 5);
slc = 3;
t1_map_3d_nt = zeros(Ny, Nx, Nz, N_nt);
card_phase_array = [15, 16, 17, 18, 19, 20];
%card_phase = 9;
resp_phase = 4;
t1_map_4d_nt = zeros(Ny, Nx, Nz, N_nt, length(card_phase_array));


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
%t1_map_3d_nt_shifted = fftshift(t1_map_3d_nt, 3);
t1_map_3d_nt_shifted = fftshift(t1_map_4d_nt, 3);

figure();
for slc = 1:size(t1_map_3d_nt_shifted, 3)
    subplot(4,4,slc);
    imagesc(sum(t1_map_3d_nt_shifted(:,:,slc,2),5));
    clim([0 2000]);
end