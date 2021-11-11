clear all;
close all;

%% load data
% Inputs: Nseg,TR,alpha_deg,alpha0_deg,TI

[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');
TE_array = [1.41, 3.38, 5.39, 7.40, 9.41, 11.42]; % 13.43, 15.44 ms
% IR_array = repmat(TE_array, [sizes(2), 1]) + repmat((0:1:(sizes(2)-1)) * params.lEchoSpacing, [sizes(5), 1]).'*1000;
IR_array = repmat(TE_array, [sizes(2), 1]) + repmat((0:1:(sizes(2)-1)) * params.lEchoSpacing, [1, 1]).'*1000;

save_f = cat(2, fid_path, 'mask.mat');
if ~exist(save_f, 'file')
    mask = zeros(Ny, Nx, Nz);
    for i = 1:Nz
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
        for j = 1:1
            for k = 1:1
                temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                cw = 0.5*max(vec(abs(temp)));
                figure();
                %imagesc(abs(temp(:,:,end,1))/cw); axis image;
                % roi = drawpolygon;
                % mask(:,:,i) = createMask(roi);
                mask(:,:,i) = roipoly(abs(temp(:,:,end,1))/cw); axis image;
            end
        end
    end
    
    close all;
    save(save_f, 'mask');
else
    load(save_f);
end

%% Dictionary generation
Nseg = 192;
TR = 0.0131; % 0.013 s?
alpha_deg = 5;
alpha0_deg = 180;
TI = 0.0105;

alpha = alpha_deg*pi/180;
e = @(R1)exp(-TR*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);
n = 1:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);
S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);

R1s = 1./logspace(log10(.1),log10(2),401);
alphas = (.5:.5:(alpha_deg*1.5))*pi/180;
% alpha = alpha_deg * pi / 180;
%   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
minB = -1
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
maxB = -0.5
Bs = linspace(minB,maxB,21);

curves=zeros(Nseg,numel(R1s),numel(alphas),numel(Bs));

for j=1:numel(R1s)
    %for k=1:numel(alphas)
        for l=1:numel(Bs)
            curves(:,j,k,l) = S(1,R1s(j),alpha,Bs(l));
        end
    %end
end

%% A new dictionary
E1 = @(t,R1) exp(-t1*R1);
M0 = 1;
M00 = -M0;
reps = 20;
M01 = M0 * (1-2*E1(TI,R1));
Mz = zeros(Nseg*reps, 1);

for r = 1:reps
   M01 = M00*E1(TI,R1) + M0 * (1-E1(TI,R1));
   for i = 1:Nseg
      Mz(i + (r-1)*Nseg) = M01 * cos(alpha) * E1(TR,R1) + M0 * (1 - E1(TR,R1));
      M01 = Mz(i+(r-1)*Nseg);
   end
   M00 = -M01;
end
%%
vec = @(x) x(:);
t1_map_3d = zeros(Ny, Nx, Nz);
t1_map = zeros(Ny, Nx, Nz, 1, 1, sizes(5));
N_nt = 3;
t1_map_2d_nt = zeros(Ny, Nx, N_nt);

%% for i = 1:Nz
%nt = 3;
for nt = 1:N_nt
    for i = 3:3
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
        for j = 1:1
            for k = 1:1
                temp = Gr\reshape(Phi(:,:,j,k,nt), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                % temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], nt);
                for m = 1:1
                    ipt = abs(reshape(temp(:,:,:,m), Ny, Nx, 1, []));
                    ipt_2d = reshape(ipt, [], Nseg);
                    mask_1d = vec(mask(:,:,i));
                    
                    tic;
                    T1Mapping_DictFit_Func;
                    toc;
                end
            end
        end
        % t1_map_3d(:,:,i) = t1_map_2d;
        t1_map_2d_nt(:,:,nt) = t1_map_2d;
    end
end
%% Save maps
save_fname = cat(2, fid_path, fid_file(1:end-4), '_Dict_T1Mapping_slc3_L1.mat');
save(save_fname, 't1_map_2d_nt');

%% Plot
figure();
for i = 1:N_nt
    subplot(2,2,i);
    imagesc(t1_map_2d_nt(:,:,i)); axis image; caxis([100 1000]);
    axis off;
    title(cat(2, 'Section ', num2str(i)));
    colormap gray;
end
%% Draw ROIs 
figure('Position', [100 100 1600 1600]);
BW = zeros(size(t1_map_2d_nt));
for i = 1:N_nt
    imagesc(t1_map_2d_nt(:,:,i)); axis image; caxis([100 1000]);
    axis off;
    title(cat(2, 'Section ', num2str(i)));
    
    BW(:,:,i) = roipoly;
end

% T1 values of ROI
t1_mean = zeros(size(BW, 3),1);
t1_sd = zeros(size(BW, 3), 1);
for i = 1:N_nt
   t1_mean(i) = mean(nonzeros(BW(:,:,i).*t1_map_2d_nt(:,:,i))); 
   t1_sd(i) = std(nonzeros(BW(:,:,i).*t1_map_2d_nt(:,:,i)));
end

metrics.t1_mean = t1_mean;
metrics.t1_sd = t1_sd;
metrics.BW = BW;
save_fname = cat(2, fid_path, fid_file(1:end-4), '_T1mapping_Analysis_L1.mat');
save(save_fname, '-struct', 'metrics');