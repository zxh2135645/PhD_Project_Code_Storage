clear all;
close all;

%% load data
% Inputs: Nseg,TR,alpha_deg,alpha0_deg,TI

[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');
TE_array = [1.41, 3.38, 5.39, 7.40, 9.41, 11.42]; % 13.43, 15.44 ms
IR_array = repmat(TE_array, [sizes(2), 1]) + repmat((0:1:(sizes(2)-1)) * params.lEchoSpacing, [sizes(5), 1]).'*1000;

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

%   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
minB = -1
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
maxB = -0.5
Bs = linspace(minB,maxB,21);

curves=zeros(Nseg,numel(R1s),numel(alphas),numel(Bs));

for j=1:numel(R1s)
    for k=1:numel(alphas)
        for l=1:numel(Bs)
            curves(:,j,k,l) = S(1,R1s(j),alphas(k),Bs(l));
        end
    end
end

%%
vec = @(x) x(:);
t1_map_3d = zeros(Ny, Nx, Nz);
t1_map = zeros(Ny, Nx, Nz, 1, 1, sizes(5));
for i = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    for j = 1:1
        for k = 1:1
            temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
            temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
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
    t1_map_3d(:,:,i) = t1_map_2d;
end




