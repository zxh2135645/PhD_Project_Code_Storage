%% T1 recovery of the current LRT parameters settings
% Updated in July 11, 2022
% Dictionary fitting looks more promising to me

Nseg = Params.lSegments;
TR = Params.lEchoSpacing ;
alpha_deg = Params.flipAngleArray;
TI = 0.0105;

alpha = alpha_deg * pi / 180;
e = @(R1) exp(-TR*R1);
Mss = @(e,alpha) (1-e) / (1-cos(alpha)*e); % Steady State
n = 1:Nseg;
Sint = @(A,e,alpha,B) A * Mss(e,alpha) * (1 - (B-1)*(e*cos(alpha)).^(n-1));
S = @(A,R1,alpha,B) Sint(A, e(R1), alpha, B);

R1 = 1 ./ logspace(log10(1),log10(1),1);

B0 = -1;
% curve = Sint(1, e(R1), alpha, B0);

%% New fitting
E1 = @(t, R1) exp(-t*R1);
M0= 1;
M00 = -M0;
reps = 100;
M01 = M0 * (1 - 2*E1(TI, R1));
Mz = zeros(Nseg*reps, 1);
R1s = 1 ./ logspace(log10(.1),log10(2000),1001);

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
Mz_dict_norm_abs_truc = Mz_dict_norm_abs(:,101:end-100);
% Load Data
% [fid_file, fid_path] = uigetfile('*.mat');
% load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx');
%%
fid_path=pwd;
% Load mask
mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Ny, Nx, Nz);
if ~exist(mask_f)
    for i = 1
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
        for j = 1:1
            for k = 1:1
                temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], Params.Necho);
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
N_nt = size(Phi,5);
slc = 1;
t1_map_3d_nt = zeros(Ny, Nx, Nz, N_nt);
card_phase =1;
resp_phase = 2;
for nt = 1:N_nt
%for nt = 1:2
    for i = 1
    %for i = slc:slc
   
        dispim = @(x,st) fftshift(x(:,:,i,:), 1);
        
        temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,nt), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Params.Necho);
        cw = 0.5*max(vec(abs(temp)));
        
        ipt_2d = abs(reshape(temp(:,:,101:end-100), [], (Nseg-200)));
        % ipt_2d = abs(reshape(temp(:,:,:), [], Nseg));
        % mask = roipoly(abs(temp(:,:,41)) / cw); axis image;
        mask_1d = vec(mask);
        
        tic;
        T1Mapping_DictFit_Func2;
        toc;
        
        t1_map_3d_nt(:,:,i,nt) = t1_map_2d;
    end
    
end
% Check images
figure();
for i = 1:size(t1_map_3d_nt, 4)
    subplot(2,5,i);
    imagesc(t1_map_3d_nt(:,:,1,i),[0 900]); axis image;colormap("jet");
end


% why!
% figure('Position', [100 100 900 900]);
% % imagesc(t1_map_3d_nt(:,:,3,1));
% mapoi = t1_map_3d_nt(:,:,slc,1);
% imagesc(mapoi); axis image;
% roi_coord = drawpolygon;
% roi = createMask(roi_coord);
% 
% mean_t1_array = zeros(size(t1_map_3d_nt, 4), 1);
% for i = 1:size(t1_map_3d_nt, 4)
%     mean_t1_array(i) = mean(nonzeros(t1_map_3d_nt(:,:,slc,i) .* roi));
% end



% N_nt = 15;
% rep_array = [1:5:192];
% 
% for nt = 1:N_nt
%     for i = slc:slc
%         dispim = @(x,st) fftshift(x(:,:,i,:), 1);
%         
%         temp = Gr\reshape(Phi(:,:,1,1,nt), L, []);
%         temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
%         cw = 0.5*max(vec(abs(temp)));
%         
%         implay(abs(temp(:,:,rep_array)) ./ cw);
%     end
% end

%% Display T1 maps
t1_map_3d_nt_shifted = fftshift(t1_map_3d_nt, 3);
figure();
for slc = 1:size(t1_map_3d_nt_shifted, 3)
    subplot(4,4,slc);
    imagesc(t1_map_3d_nt_shifted(:,:,slc,14));
end
%% Save As MAT
save_dir = cat(2, fid_path, 'dDCE_T1_Dict_Diastole/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save(cat(2, save_dir, 'T1Map_PostCon_Seg15.mat'), 't1_map_3d_nt_shifted');
%% ImageJ
% Convert mat to dicom (T1 map)
figure();
for i = 1:size(t1_map_3d_nt, 4)
    subplot(4,4,i);
    imagesc(t1_map_3d_nt(:,:,slc,i).*mask(:,:,slc)); caxis([100 1000]); axis image;
end

save_dir = cat(2, fid_path, 'DICOM_T1_Dict_Diastole/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

dicom_dir = uigetdir;

sizes = size(Phi);
Func_ConvertMat_Dicom_ImageJ(t1_map_3d_nt, fid_file, dicom_dir, save_dir, slc, sizes);


%% For Mona's Registration purpose: 10/24/2024
save_dir = cat(2, fid_path, 'T1_w_Diastole/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

dispim = @(x,st) fftshift(x(:,:,slc,:), 1);
temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,:), L, []);
temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
%cw = 0.5*max(vec(abs(temp)));
cw = 7e-04;
implay(abs(temp)/cw);

sizes = size(Phi);
T1w = reshape(abs(temp), Ny, Nx, sizes(2), sizes(5));
implay(abs(T1w())/cw);
save(cat(2,save_dir, 'T1w_slc3_Diastole.mat'), 'T1w');

mask_rect = mask(:,:,slc);
save(cat(2,save_dir, 'mask_rect.mat'), 'mask_rect');

%% T1 map values
figure('Position', [100 100 900 900]);
% imagesc(t1_map_3d_nt(:,:,3,1));
mapoi = t1_map_3d_nt(:,:,3,1);
imagesc(mapoi); axis image;
roi_coord = drawpolygon;
roi = createMask(roi_coord);

mean(nonzeros(roi .* mapoi))

% 356 %      356
% 337 - 2    365
% 394 - 3    416
% 432 - 4    458
% 441 - 5    468
% 479 - 6    504
% 488 - 7    518
% 481 - 8    496
% 509 - 9    527
% 513 - 10   523
% 507 - 11   520
%            545
%           
% 535 - 14   548
% 563 - 15   575

% Sofia_D5
% t1_dict_fit = [451.2, 544.6, 526.7, 590.6, 591.8, 602.1, 610.5, 618.7,
% 616.4, 634.8, 662.9, 643.8, 640.0, 674.6, 714.0]; % eliminate first 20
% time points
% t1_dict_fit = [426.9, 518.4, 493.2, 554.8, 551.9, 564.9, 574.8, 592.3,
% 704.1, 594.4, 624.6, 659.6, 608.8, 714.0, 704.4]; % 

%% qMRLab
%% This fitting doesn't work properly, 
% How to make contraints on one of a fitting variables?

x = [2 3 4 5 6 7 8 9 10 11 12 15 40]';
y = [337, 394, 432, 441, 479, 488, 481, 509, 513, 507, 535, 563 645]';


% x = [1 2 3 4 5 6 7 8 9 10 11 12 14 15 40]';
% y = [356, 365, 416, 458, 468, 504, 518, 496, 527, 523, 520, 545, 548, 575, 645]';

%g = fittype('a * exp(b*x) + c');
%f0 = fit(x,y,g,'StartPoint',[100;1;100]);
%f0 = fit(x,y,g);

f = @(b,x) b(1).*exp(x./b(2))+b(3);  
options = optimset('PlotFcns','optimplotfval','TolX',1e-7);

% oldopts = optimset('Display','iter','TolX',1e-6);
B = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100], options)                  % Estimate Parameters

xx = -40:40;
figure();
plot(x,y,'o')
hold on;
plot(xx,f(B,xx),'r-');
grid on;
xlabel('x')
ylabel('f(x)')
text(27, 400, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
%% This seems to be working
figure();
modelFun = @(b,x) b(1).*exp(b(2).*(x))+b(3);  
start = [-1000; -1; 100];

nlm = fitnlm(x,y,modelFun,start);
xx = linspace(0,40)';
plot(x,y,'o'); hold on;
line(xx,predict(nlm,xx),'linestyle','--','color','k')
hold off;
grid on;
xlabel('Time (min)')
ylabel('T1 (ms)')
text(20, 400, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', table2array(nlm.Coefficients(:,1))));
%set(gca,'FontSize',18)

w = [1 1 1 1 1 1 1 1 1 1 1 1 10]';
wnlm = fitnlm(x,y,modelFun,start,'Weight',w);

figure();
plot(x,y,'o'); hold on;
line(xx,predict(wnlm,xx),'color','b');
hold off;
grid on;
xlabel('Time (min)')
ylabel('T1 (ms)')
text(20, 400, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', table2array(nlm.Coefficients(:,1))));

%% Display T1 maps
img_save = cat(2, fid_path, 'T1/');
if ~exist(img_save, 'dir')
   mkdir(img_save); 
end

nsec = 15;
figure();
for i = 1:Nz
    t1_map = squeeze(flip(imrotate(t1_map_3d_nt(:,:,i,:),90),2));
    montage(t1_map); caxis([200 800]);  %colormap hot;
    
    fname = cat(2, fid_file(1:17), fid_file(end-15:end-6), num2str(nsec), '_slc', num2str(i), '_T1Map.png');
    saveas(gcf,cat(2, img_save, fname));
end
