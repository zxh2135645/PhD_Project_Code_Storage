%% FF quantification for VTR
%%
%% Test fat-water algorithms on (Peter Kellman's) data
%%
%% Author: Diego Hernando
%% Date created: August 18, 2011
%% Date last modified: February 29, 2011

% Add to matlab path
clearvars -except roi
BASEPATH = './';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);
if isunix
    addpath('/home/tianle/supporting/colormap')
    load('/mnt/LiDXXLab/Files/Tianle/cardiac_mri/volunteer_study/V_06222020/meas_MID00069_FID03822_CV_multitaskingT1T2T2star_20200702T095114/multitasking.mat');
else
    addpath('Z:\Tianle\colormap')
    load('Z:\Tianle\colormap\T1cm.mat')
    load('Z:\Tianle\colormap\T2cm.mat')
    [file,path] = (uigetfile('Z:\Tianle\cardiac_mri\*'));
%     load('Z:\Tianle\cardiac_mri\volunteer_study\V_06222020\meas_MID00069_FID03822_CV_multitaskingT1T2T2star_20200702T095114\multitasking.mat')
%     load('/mnt/LiDXXLab/Files/Tianle/cardiac_mri/volunteer_study/V_06222020/meas_MID00069_FID03822_CV_multitaskingT1T2T2star_20200702T095114/multitasking.mat');
%     load('Z:\Tianle\cardiac_mri\phantom study\phantom_09212020\FA5_6echo\meas_MID00069_FID03766_BEAT_MT_cardiac_T1T2T2star_6echo_phantom_20200923T161859\multitasking.mat')
%     load('Z:\Tianle\cardiac_mri\phantom study\phantom_10162020\meas_MID00162_FID02720_BEAT_MT_cardiac_T1T2T2star_9echo_20201018T221315\multitasking.mat')
    load(fullfile(path,file))
end



%% Load some data
if (~exist('cbin_use','var'))
    cbin_use = 11;
end
rbin_use = 1;
L_old = L;
temp = reshape(U,Ny*Nx,[]);
% for tensor recon
start_n = 1;
Phis = permute(Phi(:,:,end,cbin_use,rbin_use,end), [1 2 3 6 4 5]);  % put T2* into T1 recovery, d1:L, d2: T2*, d3: T1, d4: T2, d5: Heart d6:Respiratory
clear Phi
recon = dispim(reshape(temp*(Gr\reshape(Phis,L,[])),Ny,Nx,[]));
% recon = fliplr(recon);
I = recon(:,:,1:end);
I= reshape(I,size(I,1),size(I,2),1,1,size(I,3));
imDataParams.images = double(I);
imDataParams.FieldStrength = 3;
imDataParams.PrecessionIsClockwise = 1;
imDataParams.TE = params.alTE_seconds(1:end);
%% define mask
im_mask=abs(recon(:,:,end-mecho+1));
[im_mask,centroids]=kmeans(im_mask(:),4);
[~,air]=min(centroids);
im_mask=reshape(im_mask~=air,size(recon,1),size(recon,2));
%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 150]; % Range of R2* values, was [0,150]
algoParams.NUM_R2STARS = 65; % Numbre of R2* values for quantization, was 65
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut 
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
tic
  outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

%% Recon -- mixed fit for phase error correction
    % Initialize mixed fitting to graph cut solution
    algoParams.fieldmap = outParams.fieldmap;
    algoParams.r2starmap = outParams.r2starmap;
    algoParams.NUM_MAGN = 1;
    algoParams.THRESHOLD = 0.04;
    algoParams.range_r2star = [0 200];

    % Do mixed fitting
    %% (Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging 
    %% using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)
%     outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams );
%% display the results
water = outParams.species(1).amps;
fat = outParams.species(2).amps;
figure; imagesc(abs(water)); colormap('gray')
axis('image')
axis off 
figure; imagesc(abs(fat)); colormap('gray')
axis('image')
axis off 
R2s = outParams.r2starmap;
delta_B0 = outParams.fieldmap;
figure; imagesc(1./R2s,[0,0.06]); colormap(T2colormap); colorbar
axis('image')
axis off 
figure; imagesc(delta_B0); colormap('parula'); colorbar
axis('image')
axis off 
denom = (abs(fat + water));
denom2 = denom;
denom2(denom==0) = 1; % To avoid divide-by-zero issues
% FF = F/(F+W), F>W
%      1- W/(F+W). F<W
ff_bias_corr = 100*abs(fat)./denom2;
ff_bias_corr(abs(fat)<abs(water)) = 100*(1-abs(water(abs(fat)<abs(water)))./denom2(abs(fat)<abs(water)));
FW_ratio_bias_corr = zeros(size(ff_bias_corr));
FW_ratio_bias_corr(abs(fat)<abs(water)) = 100*(denom2(abs(fat)<abs(water))./abs(water(abs(fat)<abs(water)))-1);
FW_ratio_bias_corr(abs(fat)>=abs(water)) = 100*abs(fat(abs(fat)>=abs(water)))./(denom2(abs(fat)>=abs(water))-abs(fat(abs(fat)>=abs(water))));

denom = (abs(fat) + abs(water));
denom2 = denom;
denom2(denom==0) = 1; % To avoid divide-by-zero issues
% FF = F/(F+W), F>W
%      1- W/(F+W). F<W
ff = 100*abs(fat)./denom2;
FW_ratio = 100*abs(fat)./abs(water);
% ff(abs(fat)<abs(water)) = 100*(1-abs(water(abs(fat)<abs(water)))./denom2(abs(fat)<abs(water)));
% ff = rot90(ff,2);
figure; imagesc(ff,[0,100]); colormap('gray')
axis('image')
axis off 
ff_masked = im_mask.*ff;
T2star = 1./R2s;
save(fullfile(path,'fits_B1fromInvOnly.mat'),'T2star','ff','ff_masked','FW_ratio','FW_ratio_bias_corr','ff_bias_corr','delta_B0','-append')
RGB = ind2rgb(uint16(T2star./.06*2400),T2colormap);
% figure, imshow(RGB)
imwrite(RGB, [path 'multitaskingT2star.png']);
RGB = ind2rgb(uint16(ff./100*1024),gray(1024));
% figure, imshow(RGB)
imwrite(RGB, [path 'multitaskingFF.png']);
RGB = ind2rgb(uint16(ff_bias_corr./100*1024),gray(1024));
% figure, imshow(RGB)
imwrite(RGB, [path 'multitaskingFF_bias_corr.png']);
RGB = ind2rgb(uint16(ff_masked./100*1024),gray(1024));
% figure, imshow(RGB)
imwrite(RGB, [path 'multitaskingFF_masked.png']);
RGB = ind2rgb(uint16((delta_B0+150)./350*1024),parula(1024));
imwrite(RGB, [path 'multitaskingB0.png']);
% %%
% % alpha_deg = params.adFlipAngleDegree;
% % [curveU,curveS]=gen_curve_subspaceT1T2_VTR(Nseg,params.shortTR,params.lEchoSpacing,alpha_deg,1);
% % cL = 8;
% % curvePhi = curveU(:,1:cL);
% temp = reshape(U,Ny*Nx,[]);
% load(fullfile(path,file),'Phi','L')
% Phi_temp = squeeze(permute(Phi(:,:,:,cbin_use,rbin_use,:),[1,2,4,5,3,6]));
% clear Phi
% Phi_temp = reshape(Phi_temp,L,mecho,[]);
% % recon = dispim(reshape(temp*(Gr\Phi_temp),Ny,Nx,mecho,[]));
% %%
% Deltaf = algoParams.species(2).frequency;
% Deltaf = 42.58*3*Deltaf;
% TE = imDataParams.TE ;
% perc = algoParams.species(2).relAmps;
% Expf = exp(1i*2*pi*Deltaf.'*TE);
% Fat_spect = perc*Expf;
% imDataParams_new = imDataParams;
% 
% for L_idx = 1:Nseg*5
%     recon = dispim(reshape(temp*(Gr\Phi_temp(:,:,L_idx)),Ny,Nx,mecho,[]));
%     recon= reshape(recon,size(recon,1),size(recon,2),1,1,size(recon,3));
%     imDataParams_new.images = double(recon(:,:,:,:,1:end));
%     algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1; 
%     % r2starmap = estimateR2starGivenFieldmap( imDataParams_new, algoParams, delta_B0 );
%     amps = decomposeGivenFieldMapAndDampings( imDataParams_new,algoParams, delta_B0,R2s,R2s);
%     % amps = decomposeGivenFieldMapAndDampings( imDataParams_new,algoParams, delta_B0,r2starmap,r2starmap);
%     recon_water(:,:,L_idx) = squeeze(amps(:,:,1,:));
%     recon_fat(:,:,L_idx) = squeeze(amps(:,:,2,:));
% end
% %% fitting based on that
% % fitw=sqrt(sum(reshape(lrw(1:mecho:end).^2,Nseg*5/2,[]),2));
% lrw_rep = vec(repmat(lrw(1:mecho:end),[2,1])).';
% fitw=sqrt(sum(reshape(lrw_rep.^2,Nseg*5,[]),2));
% fitw=fitw.'/median(fitw);
% % fitw(1:10) = 0;
% isVE=1;
% %%
% TR1 = params.shortTR;
% TR2 = params.lEchoSpacing;
% e11 = @(T1)exp(-TR1./T1);
% e12 = @(T1)exp(-TR2./T1);
% M_ss_long = @(e11,e12,alpha) (1-e12+cos(alpha).*e12-cos(alpha).*e12.*e11)./(1-cos(alpha).^2.*e12.*e11);
% M_ss_short = @(e11,e12,alpha) (1-e11+cos(alpha).*e11-cos(alpha).*e12.*e11)./(1-cos(alpha).^2.*e12.*e11);
% 
% vec = @(x) x(:);
% TEs = [0 30 40 50 60]*1e-3;
% % Omega = pi/(1.22*1e-3);
% %Es = [1 2 3 4 ]*2.46*1e-3;
% e2 = @(T2)exp(-TEs/T2);
% %n = Nseg-Nseg+1:2:Nseg;
% n = 1:1:Nseg/2;
% clear Eff*
% %Calculate pct. of steady state reached
% syms T1 T2 B BT2 Eff_last alph %can't use "alpha"?
% Eff_prev=Eff_last;
% 
% j=1; %special case (IR, VE only)
% if isVE
%   Eff_sym(j)=(1 + ((exp(-(TR1+TR2)/T1)*cos(alph)^2)^(Nseg/2)).'*(cos(B*pi)*Eff_prev*exp(-TEs(j)/T2)-1)); %B is relative B1 field, not inversion efficiency
% else
%   Eff_sym(j)=(1 + ((exp(-(TR1+TR2)/T1)*cos(alph)^2)^(Nseg/2)).'*((cos(BT2*pi/2)^2-sin(BT2*pi/2)^2*exp(-TEs(j)/T2)*Eff_prev-1)));
% end
% 
% Eff_prev=Eff_sym(j);
% 
% for j=2:numel(TEs)
%   if isVE
%     Eff_sym(j)=(1 + ((exp(-(TR1+TR2)/T1)*cos(alph)^2)^(Nseg/2)).'*(cos(B*pi)*Eff_prev*exp(-TEs(j)/T2)-1)); %B is relative B1 field, not inversion efficiency
%   else
%     Eff_sym(j)=(1 + ((exp(-(TR1+TR2)/T1)*cos(alph)^2)^(Nseg/2)).'*((cos(BT2*pi/2)^2-sin(BT2*pi/2)^2*exp(-TEs(j)/T2))*Eff_prev-1)); %B is relative B1 field, not inversion efficiency
%   end
%   Eff_prev=Eff_sym(j);
% end
% Eff_last=solve(Eff_last==Eff_sym(end),Eff_last);
% Eff_sym(end)=Eff_last;
% evalstring='Eff=@(T1,T2,alph,B,BT2)cat(2';
% for j=1:numel(TEs)
%   evalstring=strcat(evalstring,sprintf(',%s',char(subs(Eff_sym(j)))));
% end
% evalstring=strcat(evalstring,');');
% eval(evalstring);
% Sint1 = @(A,e11,e12,e2,alpha,B,BT2,Eff)vec(A * M_ss_long(e11,e12,alpha) * (1 + ((cos(alpha)^2*e11*e12).^(n-1)).'*(([cos(B*pi), -sin(BT2*pi/2)^2*ones(1,4)].*e2+cos(BT2*pi/2)^2*[0 ones(1,4)]).*Eff-1))).';
% Sint2 = @(A,e11,e12,e2,alpha,B,BT2,Eff)vec((A * M_ss_short(e11,e12,alpha) + A * e11 * cos(alpha) * M_ss_long(e11,e12,alpha) * ((cos(alpha)^2*e11*e12).^(n-1)).'*(([cos(B*pi), -sin(BT2*pi/2)^2*ones(1,4)].*e2+cos(BT2*pi/2)^2*[0 ones(1,4)]).*Eff-1))).';
% Sint = @(A,e11,e12,e2,alpha,B,BT2,Eff)vec([Sint1(A,e11,e12,e2,alpha,B,BT2,Eff);Sint2(A,e11,e12,e2,alpha,B,BT2,Eff)]).';
% % S = @(A,T1,T2,alpha,B)Sint(A,e1(T1),e2(T2),alpha,B,B,circshift(Eff(T1,T2,alpha,B,B),[0, 1]));
% %right now, assume B = BT2
%   S = @(A,T1,T2,alpha,B)Sint(A*.9547,e11(T1),e12(T1),e2(T2),alpha*.9547,B,1,circshift(Eff(T1,T2,alpha*.9547,B,1),[0, 1])) ...
%     +Sint(2*A*.5513,e11(T1),e12(T1),e2(T2),alpha*.5513,B,1,circshift(Eff(T1,T2,alpha*.5513,B,1),[0, 1])) ...
%     +Sint(2*A*0.0735,e11(T1),e12(T1),e2(T2),alpha*0.0735,B,1,circshift(Eff(T1,T2,alpha*0.0735,B,1),[0, 1]));
% ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)
% %%
% im_mask=abs(recon_water(:,:,end));
% [im_mask,centroids]=kmeans(im_mask(:),4);
% [~,air]=min(centroids);
% im_mask=reshape(im_mask~=air,size(recon,1),size(recon,2));
% im_mask = true(size(recon,1),size(recon,2));
% % im_mask = createMask(impoly);
% figure,imshow(im_mask),drawnow;
% 
% [rows, cols]=ind2sub(size(im_mask),find(im_mask));
% recontemp=reshape(recon_water,size(recon,1)*size(recon,2),[]);
% % recontemp = reshape(squeeze(SimX(:,1,:,:)),size(squeeze(SimX(:,1,:,:)),1)*size(squeeze(SimX(:,1,:,:)),2),[]);
% recontemp=recontemp(im_mask(:),:);
% % recontemp = bsxfun(@times,recontemp,weights_temp.');
% % opts=[];
% % opts.MaxFunEvals = 4000;
% opts  = optimoptions(@lsqnonlin,'MaxFunEvals',1000,'Maxiter',1000);
% %%
% alpha = 3.5*pi/180; %from T1MES phantom
% x0 = [1.1, 40e-3, 1];
% xlb = [100e-3, 5e-3, .5];
% xub = [3, 150e-3, 1];
% tic
% fitmat=zeros(numel(rows),5);
% parfor j=1:numel(rows)
%     curve = double(recontemp(j,:));
%     normcurve = abs(curve(end)+1e-10);
%     curve = curve/normcurve;
%     if norm(curve)>0
%         Avp = @(T1,T2,B)ppinv(S(1,T1,T2,alpha,B).*fitw,curve.*fitw); %parameterize solution to A as function of T1,alpha,B
%         cost=@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha,x(3))-curve).*fitw;
%         [tempfit, res]=lsqnonlin(cost,x0, xlb, xub, opts);
%         tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit];
%         tempfit(5)=sqrt(res)*abs(normcurve);
%         fitmat(j,:) = tempfit;
%     end
% end
% toc
% tempfit=fitmat;
% fits=zeros(size(recon,1)*size(recon,2),5);
% fits(im_mask(:),:)=tempfit;
% fits=reshape(fits,size(recon,1),size(recon,2),5);
% %%
% %%
% T1_water = fits(:,:,2);
% T2_water = fits(:,:,3);
% T1_water(T1_water==Inf) = 0;
% T2_water(T2_water==Inf) = 0;
% 
% figure,imagesc(T1_water,[0 3]),axis('image'), axis off, colorbar, colormap(T1colormap)
% figure,imagesc(T2_water,[0 0.15]),axis('image'), axis off, colorbar, colormap(T2colormap)
% [histn,histx]=hist(abs(fits(:,:,1)),1000);
% cw=histx(find(cumsum(histn)/sum(histn)>.999,1));
% figure; imshow(abs(fits(:,:,1))/cw)
% water_M0 = (fits(:,:,1));
% save(fullfile(path,'water_fits.mat'),'T1_water','T2_water','water_M0')
% %%
% RGB = ind2rgb(uint16((T1_water)./3*2700),T1colormap);
% figure, imshow(RGB)
% imwrite(RGB, [path 'T1_water.png']);
% RGB = ind2rgb(uint16(T2_water./.15*2400),T2colormap);
% figure, imshow(RGB)
% imwrite(RGB, [path 'T2_water.png']);
% RGB = ind2rgb(uint16(abs(water_M0)./cw*1024),gray(1024));
% figure, imshow(RGB)
% imwrite(RGB, [path 'water_M0.png']);
% %%
% im_mask=abs(recon_fat(:,:,end));
% [im_mask,centroids]=kmeans(im_mask(:),4);
% [~,air]=min(centroids);
% im_mask=reshape(im_mask~=air,size(recon,1),size(recon,2));
% im_mask = true(size(recon,1),size(recon,2));
% % im_mask = createMask(impoly);
% figure,imshow(im_mask),drawnow;
% 
% [rows, cols]=ind2sub(size(im_mask),find(im_mask));
% recontemp=reshape(recon_fat,size(recon,1)*size(recon,2),[]);
% % recontemp = reshape(squeeze(SimX(:,1,:,:)),size(squeeze(SimX(:,1,:,:)),1)*size(squeeze(SimX(:,1,:,:)),2),[]);
% recontemp=recontemp(im_mask(:),:);
% % recontemp = bsxfun(@times,recontemp,weights_temp.');
% % opts=[];
% % opts.MaxFunEvals = 4000;
% opts  = optimoptions(@lsqnonlin,'MaxFunEvals',1000,'Maxiter',1000);
% %%
% alpha = 3.5*pi/180; %from T1MES phantom
% x0 = [1.1, 40e-3, 1];
% xlb = [100e-3, 10e-3, .5];
% xub = [3, 150e-3, 1];
% tic
% fitmat=zeros(numel(rows),5);
% parfor j=1:numel(rows)
%     curve = double(recontemp(j,:));
%     normcurve = abs(curve(end)+1e-10);
%     curve = curve/normcurve;
%     if norm(curve)>0
%         Avp = @(T1,T2,B)ppinv(S(1,T1,T2,alpha,B).*fitw,curve.*fitw); %parameterize solution to A as function of T1,alpha,B
%         cost=@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha,x(3))-curve).*fitw;
%         [tempfit, res]=lsqnonlin(cost,x0, xlb, xub, opts);
%         tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit];
%         tempfit(5)=sqrt(res)*abs(normcurve);
%         fitmat(j,:) = tempfit;
%     end
% end
% toc
% tempfit=fitmat;
% fits=zeros(size(recon,1)*size(recon,2),5);
% fits(im_mask(:),:)=tempfit;
% fits=reshape(fits,size(recon,1),size(recon,2),5);
% %%
% %%
% T1_fat = fits(:,:,2);
% T2_fat = fits(:,:,3);
% T1_fat(T1_fat==Inf) = 0;
% T2_fat(T2_fat==Inf) = 0;
% 
% figure,imagesc(T1_fat,[0 1]),axis('image'), axis off, colorbar, colormap(T1colormap)
% figure,imagesc(T2_fat,[0 0.15]),axis('image'), axis off, colorbar, colormap(T2colormap)
% [histn,histx]=hist(abs(fits(:,:,1)),1000);
% cw=histx(find(cumsum(histn)/sum(histn)>.999,1));
% figure; imshow(abs(fits(:,:,1))/cw)
% fat_M0 = (fits(:,:,1));
% %%
% %%
% im_mask=abs(recon_fat(:,:,end));
% [im_mask,centroids]=kmeans(im_mask(:),4);
% [~,air]=min(centroids);
% im_mask=reshape(im_mask~=air,size(recon,1),size(recon,2));
% RGB = ind2rgb(uint16(T1_fat.*im_mask*2700),T1colormap);
% figure, imshow(RGB)
% imwrite(RGB, [path 'T1_fat.png']);
% RGB = ind2rgb(uint16(T2_fat.*im_mask./.15*2400),T2colormap);
% figure, imshow(RGB)
% imwrite(RGB, [path 'T2_fat.png']);
% RGB = ind2rgb(uint16(abs(fat_M0)./cw*1024),gray(1024));
% figure, imshow(RGB)
% imwrite(RGB, [path 'fat_M0.png']);
% %%
% denom = (abs(fat_M0) + abs(water_M0));
% denom2 = denom;
% denom2(denom==0) = 1; % To avoid divide-by-zero issues
% % FF = F/(F+W), F>W
% %      1- W/(F+W). F<W
% PDFF = 100*abs(fat_M0)./denom2;
% PDFF(abs(fat_M0)<abs(water_M0)) = 100*(1-abs(water_M0(abs(fat_M0)<abs(water_M0)))./denom2(abs(fat_M0)<abs(water_M0)));
% figure; imagesc(PDFF,[0,100]); colormap('gray')
% axis('image')
% axis off 
% RGB = ind2rgb(uint16(PDFF./100*1024),gray(1024));
% % figure, imshow(RGB)
% imwrite(RGB, [path 'multitaskingPDFF.png']);
% 
% denom = (abs(fat_M0 + water_M0));
% denom2 = denom;
% denom2(denom==0) = 1; % To avoid divide-by-zero issues
% % FF = F/(F+W), F>W
% %      1- W/(F+W). F<W
% PDFF_bias_corr = 100*abs(fat_M0)./denom2;
% PDFF_bias_corr(abs(fat_M0)<abs(water_M0)) = 100*(1-abs(water_M0(abs(fat_M0)<abs(water_M0)))./denom2(abs(fat_M0)<abs(water_M0)));
% RGB = ind2rgb(uint16(PDFF_bias_corr./100*1024),gray(1024));
% % figure, imshow(RGB)
% imwrite(RGB, [path 'multitaskingPDFF_bias_corr.png']);
% %%
% water_corr = water_M0.*M_ss_short(e11(T1_water),e12(T1_water),alpha);
% denom = abs((fat) + (water_corr));
% denom2 = denom;
% denom2(denom==0) = 1; % To avoid divide-by-zero issues
% % FF = F/(F+W), F>W
% %      1- W/(F+W). F<W
% FF_corr = 100*abs(fat)./denom2;
% FF_corr(abs(fat)<abs(water_corr)) = 100*(1-abs(water_corr(abs(fat)<abs(water_corr)))./denom2(abs(fat)<abs(water_corr)));
% figure; imagesc(FF_corr,[0,100]); colormap('gray')
% axis('image')
% axis off 
% save(fullfile(path,'fits_B1fromInvOnly.mat'),'FF_corr','PDFF','PDFF_bias_corr','-append')
% %% perform single compartment T2* mapping
% TE = params.alTE_seconds;
% S = @(A,R2) abs(A*exp(-R2*TE));
% ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)
% recontemp=reshape(I,size(I,1)*size(I,2),[]);
% %%
% fits = zeros(numel(rows),2);
%     parfor j=1:numel(rows)
%         curve = abs(double(recontemp(j,:)));
% 
%         normcurve = abs(curve(end)+1000*eps);
% %         if normcurve == 0
% %             normcurve = curve(1);
% %         end
%         curve = curve/normcurve;
% 
%         Avp = @(R2)ppinv(S(1,R2),curve); %parameterize solution to A as function of R1,alpha,B
%         tempfit=lsqnonlin(@(x)abs(S(Avp(x(1)), x(1))-curve),...
%             [1/0.024],...
%             [1], [1/0.001], opts); %be careful of alpha bounds. 0.58?
%         tempfit=[Avp(tempfit(1))*normcurve, tempfit];
% 
%     %          Avp = @(R1,R2,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
%         %     tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha, x(3))-curve),...
%         %         [2/3,10,median([minB, real(curve(1)/curve(end)), maxB])],...
%         %         [1/3,1,minB], [1/.1,1/.01,maxB], opts); %be careful of alpha bounds. Lower bound for 3D is alpha/2, maybe?
%         %     tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit(1), tempfit(2), alpha, tempfit(3)];
%         %
%         fits(j,:) = tempfit;
% 
%     end
%     tempfit=fits;
%     fits=zeros(size(recon,1)*size(recon,2),2);
%     fits(:,:)=tempfit;
%     fits=reshape(fits,size(recon,1),size(recon,2),2);
%     
%     T2star_single = 1./fits(:,:,2);
%     RGB = ind2rgb(uint16(T2star_single./.06*2400),T2colormap);
% figure, imshow(RGB)
% imwrite(RGB, [path 'multitaskingT2star_single.png']);
% T2star_single(T2star_single==Inf) = 0;
% save(fullfile(path,'fits_B1fromInvOnly.mat'),'T2star_single','-append')