% %% Fit Body Unic Table of Contents
% 
% %% extract Bz B0
% close all
% clear all
% %% load files
% BzFolder = 'D:\Dropbox\ShimCoil\Fieldmaps\UnicFieldmap Results_Bz\UNIC_B0Map_05032018.mat';
% load(BzFolder);
% UNIC_BzMap = UNIC_B0Map;
% B0Folder = 'D:\Dropbox\ShimCoil\Fieldmaps\Yang_25881_2103_invovo_GRE2Echo_B0\Result\UNIC_B0Map5881_2103.mat';
% load(B0Folder);
% clear s;
% %% iterate through the struct to create Bz from files
% %BzUnic(:,:,:,1) = permute(UNIC_BzMap.B44.CH5{1},[1,3,2]);
% ChCurrent = 5;
% Bz_test = zeros(224,84,60,42); 
% %Bz has the three dimension the field,60 is the slices, 42 is the coil number
% jr = 1;
% Balun = fieldnames(UNIC_BzMap);
%  for i = 3: length(Balun)
%  ChNum = fieldnames(UNIC_BzMap.(Balun{i}))
%  for j=1: length(ChNum)
%  addon =permute(UNIC_BzMap.(Balun{i}).(ChNum{j}){ChCurrent},[1 3 2]);
%  Bz_test(:,:,:,jr) = Bz_test(:,:,:,jr) + addon;
%  jr = jr+1;
%  end
%  end
% Bz_test=Bz_test/5;
% 
% %% Load B0 and plot each slice of B0
% B0_test = UNIC_B0Map.InVIVO_BH;
% Bz_30slice = zeros(224,84,30,42);
% Bz_30slice(:,:,:,:) = Bz_test(:,:,16:45,:);
% 
% %% convert NaN to 0
% B0_test(isnan(B0_test)) = 0;
% Bz_30slice(isnan(Bz_30slice)) = 0;
% %% save files
% save('U:\Randy\Unic\Invivo studies\Code\Fitting\Bz_test.mat','Bz_test')
% save('U:\Randy\Unic\Invivo studies\Code\Fitting\Bz_30slice.mat','Bz_30slice')
% save('U:\Randy\Unic\Invivo studies\Code\Fitting\B0_test.mat','B0_test')
clear all
load('U:\Randy\Unic\Invivo studies\Data\Bz\Bz_30slice.mat')
addpath('U:\Randy\Unic\Invivo studies\Code\ReadData')
%% get B0
%choose file

B0mapFromRaw
eval( strcat('B0_test =UNIC_B0Map.Map.',string(fieldnames(UNIC_B0Map.Map)),';'));
%clearvars -except DCf stdf slice0
slice0 = [8];
DC_limit0 = [3.5];
maskpercent=0.1;
%% import
B01=B0_test;
Bz=Bz_30slice;
 B0(isnan(B0)) = 0; 
 Bz(isnan(Bz)) = 0;
  B01(isnan(B01)) = 0;
%for ss = 1:length(slice0)
 slice = slice0(1);
 slab = slice-1:slice+1;

 slab = 1:size(B01,3); %whole heart

 Bz5 = Bz;
 clear Bz Bz_f

B0 = B01(:,:,slab);
 mask = (B0 ~= 0);
 mask1 = (B01 ~= 0);
 mask2 = fix(mask); mask2(mask==0) = NaN;

 for i=1:size(mask1,3)
 mask_ss = mask1(:,:,i);
 B0_ss = B01(:,:,i)*42.57e6;
 B0_ss = B0_ss(mask_ss~=0);
 mB0_ss(i) = mean(B0_ss(:));
 stdB0_ss(i) = std(B0_ss(:));
 end
 stdB0_ss_1 = reshape(stdB0_ss,[6,5])';
 %
 so1 = 7;
 so2 = 7;

%% 
 m = 100; % upper & lower limits for the B0 maps (Hz)
 %for k = 1:2
 Bz0 = Bz5;
 [nx,ny,nz,nc] = size(Bz0);
 for i=1:size(Bz0,4)
 Bz(:,:,:,i) = reshape(Bz0(:,:,:,i), [nx ny nz]).*reshape(mask,[nx ny nz]);
 end

 for i = 1:size(Bz0,4)
 Bz1 = Bz(:,:,:,i);
 Bzf(:,i) = Bz1(:);
 end
 B0f = B0(:);
 B0f(isnan(B0f)) = 0; tic
 Bzf(isnan(Bzf)) = 0;
%% Solving
 for dc = 1:length(DC_limit0)
 lb0 = -ones(nc,1)*DC_limit0(dc);
 ub0 = ones(nc,1)*DC_limit0(dc);
 X0 = zeros(nc,1); % initial val
 qq = 5;
 options7 = optimset('Algorithm','trust-region-reflective','MaxFunEvals',1e14,'MaxIter',...
1e12,'TolFun',1e-30,'TolX',1e-16,'Disp','Iter');%,'Largescale','off',);
 end
 [X,resnorm,residual] = lsqlin(Bzf,B0f,[],[],[],[],lb0,ub0,X0,options7);
 % X0 = zeros(1,nc); % initial values
 % opt = optimset('Algorithm','levenbergmarquardt','Display','iter',...
 %'MaxIter',5000,'MaxFunEvals',2000*length(X0),'TolX',1e-5,'TolFun',1e-8);
 % X = fsolve(@(X) fun(X,Bz,B0,mask),X0,opt);
 DC = reshape(X,[nc 1]); % optimal DC current in each coil(A)
 DC
 shimming= reshape(sum(Bzf*DC,2), size(B0));
 B0shim=B0-shimming;
 %% Shimmed B0map
 B0mapFromRaw_Shimmed
 
 %%
 for i=1:size(Bzf,2)
 shimming_each(:,:,:,i) = reshape(Bzf(:,i)*DC(i),size(B0));
 end
 if length(slab) ==1
 kkk = 1;
 else
 kkk = 2;
 end
 m2 = 100;
 figure(124)
 N = size(shimming_each,4);
 for jr = 1:N+1
 if jr == N+1
 subplot(so1,so2,jr)
 imshow(flipdim(shimming(:,:,kkk),1)*42.57e6,[-m2 m2]), colormap jet
 h = title('Total');
 else
 subplot(so1,so2,jr)

 imshow(flipdim(shimming_each(:,:,kkk,jr),1)*42.57e6,[-m2 m2]),
 colormap jet
 h = title(sprintf('%d',jr));
 end
 end
 set(gcf,'color',[1 1 1])
 % Bzshim = Bz .* repmat(reshape(DC,[1 1 1 nc]),[nx ny n1]);
 % B0 maps generated by the multi-coil array with optimal DC currents (T)
 B0shim = (B0 - shimming).* mask; % B0 map after shimming(T)
% B0shimff(ss,k,dc,:,:,:) = B0shim;
% B0ff(ss,k,dc,:,:,:) = B0;
 m = 100;
 m1 = 30;
 kk = 1:length(DC);
 DC_f = [kk', DC];
 fprintf('%d:%12.4f\n',DC_f')
 DC_info = [sum(abs(DC)) mean(abs(DC)) max(abs(DC))];
 fprintf('sum:%12.4f\nmean:%12.4f\nmax:%12.4f\n',DC_info)
DCMat=reshape(DC,[7,6]);