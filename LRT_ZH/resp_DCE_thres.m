%% use last TI time for binning
clear resp*
temp = prep(U_init,st);
Phi_rt_small_init = Phi_rt_init(:,nav_indices_small);
Phi_rt_full_init = Phi_rt_init(:,Nseg*mecho:Nseg*mecho:end); % last TI
Segidx=mod(nav_indices_small-1,params.lSegments)+1;
ky = 155;
gray_im = fftshift(reshape(reshape(temp(ky,:,:),[],L_init)*(Phi_rt_full_init(:,55*mecho)),Nz,Nx,[]),1);                              
%gray_im = dispim(reshape(reshape(temp(:,1,:),[],L)*(Phi_rt_init(:,end/4)),Ny,1,Nx,[]));
gray_im = imresize(gray_im, [Nz*5 Nx]);
h=figure;
imagesc(abs(gray_im)/cw*2), colormap(gray)
title('Select Respiratory ROI')
set(gcf, 'Position', [0 10 900 400])
resp_roi=imrect;
resp_roi=createMask(resp_roi);
close(h)
drawnow

%%
[ry, rx]=size(find(sum(resp_roi,2))*find(sum(resp_roi,1)));
%U1 = vec(dispim4(reshape(temp(1,:),1,Nz,Nx,[])));     
U1 = vec(imresize(fftshift(reshape(temp(ky,:,:),Nz,Nx,[]),1),[Nz*5 Nx]));  
%U1 = vec(dispim(reshape(temp(:,1,:),Ny,1,Nx,[])));
remove = find(Segidx<=SGblock+1);  % remove the first bin
resp=reshape(reshape(U1(repmat(resp_roi(:),[L_init 1])),[],L_init)*Phi_rt_full_init,ry,rx,[]);
%resp=squeeze(mean(reshape(reshape(U1(repmat(resp_roi(:),[L 1])),[],L)*Phi_rt_init(:,1:10:end),ry,rx,[]),2));
%resp = resp(end:-1:1,:);
resp = abs(resp)/cw*5;
implay(resp)
figure, imagesc(resp(:,:,floor(end/5)));
rm = input('how many lines to remove from the top?');
resp(1:rm,:,:) = 0;
%% processs image feature
se1 = strel('disk', 15,8);
se2 = strel('disk', 3,8);
resp_rm = zeros(size(resp));
tic;
parfor i = 1:size(resp,3)
    resp_rm(:,:,i) = imopen(imopen(resp(:,:,i),se1),se2);
end
toc;
resp_rm(resp_rm>1) = 1;
resp_rm(resp_rm<0.2) = 0;
implay(resp_rm)

%% threshold
resp_th = zeros(size(resp));
tic;
parfor i = 1:size(resp,3)
    [counts,x] = imhist(resp_rm(:,:,i),4);
    T = otsuthresh(counts);
    resp_th(:,:,i) = imbinarize(resp_rm(:,:,i),T);
end
toc;
implay(resp_th)
%% get contour 
figure, imagesc(resp_th(:,:,1),[0 1]), axis('image'), colormap(gray(256))
% disp('input [center]');
% keyboard;
center = ceil(rx*3/10):floor(rx*7/10);
respmask = squeeze(mean(resp_th(:,center,:),2));
figure, imagesc(respmask)
ma = zeros(size(respmask,2),1);
respmask1 = respmask;
for j=1:numel(ma)
    ma(j)=find([respmask(:,j); 1],1);
    respmask1(ma(j):end,j)=1;
end
hold on, plot(ma,'r')
% get back those beginning lines

%% low-pass filter to resp frequency
fs=1/(params.alTR_seconds);
df=fs/size(resp,3);
hwin=2*floor((30/60)/df);
hwindow=zeros(size(resp,3),1);
hwindow(1:hwin)=hamming(hwin,'periodic');
hwindow=circshift(hwindow,[-hwin/2, 0]);
Rsig=real(ifft(fft(ma(:)).*hwindow));

%%
rbins = 6;
Rsig=round((20*Rsig)+1);
[Rn,~]=hist(Rsig,1:max(Rsig));
Rn=cumsum(Rn)/sum(Rn);
Rsig=Rn(Rsig);

[~,Ridx]=histc(Rsig,linspace(0,1,rbins+1));
Ridx(Ridx==rbins+1)=rbins;
Ridx(Ridx==0)=1;
[~,~,Ridx]=unique(Ridx.');


X = [Nseg:Nseg:Nread/mecho];
Xq = 1:Nread/mecho;
Rorg = round(interp1(X,Ridx,Xq,'nearest','extrap'));
clear X Xq
Rorg(Rorg<1)=1;
Rorg(Rorg>rbins) = rbins;
Ridx = vec(repmat(Ridx.', [mecho 1])); % Ridx, bin for last TI
Rorg = vec(repmat(Rorg.', [mecho 1])); % Rorg, bin for all imaging line
Ridx = Rorg(nav_indices);  % now the Ridx is the bin for SG line
%% 
temp = prep(U_init,st);
clear recon*
for i = 1:rbins
    Phi_idx = mean(Phi_rt_init(:,find(Rorg==i)),2);  
    recon(:,:,i) = fftshift(reshape(reshape(temp(ky,:,:),[],L_init)*Phi_idx, Nz, Nx,[]),1);
end
for i = 1:size(recon,3)
recon1(:,:,i) = imresize(recon(:,:,i),[Nz*5,Nx],'bicubic');
end
implay(abs(recon1)/cw);

%%
clear recon*
l = find((Rorg(Nseg:Nseg:end)==2));

recon = fftshift(reshape(reshape(temp(ky,:,:),[],L_init)*Phi_rt_full_init(:,l), Nz, Nx,[]),1);

for i = 1:size(recon,3)
recon1(:,:,i) = imresize(recon(:,:,i),[Nz*5,Nx],'bicubic');
end
implay(abs(recon1)/cw);

%%
clear  resp* recon*
clear magresp ma gray_im h U1 rx ry resporig Rn
clear recon_ave_new4 recon_ave_new41 recon1