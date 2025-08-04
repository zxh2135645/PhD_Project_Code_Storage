
temp = prep(U_init,st);
temp = temp/max(abs(temp(:)));
if MBfactor == 1
    temp = fftshift(temp,3);
end
if isCartesian
    temp = fftshift(temp,1);
end

if ~exist('vecOrder','var')
    try
        vecOrder = [zDir xDir yDir];
    catch
        vecOrder = [-3 2 1];
    end
end

temp  = permute(temp,[abs(vecOrder) 4]);
if sign(vecOrder(1)) > 0
    temp = flip(temp,1);
end
dim2 = size(temp,2);
temp1 = sum(sum(abs(temp(:,floor(dim2/4):floor(dim2*3/4),:,1)),1),2);
temp1(temp1<mean(temp1)) = 0;
ky = find(temp1>0,1) + ceil((find(temp1>0,1,'last') - find(temp1>0,1))/2);

Segidx = mod(navIndices-1,linesPerShot) + 1;
if alTR_seconds < 0.5
    bin_frame = linesPerShot:linesPerShot:Nread;
else
    bin_frame = navIndices;
    bin_frame(Segidx < SGBlock*2+1) = [];
end

Phi_rt_bin = Phi_rt_full_init(:,bin_frame); 

%ky = round(Ny/2)+30;
[~,frame] = max(abs(Phi_rt_bin(1,:)));
cw = max(abs(temp(:)));
gray_im = reshape(reshape(temp(:,:,ky,:),[],L_init)*(Phi_rt_bin(:,frame)),pNz,pNx,[]);                              
gray_im = imresize(gray_im, [newNz pNx]);
gray_im = abs(gray_im)/cw;
h = figure;
imagesc(gray_im), colormap(gray), axis('image');
title('Select a point within liver, under liver dome');
set(gcf, 'Position', [0 10 900 400]);
hs = impoint();
hs = round(getPosition(hs));
mask = false(size(gray_im));
mask(max(1,hs(2)-10):min(size(mask,1),hs(2)+10), max(1,hs(1)-10):min(size(mask,2),hs(1)+10)) = true;
% resp_roi=imrect;
% resp_roi=createMask(resp_roi);
% mask = roipoly;
close(h); drawnow;

%% Liver segment and identify liver dome
% % rough segments
% se1 = strel('disk', 5,8);
% gray_im = imopen(gray_im, se1);
% seed = [70,110];
% preload mask
% mask = false(size(gray_im)); 
% mask(seed(1),seed(2)) = true;
% Compute the weight array based on grayscale intensity differences.
W = graydiffweight(gray_im, mask, 'GrayDifferenceCutoff', 25);
% Segment the image using the weights.
thresh = 0.01;
[BW, D] = imsegfmm(W, mask, thresh);
figure, imshow(BW), title('contour of dome')
% % identify liver dome
start = find(sum(BW), 1, 'first');
en = find(sum(BW), 1, 'last');
curve = ones(1,en-start+1);
for i = 0:en-start
    curve(i+1) = find(BW(:,start+i),1);
end
curve = smooth(curve,10);
[y_dome,x_dome] = min(curve);
x_dome = floor(x_dome)+start-1;
y_dome = round(y_dome);
resp_roi = false(size(gray_im));
width_y = 100;
width_x = 20;
resp_roi(max(1,y_dome-50):min(size(resp_roi,1),y_dome+width_y), max(1,x_dome-width_x):min(size(resp_roi,2),x_dome+width_x)) = true;

%%
[ry, rx] = size(find(sum(resp_roi,2))*find(sum(resp_roi,1)));
U1 = vec(imresize(temp(:,:,ky,:),[newNz pNx]));  

remove = find(Segidx<=SGBlock+1);  % remove the first bin
Phi_rt_small_rm = Phi_rt_small;
Phi_rt_small_rm(:,remove) = [];
resp = reshape(reshape(U1(repmat(resp_roi(:),[L_init 1])),[],L_init)*Phi_rt_bin,ry,rx,[]);
s_norm = mean(reshape(U1(repmat(mask(:),[L_init 1])),[],L_init)*Phi_rt_bin,1);

resp = abs(resp)./reshape(abs(s_norm),1,1,[]);
% h = implay(resp);
% set(h.Parent, 'Name','dome ROI');

%% processs image feature
se1 = strel('disk',10,8);
% se1 = strel('disk', 15,8);
% se2 = strel('disk', 3,8);

resp_rm = zeros(size(resp));
resp_th = zeros(size(resp));
tic;
for i = 1:size(resp,3)
    test = imopen(resp(:,:,i),se1);
    test(test>1)   = 1;
    test(test<0.2) = 0;
    resp_rm(:,:,i) = test;
    
    % threshold
    [counts,x] = imhist(test,4);
    T = otsuthresh(counts);
    resp_th(:,:,i) = imbinarize(test,T);
    resp_th(:,:,i) = bwareafilt(logical(resp_th(:,:,i)),1);
end
toc;

h = implay([resp resp_rm resp_th]);
set(h.Parent, 'Name','dome;  opened dome;  threshold');

%% get contour 
% figure, imagesc(resp_th(:,:,1),[0 1]), axis('image'), colormap(gray(256))
% disp('input [center]');
% keyboard;
center = ceil(rx*3/10):floor(rx*7/10);
respmask = squeeze(mean(resp_th(:,center,:),2));
figure, imagesc(respmask)
ma = zeros(size(respmask,2),1);
respmask1 = respmask;
for j = 1:numel(ma)
    ma(j) = find([respmask(:,j); 1],1);
    respmask1(ma(j):end,j) = 1;
end
hold on, plot(ma,'r')
title('Track of liver dome')
% get back those beginning lines

% %% low-pass filter to resp frequency
% fs=1/(params.lEchoSpacing*SGBlock);
% df=fs/size(Phi_rt_small,2);
% hwin=2*floor((30/60)/df);
% hwindow=zeros(size(Phi_rt_small,2),1);
% hwindow(1:hwin)=hamming(hwin,'periodic');
% hwindow=circshift(hwindow,[-hwin/2, 0]);
% Rsig=real(ifft(fft(ma(:)).*hwindow));

%%
%rbins = 6;
Rsig = round((20*ma)+1);
[Rn,~] = hist(Rsig,ceil(min(Rsig)):max(Rsig));
Rn = cumsum(Rn)/sum(Rn);
Rsig = Rn(Rsig-ceil(min(Rsig))+1);
Rsig = medfilt1(Rsig,5);
Rsig(Rsig < prctile(Rsig,5))  = prctile(Rsig,5);
Rsig(Rsig > prctile(Rsig,95)) = prctile(Rsig,95);

[~,Ridx] = histc(Rsig,linspace(min(Rsig),max(Rsig),rbins+1));
Ridx(Ridx==rbins+1) = rbins;
Ridx(Ridx==0) = 1;
[~,~,Ridx] = unique(Ridx.');
% figure, plot(Ridx)
% title('Ridx')

% Ridx(Ridx==1) = 2;
% Ridx = Ridx-1;

%%
X = bin_frame;
Xq = 1:Ntpoint;
Rorg = round(interp1(X,Ridx,Xq,'nearest','extrap'));
clear X Xq
Rorg(Rorg<1) = 1;
Rorg(Rorg>rbins) = rbins;
Ridx = Rorg(navIndices);

rbins = max(Ridx);

% get back those beginning lines

% %% low-pass filter to resp frequency
% fs=1/(params.lEchoSpacing*SGBlock);
% df=fs/size(Phi_rt_small,2);
% hwin=2*floor((30/60)/df);
% hwindow=zeros(size(Phi_rt_small,2),1);
% hwindow(1:hwin)=hamming(hwin,'periodic');
% hwindow=circshift(hwindow,[-hwin/2, 0]);
% Rsig=real(ifft(fft(ma(:)).*hwindow));


%% display resp bins

clear recon
recon = zeros(pNz,pNx,rbins);
for i = 1:rbins
    Phi_idx = mean(Phi_rt_small_init(:,find(Ridx==i)),2);  
    recon(:,:,i) = reshape(reshape(temp(:,:,ky,:),[],L_init)*Phi_idx, pNz, pNx);
end
cw  = max(abs(vec(recon)));
recon = imresize(recon,[newNz,pNx],'bicubic');

h = implayZoom(abs(recon)/cw,3);
set(h.Parent, 'Name','All respiratory bins')

%%
% clear recon
% l = find((Rorg(Nseg:Nseg:end)==1));
% recon = reshape(reshape(temp(:,:,ky,:),[],L_init)*Phi_rt_full_init(:,l*Nseg), Nz, Nx,[]);
% recon = imresize(recon,[newNz,Nx],'bicubic');
% h = implayZoom(abs(recon)/cw);
% set(h.Parent, 'Name','bin = 1 for last SR times')

%%
clear recon
l = find((Ridx==1));
recon = reshape(reshape(temp(:,:,ky,:),[],L_init)*Phi_rt_small_init(:,l), pNz, pNx,[]);
recon = imresize(recon,[newNz,pNx],'bicubic');
cw  = max(abs(vec(recon)));
h = implayZoom(abs(recon)/cw);
set(h.Parent, 'Name','All frames in bin #1')
