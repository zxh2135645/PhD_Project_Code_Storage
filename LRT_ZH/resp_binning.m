%%
clear resp* bg* ma U1 Rsig
temp_im = reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L);
temp_im = reshape(temp_im*Phi_rt_small_init(:,1),Ny,Nx,[]);
cw = max(temp_im(:));

h = figure;
imshow(abs(temp_im)/abs(cw));
title('Select Respiratory ROI')
% set(gcf,'Position',[0 10 900 400]);
resp_roi = imrect;
resp_roi = createMask(resp_roi);

% bg_roi = impoly;
% bg_roi = createMask(bg_roi);
%%
[ry, rx]=size(find(sum(resp_roi,2))*find(sum(resp_roi,1)));

U1 = vec((dispim(reshape(U,Ny,Nx,Nz,[]),1)));
resp = abs(reshape(reshape(U1(repmat(resp_roi(:),[L_init 1])),[],L_init)*Phi_rt_small_init(:,1:8:end),ry,rx,[]));
implay(abs(resp/max(resp(:))));

% bg = abs(reshape(U1(repmat(bg_roi(:),[L_init 1])),[],L_init)*Phi_rt_small_init(:,1:8:end));
% bg_mean = mean(bg);
% bg_std = std(bg);
% figure,plot(bg_mean(1:1000));

% resp_bg = reshape(reshape(resp,[],size(resp,3))-bg_mean,ry,rx,[]);
% resp_bg(resp_bg<0) = 0;


% implay(abs(resp_bg/max(resp_bg(:))));
%%

se1 = strel('disk',5);

parfor i = 1:size(resp,3)
    resp_rm(:,:,i) = imopen(abs(resp(:,:,i)),se1);
    resp_rm(:,:,i) = resp_rm(:,:,i)./max(max(resp_rm(:,:,i)));
end

resp_th = zeros(size(resp));
tic;
parfor i = 1:size(resp,3)
    [counts,x] = imhist(resp_rm(:,:,i),6);
    T = otsuthresh(counts);
    resp_th(:,:,i) = imbinarize(resp_rm(:,:,i),T);
end
toc;
implay(resp_th)

%%

for i = 1:ry
    parfor j = 1:size(resp_th,3)
        resp_norm(i,j) = find(resp_th(i,:,j),1,'last');
    end
    resp_norm(i,:) = (resp_norm(i,:)-min(resp_norm(i,:)))/range(resp_norm(i,:));
end

Z = sum(resp_norm,1);
figure,plot(Z);

%%
fs=1/(params.lEchoSpacing*SGblock);
df=fs/size(resp,3);
winlp=2*floor((30/60)/df); %highest HR: 40 bpm
windowlp=zeros(size(resp,3),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((10/60)/df); %lowest HR: 3 bpm
windowhp=zeros(size(resp,3),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));

%%

Z=real(ifft(fft(Z).*hwindow));






%%

% figure, imagesc(resp_th(:,:,1),[0 1]), axis('image'), colormap(gray(256))
% % disp('input [center]');
% % keyboard;
% center = ceil(rx*3/10):floor(rx*7/10);
% respmask = squeeze(mean(resp_th(:,center,:),2));
% figure, imagesc(respmask)
% ma = zeros(size(respmask,2),1);
% respmask1 = respmask;
% for j=1:numel(ma)
%     ma(j)=find([respmask(:,j); 1],1);
%     respmask1(ma(j):end,j)=1;
% end
% 
% hold on, plot(ma,'r')
% % get back those beginning lines


%%
fs=1/(params.lEchoSpacing*SGblock);
df=fs/size(resp,3);
winlp=2*floor((30/60)/df); %highest HR: 40 bpm
windowlp=zeros(size(resp,3),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((10/60)/df); %lowest HR: 3 bpm
windowhp=zeros(size(resp,3),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));

% Rsig=real(ifft(fft(ma(:)).*hwindow));

%%
  Z = squeeze(sum(sum(resp,1),2));
  figure,plot(Z);
  Z = real(ifft(fft(Z).*hwindow));
  figure,plot(Z);
  max_Z = mean(findpeaks(Z));
  min_Z = (-1)*mean(findpeaks(Z.*(-1)));
% %   min_Z = min(Z);
  range_Z = max_Z-min_Z;
%   
  Z = ((Z-min_Z)/range_Z)*rbins;
% 
  Z=ceil(Z);
  Z(Z<1)=1;
  Z(Z>rbins)=rbins;
  Z=Z(:);
  figure,plot(Z);
  
  Ridx = Z;

