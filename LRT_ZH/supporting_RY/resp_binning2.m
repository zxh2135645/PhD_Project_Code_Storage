h=figure;
imshow(logical(abs(SEs(:,:,1))).*abs(reshape(U(1:(Ny*Nx))*Phi_rt(1,end),Ny,Nx,[])),[]);
title('Select Respiratory ROI')
resp_roi=imrect;
resp_roi=createMask(resp_roi);
close(h)
drawnow

[ry, rx]=size(find(sum(resp_roi,2))*find(sum(resp_roi,1)));
resp=(squeeze(mean(reshape(reshape(U(repmat(resp_roi(:),[L 1])),[],L)*Phi_rt,ry,rx,[]),2)));

for j=1:(Nseg/4)
minme(j)=norm(diff(abs(resp(:,j:(Nseg/4):end)),1,2));
end
[~,argmin]=min(minme);
magresp=abs(resp(:,argmin:(Nseg/4):end));

respmask=magresp>median(magresp(:));
respmask=activecontour(magresp,respmask);

ma = zeros(size(magresp,2),1);
for j=1:numel(ma)
  ma(j)=find([respmask(:,j); 1],1);
  respmask(ma(j):end,j)=1;
end
% ma2 = interp([ma(1); ma],Nseg/4);
% ma2 = ma2((Nseg/4 - j - 1):size(resp,2));
ma2 = interp1(argmin:(Nseg/4):size(resp,2),ma,1:size(resp,2),'pchip',NaN);
temp=1:size(resp,2);
ma2 = interp1(temp(~isnan(ma2)),ma2(~isnan(ma2)),1:size(resp,2),'nearest','extrap');

figure
imagesc(magresp),colormap gray;
hold all
plot(ma,'r')

figure
imagesc(abs(resp)),colormap gray;
hold all
plot(ma2,'r')

Rsig=round((10*ma2)+1);
[Rn,~]=hist(Rsig,1:max(Rsig));
Rn=cumsum(Rn)/sum(Rn);
Rsig=Rn(Rsig);

[~,Ridx]=histc(Rsig,linspace(0,1,rbins+1));
Ridx(Ridx==rbins+1)=rbins;
[~,~,Ridx]=unique(Ridx.');
