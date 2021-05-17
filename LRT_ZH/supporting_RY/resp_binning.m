h=figure;
imshow(logical(abs(SEs(:,:,1))).*abs(reshape(U(1:(Ny*Nx))*Phi_rt(1,end),Ny,Nx,[])),[]);
title('Select Respiratory ROI')
resp_roi=imrect;
resp_roi=createMask(resp_roi);
close(h)
drawnow

[ry, rx]=size(find(sum(resp_roi,2))*find(sum(resp_roi,1)));
resp=(squeeze(mean(reshape(reshape(U(repmat(resp_roi(:),[L 1])),[],L)*Phi_rt,ry,rx,[]),2)));

resporig=resp;
if ~isempty(curvePhi)
  tempnom=repmat(curvePhi-repmat(mean(curvePhi),[Nseg/2 1]),[size(resp,2)/Nseg*2 1]).';
  Nits = 2;
else
  Nits=1;
end

for its=1:Nits
  magresp=abs(resp);
  
  % xc=abs(ifft((fft(magresp)).*repmat(fft(magresp(:,end)),[1 size(magresp,2)])));
  % xc=xc./repmat(norm(magresp(:,end))*sqrt(sum(abs(magresp).^2)),[size(magresp,1) 1]);
  % [~,ma]=max(xc);
  % [~,ma]=max(medfilt1(xc,101,[],2));
  
  respmask=magresp>median(magresp(:));
  try
  respmask=activecontour(magresp,respmask,'Chan-Vese','SmoothFactor',5);
  catch
       respmask=activecontour(magresp,respmask);
  end
  ma = zeros(size(magresp,2),1);
  for j=1:numel(ma)
    ma(j)=find([respmask(:,j); 1],1);
    respmask(ma(j):end,j)=1;
  end
  
  figure(its)
  imagesc(magresp),colormap gray;
  hold all
  plot(ma,'r')
  
  % ma=medfilt1(ma,101-20*its);
  ma=fft(ma);
  if its==1
    figure,plot(log(abs(ma)))
    xlim([0 70])
    resp_cutoff=input('Cutoff? ');
  end
  ma(round(its*resp_cutoff):(end-round(its*resp_cutoff)+2))=0;
  ma=real(ifft(ma));
  
  figure(its),plot(ma,'g')
  drawnow
  
  if ~strcmp(ScanType,'Cine')
    ma=ma-min(ma);
    mar=round(ma);
    for j=1:size(resp,2)
      resp(:,j)=circshift(resporig(:,j),[-mar(j), 1]);
      resp((end-mar(j)+1):end,j)=resporig(end,j);
    end
    
    resp=(resp-resp*pinv(tempnom)*tempnom);
    for j=1:size(resp,2)
      resp((end-mar(j)+1):end,j)=resporig(1:mar(j),j);
      resp(:,j)=circshift(resp(:,j),[mar(j), 1]);
    end
  end
  
end

Rsig=round((10*ma)+1);
[Rn,~]=hist(Rsig,1:max(Rsig));
Rn=cumsum(Rn)/sum(Rn);
Rsig=Rn(Rsig);

[~,Ridx]=histc(Rsig,linspace(0,1,rbins+1));
Ridx(Ridx==rbins+1)=rbins;
[~,~,Ridx]=unique(Ridx.');

