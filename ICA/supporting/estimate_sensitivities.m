%
if Ncoils>16
%new compression
  [~,cards,cardmixer]=svde(reshape(kspace_data,[],Ncoils));
  newcoils=find(sqrt(cumsum(diag(cards).^2))/norm(diag(cards))>.99,1) %keep at least 99% of RMS. You could also manually override the new number of coils here.
  
  kspace_data=reshape(reshape(kspace_data,[],Ncoils)*cardmixer(:,1:newcoils),size(kspace_data(:,:,:,1:newcoils)));
  nav_data=reshape(reshape(nav_data,[],Ncoils)*cardmixer(:,1:newcoils),size(nav_data(:,:,:,1:newcoils)));
  Psi=cardmixer(:,1:newcoils)'*Psi*cardmixer(:,1:newcoils);
  Ncoils=newcoils;
  setup_functions;
end
%backproject
[~,~,c]=svd([vec(phantom(N)) vec(Ainv(A(phantom(N),st),st,1))],'econ');
c=real(c(1)/c(2));
[Phi_rt,~,~]=svde(nav_data(:,:));
Phi_rt1=Phi_rt(:,1)/mean(Phi_rt(:,1));
fbp_data=zeros(trajs,N,Ncoils);
for j=1:trajs
  fbp_data(j,:)=Phi_rt1(j:trajs:end)'*kspace_data(j:trajs:end,:)/norm(Phi_rt1(j:trajs:end))^2;
end
fbp=Ainv(fbp_data,st,c);

composite_fbp=sqrt(sum(abs(fbp).^2,3));
figure,imshow(dispim(composite_fbp),[]);
title('Filtered backprojection'),drawnow

%% compress coils to emphasize cardiac region? faster computation.
if false %isradial
  h=figure;
  imshow(dispim(abs(composite_fbp)),[]);
  title('Select Cardiac ROI')
  card_roi=imrect;
  card_roi=padarray(createMask(card_roi),[ovs/2 ovs/2]);
  close(h)
  drawnow
  
  [~,cards,cardmixer]=svde(reshape(fbp(repmat(card_roi,[1 1 Ncoils])),[],Ncoils));
  newcoils=ceil(Ncoils/2);
  kspace_data=reshape(reshape(kspace_data,[],Ncoils)*cardmixer(:,1:newcoils),size(kspace_data(:,:,:,1:newcoils)));
  nav_data=reshape(reshape(nav_data,[],Ncoils)*cardmixer(:,1:newcoils),size(nav_data(:,:,:,1:newcoils)));
  fbp=reshape(reshape(fbp,[],Ncoils)*cardmixer(:,1:newcoils),size(fbp(:,:,1:newcoils)));
  fbp_data=reshape(reshape(fbp_data,[],Ncoils)*cardmixer(:,1:newcoils),size(fbp_data(:,:,1:newcoils)));
  Psi=cardmixer(:,1:newcoils)'*Psi*cardmixer(:,1:newcoils);
  Ncoils=newcoils;
  setup_functions;
end



%% estimate SEs
[~,~,mixer]=svde(reshape(dispim(fbp),[],Ncoils)); %to emphasize central region
if false %if espirit
  SEs=espirit_agc(fftshift(fftshift(fft2(padarray(reshape(fbp,Ny,Nx,1,[]),[(Norig-ovs)/2 (Norig-ovs)/2 0])),1),2),7);
  % SEs=espirit_agc(fftshift(fftshift(fft2(padarray(reshape(dispim(fbp),Norig,Norig,1,[]),[Norig/2 Norig/2 0 0])),1),2),7);
  SEs=SEs((Norig-ovs)/2+(1:Ny),(Norig-ovs)/2+(1:Nx),:); %.*repmat(k<=Norig,[1 1 Ncoils]);
  SEs_mixed = reshape(reshape(SEs,[],Ncoils)*mixer(:,1),Ny,Nx);
  SEs=repmat(exp(-1i*angle(SEs_mixed)),[1 1 Ncoils]).*squeeze(SEs); %correct phase
  clear SEs_mixed
else
  SEs=squeeze(sensemaps_cbd(fbp));
end
%B1normalization
SE_corr=reshape(sqrt(sum(abs(reshape(SEs,[],Ncoils)*sqrtm(Psi)).^2,2)),Ny,Nx);
SEs=SEs./repmat(SE_corr,[1 1 Ncoils]);
pd=abs(sum(conj(SEs).*fbp,3)./sum(abs(SEs).^2,3));
  pd(pd==0)=min(pd(pd~=0));
  [x,y]=ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1);
  polyterm=[x(:).^2, y(:).^2, x(:).*y(:), x(:), y(:), ones(size(x(:)))];
  [polyterm,~,~]=svd(bsxfun(@times,polyterm,pd(:).^2),'econ');
  w=polyterm'*pd(:);
  y=zeros(size(pd(:)));
  rho=1/max(abs(polyterm*w-pd(:)));
  for it=1:10
    res=polyterm*w-pd(:);
    z=res+y/rho;
    z=sign(z).*max(abs(z)-1/rho,0);
    y=y+rho*(res-z);
    rho=rho*1.25;
    w=polyterm'*(z+pd(:)-y/rho);
  end
  [~,~,v]=svd([(polyterm*w)./pd(:), pd(:)],'econ');
  w=w*v(2)/v(1);
  bias=reshape(polyterm*w,Ny,Nx).\pd.^2;
  bias(1./bias<1/2)=2;
  bias(bias<1/4)=1/4;
  figure,imshow(abs([pd pd./bias]),[])
  SEs=bsxfun(@times,SEs,bias);
%

if false %bias correction (old method)
  pd=abs(sum(conj(SEs).*fbp,3)./sum(abs(SEs).^2,3));
  pd=pd/max(pd(:));
  wd=pd>0.1;
  pd=pd/mean(pd(wd));
  [x,y]=ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1);
  polyterm=[x(:).^2, y(:).^2, x(:).*y(:), x(:), y(:), ones(size(x(:)))].';
  w=pd(wd).'*pinv(polyterm(:,wd));
  bias=min(max(reshape(w*polyterm,size(pd)),0.5),2);
  figure,imshow(abs([pd pd./bias]),[])
  SEs=bsxfun(@times,SEs,bias);
else
  pd=abs(sum(conj(SEs).*fbp,3)./sum(abs(SEs).^2,3));
  pd(pd==0)=min(pd(pd~=0));
  [x,y]=ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1);
  polyterm=[x(:).^2, y(:).^2, x(:).*y(:), x(:), y(:), ones(size(x(:)))];
  [polyterm,~,~]=svd(bsxfun(@times,polyterm,pd(:).^2),'econ');
  w=polyterm'*pd(:);
%   norm(polyterm*w-pd(:),1)
  y=zeros(size(pd(:)));
  rho=1/max(abs(polyterm*w-pd(:)));
  for it=1:10
    res=polyterm*w-pd(:);
    z=res+y/rho;
    z=sign(z).*max(abs(z)-1/rho,0);
    y=y+rho*(res-z);
    rho=rho*1.25;
    w=polyterm'*(z+pd(:)-y/rho);
%     norm(polyterm*w-pd(:),1)
  end
  [~,~,v]=svd([(polyterm*w)./pd(:), pd(:)],'econ');
  w=w*v(2)/v(1);
  bias=reshape(polyterm*w,Ny,Nx).\pd.^2;
  bias(1./bias<1/2)=2;
  bias(bias<1/4)=1/4;
  figure,imshow(abs([pd pd./bias]),[])
  SEs=bsxfun(@times,SEs,bias);
end

if useGPU
  st.FS=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,SEs,true);
end

%% SENSE on composite
composite_init=sum(conj(SEs).*fbp,3)./sum(abs(SEs).^2,3);
composite_init(isnan(composite_init))=0;
composite=prep(pcg(@(x)AhA_sense(x,st,SEs,GhG,1),...
  Ah_sense(fbp_data,st,SEs,1),[],20,M,[],composite_init(:)),st);
figure,imshow(abs([dispim(composite_fbp) dispim(composite_init) dispim(prep(composite,st))]),[])
title('Composite Reconstructions'),drawnow;
clear st GhG
