U = U(:);

morozov=msdev^2*numel(kspace_data)

if ~iscartesian && ~strcmp(ScanType,'T2prep')
  cropTV=@(x)x(ovs/2+(1:Norig),ovs/2+(1:Norig),:,:,:);
  uncropTV=@(x)padarray(x,[ovs/2, ovs/2, 0, 0, 0]);
else
  cropTV=@(x)x;
  uncropTV=@(x)x;
end
mov = @(x)reshape(x,Ny,Nx,Nz,[]);
Wf=@(U)cropTV(cat(5,mov(U)-circshift(mov(U),[1 0 0 0]),mov(U)-circshift(mov(U),[0 1 0 0]),mov(U)-circshift(mov(U),[0 0 1 0])));
W=@(U)Wf(reshape(U,[],L)*Wt);
Whf=@(U)U(:,:,:,:,1)-circshift(U(:,:,:,:,1),[-1 0 0 0])+U(:,:,:,:,2)-circshift(U(:,:,:,:,2),[0 -1 0 0])+U(:,:,:,:,3)-circshift(U(:,:,:,:,3),[0 0 -1 0]);
Wh=@(U)vec(reshape(Whf(uncropTV(U)),[],size(Wt,2))*Wt');
WhW=@(U)vec(Whf(uncropTV(Wf(reshape(U,[],L)*WtWh))));
figmake=@(WU)sum(sqrt(sum(abs(WU(:,:,1,:,:)).^2,5)),4);

WU=W(U);
temp=sqrt(sum(abs(WU).^2,5));
lambda=2*msdev^2/mean(abs(temp(:)-median(real(temp(:)))-1i*median(imag(temp(:)))))
Y=zeros(size(WU));
alpha=max(vec(temp))
figure,hist(vec(temp),1000);
im0 = log(figmake(WU));
figure,imshow(im0,[]);
fprintf('Adjust lambda and alpha if desired, then return.\n')

keyboard;

if useGPU
  delete(gcp('nocreate')); %delete current parpool
  parpool('local',8); %open parpool with 8 workers
  gpuWorkerReset(); %reassign workers
end

rho=lambda/alpha;
tic;
for it=2:30
  WU = W(U);
  figure(100),imshow([im0, log(figmake(WU))],[]),drawnow;
  figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),drawnow;
  Z=WU+Y/rho;
  Zg = sqrt(sum(abs(Z).^2,5));
  Zg = max(abs(Zg)-alpha,0)./Zg;
  Zg(isnan(Zg))=0;
  Z=repmat(Zg,[1 1 1 1 3]).*Z;
  Y=Y+rho*(WU-Z);
  
  step = 1.25;
  alpha=alpha/step;
  rho=lambda/alpha;
  
  Uold=U;
  switch Trajectory
    case 'Radial'
      U=pcg2(@(x)AhA_ps3D(x,st,SEs,Phi_rt,sp)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
    case 'Linogram'
      U=pcg2(@(x)AhA_ps_lin(x,st,SEs,Phi_rt,A,Ah,thetas)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
    case 'Cartesian'
      U=pcg2(@(x)AhA_ps_cart(x,st,SEs,Phi_rt,Ny,Nx,Nz)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
  end
  eps = norm(U(:)-Uold(:))/norm(Uold(:))
  clear Uold;
  if (eps < 1e-3) && (eps ~= 0)
    break;
  end
end
toc;
