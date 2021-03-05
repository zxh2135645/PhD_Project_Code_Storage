U = U(:);

morozov=msdev^2*numel(kspace_data)

mov = @(x)reshape(x,Ny,Nx,Nz,[]);

wlev=4; %wavelet level
wname='sym4'; %wavelet name
padamt=ceil([Ny Nx Nz]/2^wlev)*2^wlev-[Ny Nx Nz]; %each dimension needs to be divisible by 2^wlev (16 for wlev=4)
if sum(padamt)~=0 %if needs padding
  if iscartesian
    pad=@(x)padarray(fftshift(fftshift(mov(x),3),1),[padamt 0],0,'post');
    crop=@(x)ifftshift(ifftshift(x(1:Ny,1:Nx,1:Nz,:),3),1);
  else
    pad=@(x)padarray(fftshift(mov(x),3),[padamt 0],0,'post');
    crop=@(x)ifftshift(x(1:Ny,1:Nx,1:Nz,:),3);
  end
  W=@(x)wave3d(pad(reshape(x,[],L)*Wt),wlev,wname); %could incorporate gpuArray (may be memory concerns)
  Wh=@(x)vec(reshape(crop(iwave3d(x)),[],size(Wt,2))*Wt');
else
  W=@(x)wave3d(mov(reshape(x,[],L)*Wt),wlev,wname); %could incorporate gpuArray (may be memory concerns)
  Wh=@(x)vec(reshape(iwave3d(x),[],size(Wt,2))*Wt');
end
WhW=@(x)vec(reshape(x,[],L)*WtWh);

figmake=@(WU)wavelet_figmake(WU);

WU=W(U);
temp=wavelet_collect(WU);
lambda=2*msdev^2/mean(abs(temp(:)-median(real(temp(:)))-1i*median(imag(temp(:)))))
Y=wavelet_applyfun(@(x)zeros(size(x)),WU); %Y=zeros(size(WU));
alpha=1;
figure,hist(vec(abs(temp)),1000);
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
%   figure(100),imshow([im0, log(figmake(WU))],[]),drawnow;
%   figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),drawnow;
  Z=wavelet_applyop(@plus,WU,wavelet_applyop(@rdivide,Y,rho)); %Z=WU+Y/rho;
  Z=wavelet_applyop(@times,wavelet_applyfun(@(x)sign(x),Z),...
    wavelet_applyfun(@(x)max(abs(x)-alpha,0),Z)); %Z=sign(Z).*max(abs(Z)-alpha,0);
  
  Y=wavelet_applyop(@plus,Y,wavelet_applyop(@times,wavelet_applyop(@minus,WU,Z),rho)); %Y=Y+rho*(WU-Z);
  
  step = 1.1;
  alpha=alpha/step;
  rho=lambda/alpha;
  
  Uold=U;
  switch Trajectory
    case 'Radial'
      U=pcg2(@(x)AhA_ps3D(x,st,SEs,Phi_rt,sp)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(wavelet_applyop(@minus,Z,wavelet_applyop(@rdivide,Y,rho))),[],median([5 it 10]),[],[],U); %Z-Y/rho
    case 'Linogram'
      U=pcg2(@(x)AhA_ps_lin(x,st,SEs,Phi_rt,A,Ah,thetas)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(wavelet_applyop(@minus,Z,wavelet_applyop(@rdivide,Y,rho))),[],median([5 it 10]),[],[],U);
    case 'Cartesian'
      U=pcg2(@(x)AhA_ps_cart(x,st,SEs,Phi_rt,Ny,Nx,Nz)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(wavelet_applyop(@minus,Z,wavelet_applyop(@rdivide,Y,rho))),[],median([5 it 10]),[],[],U);
  end
  eps = norm(U(:)-Uold(:))/norm(Uold(:))
  clear Uold;
  if (eps < 1e-3) && (eps ~= 0)
    break;
  end
end
toc;
