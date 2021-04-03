U = U(:);

morozov=msdev^2*numel(kspace_data)

cropTV=@(x)x(ovs/2+(1:Norig),ovs/2+(1:Norig),:,:);
uncropTV=@(x)padarray(x,[ovs/2, ovs/2, 0, 0]);

mov = @(x)reshape(x,Ny,Nx,[]);
Wf=@(U)cropTV(cat(4,mov(U)-circshift(mov(U),[1 0 0]),mov(U)-circshift(mov(U),[0 1 0])));
W=@(U)Wf(reshape(U,[],L)*Wt);
Whf=@(U)U(:,:,:,1)-circshift(U(:,:,:,1),[-1 0 0])+U(:,:,:,2)-circshift(U(:,:,:,2),[0 -1 0]);
Wh=@(U)vec(reshape(Whf(uncropTV(U)),[],size(Wt,2))*Wt');
WhW=@(U)vec(Whf(uncropTV(Wf(reshape(U,[],L)*WtWh))));
figmake=@(WU)sum(sqrt(sum(abs(WU).^2,4)),3);

MS=@(x)vec(bsxfun(@rdivide,mov(x),sum(abs(SEs).^2,3)));

WU=W(U);
% lambda=morozov/norm(vec(sqrt(sum(abs(WU).^2,4))),1)
temp=sqrt(sum(abs(WU).^2,4));
lambda=2*msdev^2/mean(abs(temp(:)-median(real(temp(:)))-1i*median(imag(temp(:)))))
Y=zeros(size(WU));
% Z = WU;
alpha=max(vec(temp))
figure,hist(vec(temp),1000);
im0 = log(figmake(WU));
figure,imshow(im0,[]);
fprintf('Adjust lambda and alpha if desired, then return.\n')
keyboard;
rho=lambda/alpha;
tic;
for it=2:30
  WU = W(U);
  figure(100),imshow([im0, log(figmake(WU))],[]),drawnow;
  figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),drawnow;
  %   Zold = Z;
  Z=WU+Y/rho;
  Zg = sqrt(sum(abs(Z).^2,4));
  Zg = max(abs(Zg)-alpha,0)./Zg;
  Zg(isnan(Zg))=0;
  Z=repmat(Zg,[1 1 1 2]).*Z;
  Y=Y+rho*(WU-Z);
  %   pr=norm(WU(:)-Z(:));
  %   dr=rho*norm(vec(Wh(Z-Zold)));
  %   step=pr/dr
  clear Zold;
  
  %   if step < 0.01
  %       step = 1
  %   else
  %       step = median([1.25 step 2])
  %   end
  step = 1.25;
  alpha=alpha/step;
  rho=lambda/alpha;
  
  Uold=U;
  U=pcg(@(x)AhA_ps(x,st,SEs,Phi_rt,sp)+rho/2*WhW(x),...
    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),MS,[],U);
  eps = norm(U(:)-Uold(:))/norm(Uold(:))
  clear Uold;
%   if (eps < 1e-3) && (eps ~= 0)
%     break;
%   end
end
toc;