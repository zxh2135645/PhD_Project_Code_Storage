%%
load AC_recon Gr Phi U L Wt sizes ScanType Ny Nx dispim vec iscartesian N ovs params Nseg alpha0_deg mainpath file_string curveU cL

if ~iscartesian
  recon = Gr\reshape(Phi(:,:,1,1,ceil(end/2)),L,[]);
else
  recon = Gr\reshape(Phi(:,1:2:end,1,1,ceil(end/2)),L,[]);
end
recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
recon=reshape(recon,N-ovs,N-ovs,[],1);
recon=recon(41:105,43:107,:);

% recon=recon(39:44,39:44,:);
[transf,~,~]=svd(curveU(1:2:end,:),'econ');

%%
recosize=size(recon);
imr = @(x)reshape(x,recosize(1:2));
fullrep = @(x)repmat(x,[1 1 recosize(3)]);

alpha = params.adFlipAngleDegree*pi/180;
e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) ./ (1-cos(alpha)*e);
n = repmat(reshape(1:2:Nseg,1,1,[]),[recosize(1:2) 1]);
Sint = @(A,e,alpha,B)fullrep(A .* Mss(e,alpha)) .* (1 + fullrep(B-1).*fullrep(e*cos(alpha)).^(n-1)) * sin(alpha);
S = @(A,R1,alpha,B)reshape(reshape(Sint(A,e(R1),alpha,B),[],recosize(3))*transf,recosize(1),recosize(2),[]);

ppinv=@(x,y)sum(y.*conj(x),2)./sum(abs(x).^2,2);

if alpha0_deg == 90
  minB = 0;
  maxB = .5;
elseif alpha0_deg == 180
  minB = -1;
  maxB = 1;
end

% y=double(recon./repmat(recon(:,:,end),[1 1 recosize(3)]));
initB = min(max(minB,real(recon(:,:,end)./recon(:,:,1))),maxB);
y = reshape(reshape(recon/max(vec(recon(:,:,end))),[],recosize(3))*transf,recosize(1),recosize(2),[]);

Avp = @(R1,alpha,B)imr(ppinv(reshape(S(ones(recosize(1:2)),R1,alpha,B),[],cL),reshape(y,[],cL)));

cost=@(x)abs(vec(S(Avp(imr(x(1+(1:prod(recosize(1:2))))),x(1),imr(x(prod(recosize(1:2))+1+(1:prod(recosize(1:2)))))),...
  imr(x(1+(1:prod(recosize(1:2))))),x(1),imr(x(prod(recosize(1:2))+1+(1:prod(recosize(1:2))))))-y));

opts=[];
opts.MaxFunEvals = 1e6;
opts.MaxIter = 500;
tic;
tempfit=lsqnonlin(cost,cat(1,alpha/2,vec(ones(recosize(1:2))),initB(:)).',cat(1,0,1/3*vec(ones(recosize(1:2))),minB*vec(ones(recosize(1:2)))).',cat(1,alpha,1/.01*vec(ones(recosize(1:2))),maxB*vec(ones(recosize(1:2)))).',opts);
toc;

tempfit(1)*180/pi
figure,imagesc(abs(Avp(imr(tempfit(1+(1:prod(recosize(1:2))))),tempfit(1),imr(tempfit(prod(recosize(1:2))+1+(1:prod(recosize(1:2)))))))),colormap gray,caxis([0 max(caxis)])
figure,imagesc(imr(tempfit((recosize(1)*recosize(2))+1+(1:(recosize(1)*recosize(2))))),[-1 1])
figure,imagesc(1./imr(tempfit(1+(1:(recosize(1)*recosize(2))))),[0 2.3])