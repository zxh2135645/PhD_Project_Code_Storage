%v0.2.1

Segidx=mod((1:size(Phi_rt_init,2))-1,params.lSegments/2)+1;
ccL=3;

its=35;

fs=1/(params.lEchoSpacing*2);
df=fs/size(Phi_rt_init,2);
winlp=2*floor((150/60)/df);
windowlp=zeros(size(Phi_rt_init,2),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((40/60)/df);
windowhp=zeros(size(Phi_rt_init,2),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));

alpha=cbins; %cbins;

%%
Phicard=zeros(rbins,ceil(cbins/2),L_init,ccL);
bestres=inf;
for outerit=1:10
  outerit
  Z=ceil(rand(size(Phi_rt_init,2),1)*ceil(cbins/2));
  for it=1:its
    
    for j=1:rbins
      for k=1:ceil(cbins/2)
        Phicard(j,k,:,:)=Phi_rt_init(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
      end
    end
    %switch systole/diastole?
    nucnorm=zeros(2^(rbins-1),1);
    for j=1:2^(rbins-1) %never switch first one
      reverse_these=(dec2bin(j-1,rbins)=='1');
      Phicardtemp=Phicard;
      Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
      Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
      nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
    end
    [~,j]=min(nucnorm);
    reverse_these=(dec2bin(j-1,rbins)=='1');
    Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);
    
    parfor j=1:size(Phi_rt_init,2);
      [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_init(:,j),reshape(reshape(Phicard(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
    end
 
    Z=fft(Z/sqrt(numel(Z))).*hwindow;
    if (alpha>0) && (it > 10)
      Z=sign(Z).*max(abs(Z)-alpha,0);
    end
    Z=real(ifft(Z));
    Z=Z./abs(hilbert(Z));
    Z=ceil((Z+1)*ceil(cbins/2)/2);
    Z(Z<1)=1;
    Z(Z>ceil(cbins/2))=ceil(cbins/2);
    
  end
  
  for j=1:rbins
    for k=1:ceil(cbins/2)
      Phicard(j,k,:,:)=Phi_rt_init(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
    end
  end
  %switch systole/diastole?
  nucnorm=zeros(2^(rbins-1),1);
  for j=1:2^(rbins-1) %never switch first one
    reverse_these=(dec2bin(j-1,rbins)=='1');
    Phicardtemp=Phicard;
    Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
    Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
    nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
  end
  [~,j]=min(nucnorm);
  reverse_these=(dec2bin(j-1,rbins)=='1');
  Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);
  
  res=zeros(size(Z));
  parfor j=1:size(Phi_rt_init,2);
    res(j)=norm(Phi_rt_init(:,j)-squeeze(Phicard(Ridx(j),Z(j),:,:))*curvePhi(Segidx(j),1:ccL).');
  end
  res=norm(res)
  
  if res<bestres
    Phicardbest=Phicard;
    Hidx=Z;
    bestres=res;
  end
end

%%
parfor j=1:size(Phi_rt_init,2);
  [~,Hidx(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_init(:,j),reshape(reshape(Phicardbest(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
end

Hidx=fft(Hidx/sqrt(numel(Hidx))).*hwindow;
if (alpha>0)
  Hidx=sign(Hidx).*max(abs(Hidx)-alpha,0);
end
Hidx=real(ifft(Hidx));
Hidx=angle(hilbert(Hidx));
Hidx=ceil((Hidx/pi+1)*ceil(cbins/2));
Hidx(Hidx==0)=1;
cbins=ceil(cbins/2)*2;

%% View
for j=1:rbins
  for k=1:cbins
    Phicard(j,k,:,:)=Phi_rt_init(:,(Ridx==j)&(Hidx==k))*pinv(curvePhi(Segidx((Ridx==j)&(Hidx==k)),1:ccL).');
  end
end

bestres
try
  temp=dispim(reshape(reshape(U_init,[],L_init)*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL)*curvePhi(120,1:ccL).',[],L_init).',Ny,Nx,[]));
  implay(2*abs(temp)/max(abs(temp(:))))
catch
  temp=pinv(Phi_rt_init.')*nav_data(:,:);
  temp=temp.'*reshape(reshape(Phicard,[],ccL)*curvePhi(end,1:ccL).',[],L_init).';
  imagesc(abs(temp))
end
  HRinterv=5;
figure,hist(diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*2000,10)
xlabel('RR interval (ms)')
HRintv=diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*2000;
figure,plot(diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*2000)
figure;plot(HRintv(1:HRinterv:end))
ylabel('RR interval (ms)')
%%
% temp=dispim(reshape(reshape(U_init,[],L_init)*Phi_rt_init(:,Ridx==1),Ny,Nx,[]));
% implay(abs(temp/max(abs(temp(:)))))