%v0.2

Segidx=mod((1:size(Phi_rt_init,2))-1,params.lSegments/2)+1;
% cL=5;

its=35;

 fs=1/(params.lEchoSpacing*2);
 df=fs/size(Phi_rt_init,2);
% hwin=2*floor((50/60)/df);%
%hwin=2*floor((150/60)/df);%
% hwindow=zeros(size(Phi_rt_init,2),1);
% hwindow(1:hwin)=hamming(hwin,'periodic');
% hwindow=circshift(hwindow,[-hwin/2, 0]);
%%
%Z=Phi_rt_init(2,:)';
%Z=imag(ifft(fft(Z).*hwindow));
%    Z=Z-min(Z)+1;
    
%%

Phiresp=zeros(rbins,L_init,cL);
bestres=inf;
sect=ceil(180*fs);%sliding window size(sec)
Sectn=ceil(length(Phi_rt_init)/sect);
for tsect=1:Sectn*2-1
    tsect
    if (ceil((tsect-1)*sect/2)+sect)<=length(Phi_rt_init)
    Secwindow=(ceil((tsect-1)*sect/2)+1:ceil((tsect-1)*sect/2)+sect);
    else
    Secwindow=ceil((tsect-1)*sect/2)+1:length(Phi_rt_init);  
    end
    Phi_rt_init_temp=Phi_rt_init(:,Secwindow);
    %fs=1/(params.lEchoSpacing*2);
    df=fs/size(Phi_rt_init_temp,2);
    hwin=2*floor((50/60)/df);%
    hwindow=zeros(size(Phi_rt_init_temp,2),1);
    hwindow(1:hwin)=hamming(hwin,'periodic');
    hwindow=circshift(hwindow,[-hwin/2, 0]);
for outerit=1:5
  outerit
  if tsect==1
      Z=ceil(rand(size(Phi_rt_init_temp,2),1)*rbins);
  else
      lastZ=Z;
      Z=ceil(rand(size(Phi_rt_init_temp,2),1)*rbins);
      Z(1:ceil(length(Secwindow)/2))=lastZ(length(lastZ)-ceil(length(Secwindow)/2)+1:end);
  end
  for it=1:its
    for j=1:rbins
      Phiresp(j,:,:)=Phi_rt_init_temp(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
    end
    
    parfor j=1:size(Phi_rt_init_temp,2);
      [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_init_temp(:,j),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2));
    end
    
    Z=real(ifft(fft(Z).*hwindow));
    Z=Z-min(Z)+1;
    
    [Zn,~]=hist(Z*10,1:round(max(Z)*10));
    Zn=cumsum(Zn)/sum(Zn);
    Z=Zn(round(Z*10));
    Z=ceil(Z*rbins);
    Z(Z<1)=1;
    Z(Z>rbins)=rbins;
    Z=Z(:);
    %figure;plot(Z)
  end
  
  for j=1:rbins
    Phiresp(j,:,:)=Phi_rt_init_temp(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
  end
  
  res=zeros(size(Z));
  parfor j=1:size(Phi_rt_init_temp,2);
    res(j)=norm(Phi_rt_init_temp(:,j)-squeeze(Phiresp(Z(j),:,:))*curvePhi(Segidx(j),:).');
  end
  res=norm(res)
  
  if res<bestres
    Ridx(Secwindow)=Z;
    bestres=res;
  end

end
end

%%
for j=1:rbins
  Phiresp(j,:,:)=Phi_rt_init(:,Ridx==j)*pinv(curvePhi(Segidx(Ridx==j),:).');
end
res=zeros(size(Ridx));
parfor j=1:size(Phi_rt_init,2);
  res(j)=norm(Phi_rt_init(:,j)-squeeze(Phiresp(Ridx(j),:,:))*curvePhi(Segidx(j),:).');
end
bestres
temp=dispim(reshape(reshape(U_init,[],L_init)*reshape(reshape(Phiresp,[],cL)*curvePhi(end,:).',[],L_init).',Ny,Nx,[]));
implay(abs(temp)/max(abs(temp(:))))
figure;plot(diff(find((Ridx(2:end)==rbins)&(diff(Ridx)==1))))
xlabel('Respiratory interval (ms)')
%%
% temp=dispim(reshape(reshape(U_init,[],L_init)*Phi_rt_init(:,Ridx==1),Ny,Nx,[]));
% implay(abs(temp/max(abs(temp(:)))))