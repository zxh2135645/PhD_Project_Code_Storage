%v0.2.1

Segidx(:) = 1;
%Segidx=mod(nav_indices-1,params.lSegments)+1;
ccL=1;

its=35;

fs=1/(params.lEchoSpacing*SGblock);
df=fs/size(Phi_rt_small_init,2);
winlp=2*floor((150/60)/df); %highest HR: 150 bpm
windowlp=zeros(size(Phi_rt_small_init,2),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((60/60)/df); %lowest HR: 40 bpm
windowhp=zeros(size(Phi_rt_small_init,2),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));

alpha=cbins/2; %cbins;


% cFilt = designfilt('bandpassfir','FilterOrder',20,...
% 'CutoffFrequency1',60/60,'CutoffFrequency2',150/60,...
% 'SampleRate',29.2227);


% save_Ridx = Ridx;
% Ridx = Z_ref;


%%
bestres=inf;
for outerit=1:10
  outerit
  Phicard=zeros(rbins,ceil(cbins/2),L_init,ccL);
%   Z=ceil(rand(size(Phi_rt_small_init,2),1)*ceil(cbins/2));
  
  Z=abs(Phi_rt_small_init(outerit,:)).';
  
  Z=real(ifft(fft(Z).*hwindow));

%   Z = sgolayfilt(Z,0,5);
%   Z = filtfilt(bpFilt,Z);

  Z=(Z-min(Z))/range(Z)*(rbins-1)+1;
  [Zn,~]=hist(Z*10,1:round(max(Z)*10));
  Zn=cumsum(Zn)/sum(Zn);
  Z=Zn(round(Z*10));
  Z=ceil(Z*rbins);
  Z(Z<1)=1;
  Z(Z>rbins)=rbins;
  Z=Z(:);

  
  
  for it=1:its
    
    for j=1:rbins
      for k=1:ceil(cbins/2)
%         Phicard(j,k,:,:)=Phi_rt_small_init(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
        Phicard(j,k,:,:)= mean(Phi_rt_small_init(:,(Ridx==j)&(Z==k)),2);
      end
    end
    %switch systole/diastole?
    nucnorm=zeros(2^(rbins-1),1);
    parfor j=1:2^(rbins-1) %never switch first one
      reverse_these=(dec2bin(j-1,rbins)=='1');
      Phicardtemp=Phicard;
      Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
      Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
      Phicardtemp(isnan(Phicardtemp)) = 0;
      Phicardtemp(isinf(Phicardtemp)) = 0;
      nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
    end
    [~,j]=min(nucnorm);
    reverse_these=(dec2bin(j-1,rbins)=='1');
    Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);
    
    parfor j=1:size(Phi_rt_small_init,2);
%       [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_small_init(:,j),reshape(reshape(Phicard(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
      [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_small_init(:,j),reshape(Phicard(Ridx(j),:,:,:),ceil(cbins/2),[]).')).^2));
    end
    
    Z=fft(Z/sqrt(numel(Z))).*hwindow;
%     Z = filtfilt(cFilt,Z/sqrt(numel(Z)));
    if (alpha>0) && (it > 10) && (max(abs(Z))>alpha)
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
%       Phicard(j,k,:,:)=Phi_rt_small_init(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
      Phicard(j,k,:,:)=mean(Phi_rt_small_init(:,(Ridx==j)&(Z==k)),2);
    end
  end
  %switch systole/diastole?
  nucnorm=zeros(2^(rbins-1),1);
  parfor j=1:2^(rbins-1) %never switch first one
    reverse_these=(dec2bin(j-1,rbins)=='1');
    Phicardtemp=Phicard;
    Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
    Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
    nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
  end
  [~,j]=min(nucnorm);
  reverse_these=(dec2bin(j-1,rbins)=='1');
  Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);
  
  % Do different Hilbert processing to find true residual
  parfor j=1:size(Phi_rt_small_init,2)
%     [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_small_init(:,j),reshape(reshape(Phicard(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
    [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_small_init(:,j),reshape(Phicard(Ridx(j),:,:,:),ceil(cbins/2),[]).')).^2));
  end
  
  Z=fft(Z/sqrt(numel(Z))).*hwindow;
%   Z = filtfilt(cFilt,Z/sqrt(numel(Z)));
  if (alpha>0)
    Z=sign(Z).*max(abs(Z)-alpha,0);
  end
  Z=real(ifft(Z));
  Z=angle(hilbert(Z));
  Z=ceil((Z/pi+1)*ceil(cbins/2));
  Z(Z==0)=1;
  cbins2=ceil(cbins/2)*2;
  
  Phicard=zeros(rbins,cbins2,L_init,ccL);
  for j=1:rbins
    for k=1:cbins2
%       Phicard(j,k,:,:)=Phi_rt_small_init(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
      Phicard(j,k,:,:)= mean(Phi_rt_small_init(:,(Ridx==j)&(Z==k)),2);
    end
  end
  
  res=zeros(size(Z));
  parfor j=1:size(Phi_rt_small_init,2);
%     res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phicard(Ridx(j),Z(j),:,:))*curvePhi(Segidx(j),1:ccL).');
    res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phicard(Ridx(j),Z(j),:,:)));
  end
  res=norm(res)
  
  if res<bestres
    Phicardbest=Phicard;
    Hidx=Z;
    bestres=res;
  end
end



%%
Phicard=Phicardbest;
bestres
try
  temp1=reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL)*curvePhi(120,1:ccL).',[],L_init).',Ny,Norig,[]);
  temp1 = ifftshift(temp1,1);
  implay(2*abs(temp1)/max(abs(temp1(:))))
% catch
%   temp1=pinv(Phi_rt_small.')*nav_data(:,:);
%   temp1=temp1.'*reshape(reshape(Phicard,[],ccL)*curvePhi(end,1:ccL).',[],L_init).';
%   imagesc(abs(temp1))
end
  
figure,hist(diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*SGblock*1000,10)
xlabel('RR interval (ms)')
%%
% temp=dispim(reshape(reshape(U_init,[],L_init)*Phi_rt_small_init(:,Ridx==1),Ny,Nx,[]));
% implay(abs(temp/max(abs(temp(:)))))