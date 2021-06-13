%v0.2
% =======================================================================>>
% XZ - 1 (Will get averaged T1 weighting)
%indices_temp = 1:1:length(nav_indices);
%Segidx = mod(indices_temp-1,params.lSegments)+1;
%Segidx = vec(repmat(Segidx, [params.NEco, 1]));
%Segidx = Segidx(1:length(nav_indices), 1);
% =========================================================================
% =========================================================================
% XZ - 2 (will have temporal resolved T1 weighting)
num_time_interval = params.Averages / seg;
time_interval_seg = (nav_indices(end)+cutoff+1) / num_time_interval;

seg_multiplier = fix((nav_indices-1+cutoff)/time_interval_seg)+1;

Segidx_temp = mod(nav_indices-1,params.lSegments)+1;
Segidx_temp = Segidx_temp(1:length(Segidx_temp)/params.NEco);

seg_multiplier_temp = [];
for i = 1:num_time_interval
   len_temp = sum(seg_multiplier == i);
   seg_multiplier_temp =  [seg_multiplier_temp, ones(1, len_temp / params.NEco)];
end

Segidx = Segidx_temp + (seg_multiplier_temp - 1).* params.lSegments;
Segidx = vec(repmat(Segidx, [params.NEco, 1]));
% =======================================================================<<

%Segidx = interp1(nav_indices, Segidx(:), 1:Nread, 'nearest', 'extrap');
%Segidx(:) = 1;
%cL=5;

its=35;


% XG
clear Phi_rt_FirstEcho
save_Phi_rt_small_init = Phi_rt_small_init;
Phi_rt_FirstEcho(:,:) = Phi_rt_small_init(:,1:8:end);
% Phi_rt_small_init = sgolayfilt(Phi_rt_FirstEcho.',0,5);
% clear Phi_rt_small_init;
Phi_rt_small_init = Phi_rt_FirstEcho;


% fs=1/(params.lEchoSpacing*SGblock);
% df=fs/size(Phi_rt_small_init,2);
% hwin=2*floor((40/60)/df);
% hwindow=zeros(size(Phi_rt_small_init,2),1);
% hwindow(1:hwin)=hamming(hwin,'periodic');
% hwindow=circshift(hwindow,[-hwin/2, 0]);



% set lower limit
fs=1/(params.lEchoSpacing*SGblock);
df=fs/size(Phi_rt_small_init,2);
winlp=2*floor((30/60)/df); %highest HR: 40 bpm
windowlp=zeros(size(Phi_rt_small_init,2),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((10/60)/df); %lowest HR: 3 bpm
windowhp=zeros(size(Phi_rt_small_init,2),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
%

bpFilt = designfilt('bandpassfir','FilterOrder',20,...
'CutoffFrequency1',10/60,'CutoffFrequency2',40/60,...
'SampleRate',29.2227);

%

winwidth=0.1/(params.lEchoSpacing*SGblock); %100 ms window width
winwidth=ceil((winwidth-1)/2)*2+1; %make odd


%% Original
% What is curvePhi here?
%curvePhi = ceil(rand(size(Phi_rt_small_init,2),1));
%curvePhi(:) = 1;
%Segidx = curvePhi.';

Phiresp=zeros(rbins,L_init,cL);
bestresnorm=inf;

Ridx = [];

A = sgolayfilt(Phi_rt_small_init.',0, 5);

for outerit=1:min(15,L_init)
  outerit
  
%   Z=realify(Phi_rt_small_init(outerit,:)).';
%   XG
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


%   Z=abs(Phi_rt_small_init(outerit,:)).';
%   Z=real(ifft(fft(Z).*hwindow));
%   Z = filtfilt(bpFilt,Z);

% 
%   max_Z = mean(findpeaks(Z));
%   min_Z = (-1)*mean(findpeaks(Z.*(-1)));
% % %   min_Z = min(Z);
%   range_Z = max_Z-min_Z;
% %   
%   Z = ((Z-min_Z)/range_Z)*rbins;
% % 
  Z=ceil(Z);
  Z(Z<1)=1;
  Z(Z>rbins)=rbins;
  Z=Z(:);




  
%   figure, plot(Z(1:100));
  
%   figure,plot(Z(1:1000));
  
  resmat=zeros(rbins,numel(Z));
  for it=1:its
    for j=1:rbins
      Phiresp(j,:,:)=Phi_rt_small_init(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
      % Phiresp(j,:,:)=mean(abs(Phi_rt_small_init(:,Z==j)),2);
      % Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Z==j),2); % XG
    end

    parfor j=1:size(Phi_rt_small_init,2)
      resmat(:,j)=sum(abs(bsxfun(@minus,(Phi_rt_small_init(:,j)),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
%     resmat(:,j)=sum(abs(bsxfun(@minus,abs(Phi_rt_small_init(:,j)),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
    end
    resmat=bsxfun(@minus,resmat,mean(resmat));
    resmat_filt=sgolayfilt(resmat.',0,winwidth).';
    [~,Z]=min(resmat_filt);
    Z=Z(:);
    
    Z=real(ifft(fft(Z).*hwindow));
%     Z = filtfilt(bpFilt,Z);
%     max_Z = mean(findpeaks(Z));
%     min_Z = (-1)*mean(findpeaks(Z.*(-1)));
%     range_Z = max_Z-min_Z;
%     Z=((Z-min(Z))/range(Z))*rbins;    
    Z=Z-min(Z)+1;
    
    [Zn,~]=hist(Z*10,1:round(max(Z)*10));
    Zn=cumsum(Zn)/sum(Zn);
    Z=Zn(round(Z*10));
    Z=ceil(Z*rbins);
    Z(Z<1)=1;
    Z(Z>rbins)=rbins;
    Z=Z(:);
%    
%     Z=ceil(Z);
%     Z(Z<1)=1;
%     Z(Z>rbins)=rbins;
%     Z=Z(:);
%     


    
  end
  


  
%   plot(Z(1:100));
  
%   figure,plot(Z(1:1000));
  
%   resmat=zeros(rbins,numel(Z));
%   for it=1:its
%     for j=1:rbins
% %       Phiresp(j,:,:)=Phi_rt_small_init(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
% %       Phiresp(j,:,:)=mean(abs(Phi_rt_small_init(:,Z==j)),2);
%       Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Z==j),2);
%     end
% 
%     parfor j=1:size(Phi_rt_small_init,2)
%       resmat(:,j)=sum(abs(bsxfun(@minus,(Phi_rt_small_init(:,j)),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
% %       resmat(:,j)=sum(abs(bsxfun(@minus,abs(Phi_rt_small_init(:,j)),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
%     end
%     resmat=bsxfun(@minus,resmat,mean(resmat));
%     resmat_filt=sgolayfilt(resmat.',0,winwidth).';
%     [~,Z]=min(resmat_filt);
%     Z=Z(:);
%     
%     Z=real(ifft(fft(Z).*hwindow));
% %     Z = filtfilt(bpFilt,Z);
%     max_Z = mean(findpeaks(Z));
%   min_Z = (-1)*mean(findpeaks(Z.*(-1)));
%     Z=(Z-min(Z))/range(Z)*rbins;    
% %     Z=Z-min(Z)+1;
%     
% %     [Zn,~]=hist(Z*10,1:round(max(Z)*10));
% %     Zn=cumsum(Zn)/sum(Zn);
% %     Z=Zn(round(Z*10));
% %     Z=ceil(Z*rbins);
% %     Z(Z<1)=1;
% %     Z(Z>rbins)=rbins;
% %     Z=Z(:);
% %     Z=real(ifft(fft(Z).*hwindow));
% %     Z=Z-min(Z);
% % 
%     Z=ceil(Z);
%     Z(Z<1)=1;
%     Z(Z>rbins)=rbins;
%     Z=Z(:);
% 
% 
%     
%   end
  
  for j=1:rbins
    Phiresp(j,:,:)=Phi_rt_small_init(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
    % Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Z==j),2);
  end
  
  res=zeros(size(Z));
  parfor j=1:size(Phi_rt_small_init,2)
    res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phiresp(Z(j),:,:))*curvePhi(Segidx(j),:).');
%     res(j)=norm(Phi_rt_small_init(:,j)-Phiresp(Z(j),:,:).');
  end
  resnorm=norm(res)
  
  if resnorm<bestresnorm
    Ridx=Z;
    bestres=res;
    bestresnorm=resnorm;
    Z_num = outerit;
    Z_ref = Z;
  end
end


% xg
% Z_matrix = [];
% 
% parfor outerit=1:L_init
%     outerit
%     if outerit ~= Z_num
% %       Z=realify(Phi_rt_small_init(outerit,:)).';
%       Z=abs(Phi_rt_small_init(outerit,:)).';
%       Z=real(ifft(fft(Z).*hwindow));
% %       Z=(Z-min(Z))/range(Z)*(rbins-1)+1;
% %       [Zn,~]=hist(Z*10,1:round(max(Z)*10));
% %       Zn=cumsum(Zn)/sum(Zn);
% %       Z=Zn(round(Z*10));
% %       Z=ceil(Z*rbins);
% %       Z(Z<1)=1;
% %       Z(Z>rbins)=rbins;
% %       Z=Z(:);
%       [c,lags] = xcorr(Z_ref,Z);
%       c_max_ind = find(~(c-max(c)));
%       Z = circshift(Z,lags(c_max_ind));
%       Z_matrix(outerit,:) = Z;
%     end
%     if outerit == Z_num
%       Z_matrix(outerit,:) = Z_ref;
%     end
% end
% 
% Z_final = mean(Z_matrix);
% Z_final = smooth(Z_final,10,'lowess');

% Ridx(Z_final<3.5) = floor(Z_final(Z_final<3.5));
% Ridx(Z_final>=3.5) = ceil(Z_final(Z_final>=3.5));

for j=1:rbins
    Phiresp(j,:,:)=Phi_rt_small_init(:,Ridx==j)*pinv(curvePhi(Segidx(Ridx==j),:).');
    % Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Ridx==j),2);
  end
  
  res=zeros(size(Z));
  parfor j=1:size(Phi_rt_small_init,2)
    res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phiresp(Ridx(j),:,:))*curvePhi(Segidx(j),:).');
%     res(j)=norm(Phi_rt_small_init(:,j)-Phiresp(Ridx(j),:,:).');
  end
  resnorm=norm(res)
% 
% Ridx = Z_ref;
disp('Finished resp binning.')
%% XG

% curvePhi = ceil(rand(size(Phi_rt_small_init,2),1));
% curvePhi(:) = 1;
% Segidx = curvePhi.';
% 
% Phiresp=zeros(rbins,L_init,cL);
% bestresnorm=inf;
% 
% Ridx = [];
% 
% filtered_Phi = filtfilt(bpFilt,Phi_rt_small_init.');
% 
% for outerit=1:min(10,L_init)
%   outerit
%   
% %   XG
% %   Z=realify(filtered_Phi(:,outerit));
%   Z=abs(filtered_Phi(:,outerit));
%   max_Z = mean(findpeaks(Z));
%   min_Z = (-1)*mean(findpeaks(Z.*(-1)));
% %   min_Z = min(Z);
%   range_Z = max_Z-min_Z;
%   Z = ((Z-min_Z)/range_Z)*rbins;
%   Z = ceil(Z);
%   
%   Z(Z<1) = 1;
%   Z(Z>rbins) = rbins;
%   
%   Z1 = Z;
%   
%   resmat=zeros(rbins,numel(Z));
%   for it=1:its
%     for j=1:rbins
% %       Phiresp(j,:,:)=Phi_rt_small_init(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
%       Phiresp(j,:,:)=mean(abs(filtered_Phi(Z==j,:)));
% %       Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Z==j),2);
%     end
% 
%     parfor j=1:size(Phi_rt_small_init,2)
%       resmat(:,j)=sum(abs(bsxfun(@minus,abs(filtered_Phi(j,:).'),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
% %       resmat(:,j)=sum(abs(bsxfun(@minus,Phi_rt_small_init(:,j),reshape(reshape(Phiresp,[],cL)*curvePhi(Segidx(j),:).',rbins,[]).')).^2);
%     end
%     resmat=bsxfun(@minus,resmat,mean(resmat));
%     resmat_filt=sgolayfilt(resmat.',0,winwidth).';
%     [B,Z]=min(resmat_filt);
% %     [B,Z]=min(resmat);
%     Z=Z(:);
%    
%     Z = filtfilt(bpFilt,Z);
% %     Z=Z-min(Z)+1;
% %     
% %     [Zn,~]=hist(Z*10,1:round(max(Z)*10));
% %     Zn=cumsum(Zn)/sum(Zn);
% %     Z=Zn(round(Z*10));
% %     Z=ceil(Z*rbins);
% %     Z(Z<1)=1;
% %     Z(Z>rbins)=rbins;
% %     Z=Z(:);
%     max_Z = mean(findpeaks(Z));
%     min_Z = (-1)*mean(findpeaks(Z.*(-1)));
% %   min_Z = min(Z);
%     range_Z = max_Z-min_Z;
%     Z = ((Z-min_Z)/range_Z)*rbins+0.5;
%     Z = round(Z);
%   
%     Z(Z<1) = 1;
%     Z(Z>rbins) = rbins;
%     
%     
%     
%   end
%   
%   
%   for j=1:rbins
% %     Phiresp(j,:,:)=Phi_rt_small_init(:,Z==j)*pinv(curvePhi(Segidx(Z==j),:).');
% %     Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Z==j),2);
%     Phiresp(j,:,:)=mean(filtered_Phi(Z==j,:));
%   end
%   
%   res=zeros(size(Z));
%   parfor j=1:size(Phi_rt_small_init,2)
%     res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phiresp(Z(j),:,:))*curvePhi(Segidx(j),:).');
% %     res(j)=norm(Phi_rt_small_init(:,j)-Phiresp(Z(j),:,:).');
%   end
%   resnorm=norm(res)
%   
%   if resnorm<bestresnorm
%     Ridx=Z;
%     bestres=res;
%     bestresnorm=resnorm;
%     Z_num = outerit;
%     Z_ref = Z;
%   end
% end
% 
% 
% % xg
% % Z_matrix = [];
% % 
% % parfor outerit=1:L_init
% %     outerit
% %     if outerit ~= Z_num
% % %       Z=realify(Phi_rt_small_init(outerit,:)).';
% %       Z=abs(Phi_rt_small_init(outerit,:)).';
% %       Z=real(ifft(fft(Z).*hwindow));
% % %       Z=(Z-min(Z))/range(Z)*(rbins-1)+1;
% % %       [Zn,~]=hist(Z*10,1:round(max(Z)*10));
% % %       Zn=cumsum(Zn)/sum(Zn);
% % %       Z=Zn(round(Z*10));
% % %       Z=ceil(Z*rbins);
% % %       Z(Z<1)=1;
% % %       Z(Z>rbins)=rbins;
% % %       Z=Z(:);
% %       [c,lags] = xcorr(Z_ref,Z);
% %       c_max_ind = find(~(c-max(c)));
% %       Z = circshift(Z,lags(c_max_ind));
% %       Z_matrix(outerit,:) = Z;
% %     end
% %     if outerit == Z_num
% %       Z_matrix(outerit,:) = Z_ref;
% %     end
% % end
% % 
% % Z_final = mean(Z_matrix);
% % Z_final = smooth(Z_final,10,'lowess');
% 
% % Ridx(Z_final<3.5) = floor(Z_final(Z_final<3.5));
% % Ridx(Z_final>=3.5) = ceil(Z_final(Z_final>=3.5));
% 
% for j=1:rbins
% %     Phiresp(j,:,:)=Phi_rt_small_init(:,Ridx==j)*pinv(curvePhi(Segidx(Ridx==j),:).');
%     Phiresp(j,:,:)=mean(Phi_rt_small_init(:,Ridx==j),2);
%   end
%   
%   res=zeros(size(Z));
%   parfor j=1:size(Phi_rt_small_init,2)
%     res(j)=norm(Phi_rt_small_init(:,j)-squeeze(Phiresp(Ridx(j),:,:))*curvePhi(Segidx(j),:).');
% %     res(j)=norm(Phi_rt_small_init(:,j)-Phiresp(Ridx(j),:,:).');
%   end
%   resnorm=norm(res)
% % 
% % Ridx = Z_ref;
% disp('Finished resp binning.')





%%
for j=1:rbins
  Phiresp(j,:,:)=Phi_rt_small_init(:,Ridx==j)*pinv(curvePhi(Segidx(Ridx==j),:).');
end 
bestresnorm

temp=ifftshift(reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(Phiresp,[],cL)*curvePhi(end,:).',[],L_init).',Ny,Norig,[]),1);

implay(abs(temp)/max(abs(temp(:))))



%%
%gray_im = fftshift(reshape(U_reshap*Phi_small,Nz,Nx,[]),1);






















