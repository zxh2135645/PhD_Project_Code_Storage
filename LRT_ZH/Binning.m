%% preparation
clear temp* Bin* resp* card*
%%
temp_Phi_rt_small = Phi_rt_small_init(:,1:params.NEco:end);
% temp_Nm = params.NRepMeas %total measurements
temp_Nm = 1;
temp_NLPM = (size(kspace_data,1)+cutoff)/params.NEco/SGblock/temp_Nm % Nav lines per measurement
temp_NLPM_2 = params.lSegments*Nz/SGblock;
temp_Lin_check = temp_NLPM - temp_NLPM_2


temp_im = reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L)*temp_Phi_rt_small(:,1),Ny,Nx,[]);
temp_im = abs(temp_im(:,:,1,1)/max(temp_im(:)));

%Draw resp roi
figure('Name','resp roi'),imshow(temp_im);
roi_1 = impoly;
temp_resp_mask = createMask(roi_1);

%Draw card roi
% figure('Name','card roi'),imshow(temp_im);
% roi_2 = impoly;
% temp_card_mask = createMask(roi_2);
% temp_card_mask = ifftshift(temp_card_mask,1);
% 
% temp_card_mask = repmat(temp_card_mask, 1,1, Nz);
% temp_card_mask = vec(temp_card_mask);
%% Filter design

fs=1/(params.lEchoSpacing*SGblock);

% L_filter = designfilt('lowpassfir', ...
%         'PassbandFrequency',30/60, ...
%         'StopbandFrequency',50/60, ...
%         'SampleRate',fs);
% 
% B_filter = designfilt('bandpassfir','FilterOrder',10, ...
%         'CutoffFrequency1',5/60, ...
%         'CutoffFrequency2',30/60, ...
%         'SampleRate',fs);


%% 

Ridx = [];
Hidx = [];
its=35;
ccL = 1;
temp_step = 0;

%%
for temp_step = 0:(temp_Nm-1)
    %%
    temp_step
    if temp_step == 0
%         clear Bin* resp* card*
        Bin_start = 1;
        Bin_cutoff = cutoff/params.NEco/SGblock;
        Bin_end = Bin_start + temp_NLPM -1 - Bin_cutoff;
    else
        
        Bin_start = temp_step * temp_NLPM  - Bin_cutoff +1;
        Bin_end = Bin_start + temp_NLPM -1;
    end
   
    resp_Phi_rt=sgolayfilt(double(temp_Phi_rt_small(:,Bin_start:Bin_end).'),0,5).';

    resp_recon = abs(reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L)*resp_Phi_rt,Ny,Nx,[]));
    resp_mask_rep = repmat(temp_resp_mask,[1,1,size(resp_recon,3)]);
    resp_sig = resp_recon(resp_mask_rep);
    %%
    resp_input = reshape(resp_sig,sum(temp_resp_mask(:)),[]);
    resp_input = abs(resp_input/max(resp_input(:)));
    
    q = 10;
    resp_Mdl = rica(resp_input',q,'NonGaussianityIndicator',ones(q,1));
    resp_ICA = transform(resp_Mdl,resp_input');

    %
    figure(102);
    set(gcf,'Name','ICA results');
    for n=1:q
     figure(102); subplot(q+1,1,n); plot(resp_ICA(:,n));
    end

    %% Spectrum analysis
    
%     fs=1/(params.lEchoSpacing*SGblock);
    df=fs/size(resp_Phi_rt,2);
    winhp=floor((20/60)/df); %highest HR: 30 bpm
    winlp=floor((5/60)/df); %lowest HR: 10 bpm
    l_hwindow = size(resp_Phi_rt,2);
    hwindow=zeros(size(resp_Phi_rt,2),1);
    hwindow(winlp:winhp) = 1;
    hwindow(l_hwindow-winhp+2:l_hwindow-winlp+2) = 1;
    
    for p = 1:q
        resp_x = fft(resp_ICA(:,p));
        resp_x(1) = 0;
        resp_q(p) = norm(resp_x.*hwindow)/norm(resp_x);
    end
    
    resp_Np = find(resp_q == max(resp_q))
%     resp_Np = input(sprintf('number of q: '));
%     close(102);
    %%  
    resp_Z = sgolayfilt(resp_ICA(:,resp_Np),0,7);
    
    if temp_step ~= 0
        close(103);
    end
    figure(103); set(gcf,'Name','resp_Bin');
    figure(103), subplot(3,1,1), plot(resp_Z);
    
%     resp_Z= filtfilt(L_filter,resp_Z);  
%     resp_Z= filtfilt(B_filter,resp_Z); 
    
    resp_Z = fft(resp_Z/sqrt(numel(resp_Z)));
    figure(105),subplot(3,1,1),plot(abs(resp_Z));
    
    resp_Z = resp_Z.*hwindow;
    figure(105), subplot(3,1,2), plot(abs(resp_Z));
 
    fprintf('Delete resp frequences.\n')
    keyboard;
    
%     resp_Z([54, 6590]) = 0;
    figure(105),subplot(3,1,3),plot(abs(resp_Z(2:end)));
    
    resp_Z=real(ifft(resp_Z));
    figure(103), subplot(3,1,2), plot(resp_Z);
  
    max_Z = mean(findpeaks(resp_Z));
    min_Z = (-1)*mean(findpeaks(resp_Z.*(-1)));
    range_Z = max_Z-min_Z;
    resp_Z = ((resp_Z-min_Z)/range_Z)*rbins;
%     resp_Z = ((resp_Z-min(resp_Z))/range(resp_Z))*(rbins-1)+1;
    
    
%     [Zn,~]=hist(resp_Z*10,1:round(max(resp_Z)*10));
%     Zn=cumsum(Zn)/sum(Zn);
%     resp_Z=Zn(round(resp_Z*10));
%     resp_Z=ceil(resp_Z*rbins);
%     resp_Z(resp_Z<1)=1;
%     resp_Z(resp_Z>rbins)=rbins;
%     resp_Z=resp_Z(:);
%  
    resp_Z=ceil(resp_Z);
    resp_Z(resp_Z<1)=1;
    resp_Z(resp_Z>rbins)=rbins;
    resp_Z=resp_Z(:);
    figure(103), subplot(3,1,3), plot(resp_Z);
%     
    
    Ridx(Bin_start:Bin_end) = resp_Z;
    %% Cardiac binning preparation
    
    fs=1/(params.lEchoSpacing*SGblock);
    
    df=fs/size(resp_Phi_rt,2);
    winhp=floor((150/60)/df); %highest HR: 150 bpm
    winlp=floor((80/60)/df); %lowest HR: 40 bpm
    l_hwindow = size(resp_Phi_rt,2);
    hwindow=zeros(size(resp_Phi_rt,2),1);
    hwindow(winlp:winhp) = 1;
    hwindow(l_hwindow-winhp+2:l_hwindow-winlp+2) = 1;
    
    winwidth=0.035/(params.lEchoSpacing*SGblock); %35 ms window width
    winwidth=ceil((winwidth-1)/2)*2+1; %make odd
    
   %% ICA cardiac binning
   
    for p = 1:q
        card_x = fft(resp_ICA(:,p));
        card_x(1) = 0;
        card_q(p) = norm(card_x.*hwindow)/norm(card_x);
    end
    
    card_Np = find(card_q == max(card_q))
    
    card_Z = sgolayfilt(resp_ICA(:,card_Np),0,winwidth);
    
    figure(104); set(gcf,'Name','card_Bin');
    figure(104), subplot(3,1,1), plot(card_Z(1:2000));
    
    card_Z = fft(card_Z/sqrt(numel(card_Z)));
    figure(106),subplot(3,1,1),plot(abs(card_Z));
    
    card_Z = card_Z.*hwindow;
    figure(106),subplot(3,1,2),plot(abs(card_Z));
    
    RR_int = 1000/((find(card_Z(2:end) == max(card_Z(2:end/2)),1)-1)*df)  %theoretical R-R interval seconds
    cbins = floor(RR_int/1000/params.lEchoSpacing/SGblock/2)*2
  
    
    
    fprintf('Delete resp frequences.\n')
    keyboard;
%     card_Z([319 517 8509 8707]) = 0;
    figure(106),subplot(3,1,3),plot(abs(card_Z));
    card_Z = real(ifft(card_Z));
    
%     card_Z=real(ifft(fft(card_Z).*hwindow));
    figure(104), subplot(3,1,2), plot(card_Z(1:2000));
    
    card_Z=card_Z./abs(hilbert(card_Z));
    keyboard;
    
    card_Z=angle(hilbert(card_Z));
    card_Z=ceil((card_Z/pi+1)*cbins/2);
    card_Z(card_Z==0)=1;
    
%     card_Z=(card_Z-min(card_Z))/range(card_Z)*(cbins-1)+1;
%     [Zn,~]=hist(card_Z*10,1:round(max(card_Z)*10));
%     Zn=cumsum(Zn)/sum(Zn);
%     card_Z=Zn(round(card_Z*10));
%     card_Z=ceil(card_Z*cbins);
%     card_Z(card_Z<1)=1;
%     card_Z(card_Z>cbins)=cbins;
%     card_Z=card_Z(:);
    figure(104), subplot(3,1,3), plot(card_Z(1:2000));
    
    
    
    
%     card_Z=fft(card_Z/sqrt(numel(card_Z))).*hwindow;
%     card_Z=real(ifft(card_Z));
%     card_Z=card_Z./abs(hilbert(card_Z));
%    
%     
%     card_Z=ceil((card_Z+1)*cbins/2/2);
%     card_Z(card_Z<1)=1;
%     card_Z(card_Z>cbins/2)=cbins/2; 
%     figure(104), subplot(3,1,3), plot(card_Z(1:2000));
%     
    
    card_Z_ICA = card_Z(:)';
    
    
    
    %% Phi_rt Cardiac binning

    card_Phi_rt= abs(temp_Phi_rt_small(:,Bin_start:Bin_end));

    
    for idx = 1:L_init
        card_Z = abs(card_Phi_rt(idx,:)).';
        card_Z = fft(card_Z);
        card_Z(1) = 0;
%         figure(105),subplot(8,4,idx),plot((abs(card_Z(1:2000))));
        card_idx(idx) = norm(card_Z.*hwindow)/norm(card_Z);
    end   
    [~,card_sorted] = sort(card_idx,'descend')
    
    alpha=cbins/2; %cbins;

    bestres=inf;
    
    for outerit=1:10
%       outerit
      Phicard=zeros(rbins,ceil(cbins/2),L_init,ccL);

      card_Z=abs(card_Phi_rt(card_sorted(outerit),:)).';
%       card_Z=abs(card_Phi_rt(outerit,:)).';
%       figure(104),subplot(5,3,(outerit-1)*3+1),plot(card_Z(1:1000));
%       figure(104),subplot(5,3,(outerit-1)*3+2),plot(abs(fft(card_Z)));
      
      card_Z=fft(card_Z/sqrt(numel(card_Z))).*hwindow;
      card_Z=real(ifft(card_Z));
      card_Z=card_Z./abs(hilbert(card_Z));
      card_Z=ceil((card_Z+1)*cbins/2/2);
      card_Z(card_Z<1)=1;
      card_Z(card_Z>cbins/2)=cbins/2; 
%       figure(104),subplot(5,3,(outerit-1)*3+3),plot(card_Z(1:1000));
      


      for it=1:its

        for j=1:rbins
          for k=1:ceil(cbins/2)
    %         Phicard(j,k,:,:)=Phi_rt_temp(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
            Phicard(j,k,:,:)= mean(card_Phi_rt(:,(resp_Z==j)&(card_Z==k)),2);
          end
        end
        %switch systole/diastole?
        nucnorm=zeros(2^(rbins-1),1);
        parfor j=1:2^(rbins-1) %never switch first one
          reverse_these=(dec2bin(j-1,rbins)=='1');
          Phicardtemp=Phicard;
          Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
          Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
          Phicardtemp(isnan(Phicardtemp))=0;  
          nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
        end
        [~,j]=min(nucnorm);
        reverse_these=(dec2bin(j-1,rbins)=='1');
        Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);
        
        card_resmat=zeros(cbins/2,numel(card_Z));
        parfor j=1:size(card_Phi_rt,2)
    %       [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_temp(:,j),reshape(reshape(Phicard(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
%           [~,card_Z(j)]=min(sum(abs(bsxfun(@minus,card_Phi_rt(:,j),reshape(Phicard(resp_Z(j),:,:,:),ceil(cbins/2),[]).')).^2));
          card_resmat(:,j) = sum(abs(bsxfun(@minus,card_Phi_rt(:,j),reshape(Phicard(resp_Z(j),:,:,:),ceil(cbins/2),[]).')).^2);
        end
        
        card_resmat=bsxfun(@minus,card_resmat,mean(card_resmat));
        resmat_filt=sgolayfilt(card_resmat.',0,winwidth).';
        [~,card_Z]=min(resmat_filt);
        card_Z=card_Z(:);

        card_Z=fft(card_Z/sqrt(numel(card_Z))).*hwindow;
    %     Z = filtfilt(cFilt,Z/sqrt(numel(Z)));
        if (alpha>0) && (it > 10) && (max(abs(card_Z))>alpha)
          card_Z=sign(card_Z).*max(abs(card_Z)-alpha,0);
        end
        card_Z=real(ifft(card_Z));
        
        if it < its
          card_Z=card_Z./abs(hilbert(card_Z));
          card_Z=ceil((card_Z+1)*cbins/2/2);
          card_Z(card_Z<1)=1;
          card_Z(card_Z>cbins/2)=ceil(cbins/2);
        else
          card_Z=angle(hilbert(card_Z));
          card_Z=ceil((card_Z/pi+1)*cbins/2);
          card_Z(card_Z==0)=1;
        end
        
%         card_Z=card_Z./abs(hilbert(card_Z));
%         card_Z=ceil((card_Z+1)*(cbins/2)/2);
%         card_Z(card_Z<1)=1;
%         card_Z(card_Z>ceil(cbins/2))=ceil(cbins/2);

      end

      for j=1:rbins
        for k=1:ceil(cbins/2)
    %       Phicard(j,k,:,:)=Phi_rt_temp(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
          Phicard(j,k,:,:)=mean(card_Phi_rt(:,(resp_Z==j)&(card_Z==k)),2);
        end
      end
      %switch systole/diastole?
      nucnorm=zeros(2^(rbins-1),1);
      parfor j=1:2^(rbins-1) %never switch first one
        reverse_these=(dec2bin(j-1,rbins)=='1');
        Phicardtemp=Phicard;
        Phicardtemp(reverse_these,:,:,:)=flip(Phicardtemp(reverse_these,:,:,:),2);
        Phicardtemp=permute(Phicardtemp,[2 1 3 4]);
        Phicardtemp(isnan(Phicardtemp))=0;  
        nucnorm(j)=sum(svd(Phicardtemp(:,:).','econ'));
      end
      [~,j]=min(nucnorm);
      reverse_these=(dec2bin(j-1,rbins)=='1');
      Phicard(reverse_these,:,:,:)=flip(Phicard(reverse_these,:,:,:),2);

      % Do different Hilbert processing to find true residual
      parfor j=1:size(card_Phi_rt,2)
    %     [~,Z(j)]=min(sum(abs(bsxfun(@minus,Phi_rt_temp(:,j),reshape(reshape(Phicard(Ridx(j),:,:,:),[],ccL)*curvePhi(Segidx(j),1:ccL).',ceil(cbins/2),[]).')).^2));
        [~,card_Z(j)]=min(sum(abs(bsxfun(@minus,card_Phi_rt(:,j),reshape(Phicard(resp_Z(j),:,:,:),ceil(cbins/2),[]).')).^2));
      end

      card_Z=fft(card_Z/sqrt(numel(card_Z))).*hwindow;
    %   Z = filtfilt(cFilt,Z/sqrt(numel(Z)));
      if (alpha>0)
        card_Z=sign(card_Z).*max(abs(card_Z)-alpha,0);
      end
      card_Z=real(ifft(card_Z));
      card_Z=angle(hilbert(card_Z));
      card_Z=ceil((card_Z/pi+1)*(cbins/2));
      card_Z(card_Z==0)=1;
      cbins2=ceil(cbins/2)*2;

      Phicard=zeros(rbins,cbins2,L_init,ccL);
      for j=1:rbins
        for k=1:cbins2
    %       Phicard(j,k,:,:)=Phi_rt_temp(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
          Phicard(j,k,:,:)= mean(card_Phi_rt(:,(resp_Z==j)&(card_Z==k)),2);
        end
      end

      res=zeros(size(card_Z));
      parfor j=1:size(card_Phi_rt,2);
    %     res(j)=norm(Phi_rt_temp(:,j)-squeeze(Phicard(Ridx(j),Z(j),:,:))*curvePhi(Segidx(j),1:ccL).');
        res(j)=norm(card_Phi_rt(:,j)-squeeze(Phicard(resp_Z(j),card_Z(j),:,:)));
      end
      res=norm(res)

      if res<bestres
        Phicardbest=Phicard;
        Hidx(Bin_start:Bin_end) = card_Z;
        bestres=res
        final_iter = card_sorted(outerit)
      end
    end
%}
end


%%
figure,plot(card_Z_ICA(1:1000));hold on;plot(Hidx(1:1000));

size(Ridx)
size(Hidx)

save_Ridx = Ridx;
save_Hidx = Hidx;
figure,subplot(2,1,1),plot(Ridx(1:2000));subplot(2,1,2),plot(Hidx(1:2000));

figure,hist(diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*SGblock*1000,10)
xlabel('RR interval (ms)')


%% Check binning results
% Phicard=Phicardbest;
% bestres

% try
%   temp1=reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL),[],L_init).',Ny,Nx,[]);
% %   temp1=reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard,[2 1 3 4]),[],ccL)*curvePhi(120,1:ccL).',[],L_init).',Ny,Norig,[]);
%   temp1 = ifftshift(temp1,1);
%   implay(2*abs(temp1)/max(abs(temp1(:))))
% catch
%   temp1=pinv(Phi_rt_small.')*nav_data(:,:);
%   temp1=temp1.'*reshape(reshape(Phicard,[],ccL)*curvePhi(end,1:ccL).',[],L_init).';
%   imagesc(abs(temp1))
% end
  
temp_Phi_rt_small_2 = Phi_rt_small_init(:,8:params.NEco:end);


for j=1:rbins
  Phiresp(j,:,:)=mean(temp_Phi_rt_small(:,Ridx==j),2);
end 

temp=reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(Phiresp,[],cL),[],L_init).',Ny,Norig,[]);

h = implay(abs(temp)/max(abs(temp(:))));
set(h.Parent,'Name','resp_binning');


Phicard=zeros(rbins,cbins,L_init,ccL);
for j=1:rbins
    for k=1:cbins
    %     Phicard(j,k,:,:)=Phi_rt_temp(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
      Phicard(j,k,:,:)= mean(temp_Phi_rt_small_2(:,(resp_Z==j)&(card_Z==k)),2);
    end
end


temp2=reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard(1,:,:,:),[2 1 3 4]),[],ccL),[],L_init).',Ny,Nx,[]);
% temp2 = ifftshift(temp1,1);
h = implay(abs(temp2)/max(abs(temp2(:))));
set(h.Parent,'Name','card_binning');


%%

    Ridx = interp1(1:params.NEco:(params.NEco*length(Ridx)),Ridx,1:(params.NEco*length(Ridx)),'previous','extrap');
    Hidx = interp1(1:params.NEco:(params.NEco*length(Hidx)),Hidx,1:(params.NEco*length(Hidx)),'previous','extrap');
    
%     Ridx = reshape(repmat(Ridx,1,8).',[],1);
    Segidx = Ridx;
    Segidx(:) = 1;
%     Hidx = reshape(repmat(Hidx,1,8).',[],1);
    wall_clock = vec(repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco])).';
    
    
    
    

    
    
    
    