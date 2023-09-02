%Orthonormalize U (introduces image weighting into Phi_rt, so that binning
%it is identical to directly binning the real-time images)
if exist('U','var')
  Wti=reshape(U_init,[],L);

  %   Wti=Wti(roi_mask,:);
  Wti=Wti'*Wti;
  Wti = sqrtm(Wti); %actually inverse of Wt
  U_init=vec(reshape(U_init,[],L)/Wti);
  U=U_init;
else
  Wti = inv(S_rt(1:L,1:L)); %untested
  U=reshape(nav_data(:,:).'*pinv(Phi_rt_small)/Wti,[],Ncoils,L);
  U=padarray(complex(imresize(real(U),[Norig Norig],'nearest'),imresize(imag(U),[Norig Norig],'nearest')),[ovs/2 ovs/2 0]);
  U=vec(repmat(reshape(U,Ny,Nx,1,L),[1 1 Nz 1]));
  U_init=U;
end
Phi_rt_init=Wti*Phi_rt_init;
Phi_rt=Phi_rt_init;
Phi_rt_full_init=Wti*Phi_rt_full_init;
Phi_rt_full=Phi_rt_full_init;
Phi_rt_small_init=Wti*Phi_rt_small_init;
Phi_rt_small=Phi_rt_small_init;

% Create mask from Binning.m
temp_Phi_rt_small = Phi_rt_small_init(:,1:params.NEco:end); % By only looking at first echo
% temp_Nm = params.NRepMeas %total measurements
temp_Nm = 1;
temp_NLPM = (size(kspace_data,1)+cutoff)/params.NEco/SGblock/temp_Nm % Nav lines per measurement
temp_NLPM_2 = params.lSegments*Nz/SGblock;
temp_Lin_check = temp_NLPM - temp_NLPM_2 % Why? What does this mean?

dispim2 = @(x,st) mean(fftshift(x(:,:,24:26,:),1),3);
temp_im = reshape(reshape(dispim2(reshape(U_init,Ny,Nx,Nz,[])),[],L)*temp_Phi_rt_small(:,1),Ny,Nx,[]);
temp_im = abs(temp_im(:,:,1,1)/max(temp_im(:)));

%Draw resp roi
figure('Name','resp roi'),imshow(temp_im);
roi_1 = impoly;
temp_resp_mask = createMask(roi_1);

mask=ApplyAdjvolumeMask_LRT(twix_obj);
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

Ridx = [];
Hidx = [];
its=35;
ccL = 1;
temp_step = 0;

%%
for temp_step = 0:(temp_Nm-1)
    %
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
    %
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

    % Spectrum analysis
    
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
    %  
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
end

%% Cardiac binning preparation
    
    %fs=1/(params.lEchoSpacing*SGblock);
    
    df=fs/size(resp_Phi_rt,2);
    winhp=floor((150/60)/df); %highest HR: 150 bpm
    winlp=floor((60/60)/df); %lowest HR: 40 bpm
    l_hwindow = size(resp_Phi_rt,2);
    hwindow=zeros(size(resp_Phi_rt,2),1);
    hwindow(winlp:winhp) = 1;
    hwindow(l_hwindow-winhp+2:l_hwindow-winlp+2) = 1;
    
    winwidth=0.035/(params.lEchoSpacing*SGblock); %35 ms window width
    winwidth=ceil((winwidth-1)/2)*2+1; %make odd
    
  % ICA cardiac binning
   
    for p = 1:q
        card_x = fft(resp_ICA(:,p));
        card_x(1) = 0;
        card_q(p) = norm(card_x.*hwindow)/norm(card_x);
    end
    
    card_Np = find(card_q == max(card_q))
    
    card_Z = sgolayfilt(resp_ICA(:,card_Np),0,winwidth); % smoothing
    
    figure(104); set(gcf,'Name','card_Bin');
    figure(104), subplot(3,1,1), plot(card_Z(1:2000));
    
    card_Z = fft(card_Z/sqrt(numel(card_Z)));
    figure(106),subplot(3,1,1),plot(abs(card_Z));
    
    card_Z = card_Z.*hwindow;
    figure(106),subplot(3,1,2),plot(abs(card_Z));
    
    RR_int = 1000/((find(card_Z(2:end) == max(card_Z(2:end/2)),1)-1)*df)  %theoretical R-R interval seconds
    %cbins = floor(RR_int/1000/params.lEchoSpacing/SGblock/2)*2;
    %cbins = 10;
  
    
    
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
    
    % Phi_rt Cardiac binning

    card_Phi_rt= abs(temp_Phi_rt_small(:,Bin_start:Bin_end));
    
    Hidx(Bin_start:Bin_end) = card_Z_ICA;


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
  
temp_Phi_rt_small_2 = Phi_rt_small_init(:,1:params.NEco:end);

clear Phiresp;
for j=1:rbins
  Phiresp(j,:,:)=mean(temp_Phi_rt_small(:,Ridx==j),2);
end 

temp=reshape(reshape(dispim2(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(Phiresp,[],1),[],L_init).',Ny,Norig,[]);

h = implay(abs(temp)/max(abs(temp(:))));
set(h.Parent,'Name','resp_binning');


Phicard=zeros(rbins,cbins,L_init,ccL);
% card_Z = Z;
for j=1:rbins
    for k=1:cbins
    %     Phicard(j,k,:,:)=Phi_rt_temp(:,(Ridx==j)&(Z==k))*pinv(curvePhi(Segidx((Ridx==j)&(Z==k)),1:ccL).');
      Phicard(j,k,:,:)= mean(temp_Phi_rt_small_2(:,(resp_Z==j)&(card_Z==k)),2);
    end
end


temp2=reshape(reshape(dispim2(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)*reshape(reshape(permute(Phicard(1,:,:,:),[2 1 3 4]),[],ccL),[],L_init).',Ny,Nx,[]);
% temp2 = ifftshift(temp1,1);
h = implay(abs(temp2)/max(abs(temp2(:))));
set(h.Parent,'Name','card_binning');

Ridx_init = Ridx;
Hidx_init = Hidx;
%%
num_time_interval = seg;
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

Ridx = interp1(1:params.NEco:(params.NEco*length(Ridx)),Ridx,1:(params.NEco*length(Ridx)),'previous','extrap');
Hidx = interp1(1:params.NEco:(params.NEco*length(Hidx)),Hidx,1:(params.NEco*length(Hidx)),'previous','extrap');
    
%     Ridx = reshape(repmat(Ridx,1,8).',[],1);
    %Segidx = Ridx;
    %Segidx(:) = 1;
%     Hidx = reshape(repmat(Hidx,1,8).',[],1);
wall_clock = vec(repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco])).';

figure, 
subplot(4,1,1); plot(Hidx(1:10000));
subplot(4,1,2); plot(Ridx(1:10000));
subplot(4,1,3); plot(Segidx(1:10000));
subplot(4,1,4); plot(seg_multiplier(1:10000));