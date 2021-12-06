% Data binningwith ICA analysis
% Cedars Sinai 
% Randy Yang 2018 Oct
% Hsin-Jung.Yang@cshs.org

%% setup variables
clear masksig masksigC masksigR signalSeg Zseg f
ntemp=ceil(total_time/400);
Phi_rt_temp=sgolayfilt(double(Phi_rt(:,SGblock:SGblock:end).'),0,5).';% only use navdata for gatting
segtemp=floor(size(Phi_rt_temp,2)/ntemp);
cstep=0;
% use the adjustvolume as the cardiac ROI
mask=ApplyAdjvolumeMask_LRT_cart(twix_obj);
systolep=[];
diastolep=[];
systolepd=[];
diastolepd=[];
systoleprom=[];
diastoleprom=[];
systolewidth=[];
diastolwidth=[];
expp=[];
insp=[];
exppd=[];
inspd=[];
Diaframxy=0;
SliceforMotionROI=floor(Nz/3)-1:floor(Nz/3)+1;
Uslice=ifftshift(1:Nz);
dispimICA = @(x) mean(fftshift(fftshift(x(:,:,Uslice(SliceforMotionROI),:),1),3),3);
dispimICA3D = @(x) (fftshift(fftshift(x(:,:,:,:),1),3));

sample_im = dispimICA(reshape(U,Ny,Nx,Nz,[]));
sample_im = abs(sample_im(:,:,:,1)/max(sample_im(:)));

implay(sample_im);
disp('cardiac ROI')
figure(110),imshow(sample_im);
title('Cardiac ROI')
if exist(Croi)
    %load('Croi.mat')
else
    Croi =roipoly; %RY draw ROI
end
    %load('Croi.mat')
disp('cardiac ROI Done')
%  
%  disp('Respiratory ROI')
% %figure,imshow(sample_im);
% title('Respiratory ROI')
% Rroi = roipoly;
 
mask_s1 = mask(:,:,Nz/2);
mask_s1 = flip(mask_s1,2)';
mask_s1 = flip(mask_s1,2);
%%
for ctemp=1:ntemp
    
    %% segmentwised realtime recon to get the ROI of heart in 2D slices
recon = (reshape(reshape(dispimICA(reshape(U,Ny,Nx,Nz,[])),[],L)*Phi_rt_temp(:,(ctemp-1)*segtemp+1:(ctemp)*segtemp),Ny,Nx,[]));
cw=abs(max(recon(:)));
%recon=mean(recon3D(:,:,SliceforCardiacMotion,:),3);
% recon=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,(ctemp-1)*segtemp+1:(ctemp)*segtemp),Ny,Nx,[])));
cstep=cstep+1;
    tempmask=repmat(Croi,[1,1,size(recon,3)]);
    tempmasksig=recon(tempmask);
    masksig((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3),:)=reshape(tempmasksig,sum(Croi(:)),[])';
%      tempmask=repmat(Rroi,[1,1,size(recon,3)]);
%     tempmasksig=recon(tempmask);
%     masksigR((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3),:)=reshape(tempmasksig,sum(Rroi(:)),[])';

%comment for debug
clear Phi_rt_temp tempmask tempmasksig
%% ICA analysis

%input=masksig;
input= masksig((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3),:);
%
plotrange=1:10000;%size(input,1);
TR=params.lEchoSpacing;
Nav_interval=params.lEchoSpacing*2; 

input=abs(input);
input=input/max(input(:));
%nComps=4;
%[icasig, A, W] = fastica(input,'numOfIC', nComps);
q=8;
Mdl = rica(input,q,'NonGaussianityIndicator',ones(q,1));
%
unmixed = transform(Mdl,input);


%%
figure(102);
set(gcf,'Name','ICA results');
for n=1:q

 figure(102); subplot(q+1,1,n); plot((1:length(unmixed(plotrange,n))),unmixed(plotrange,n));  hold on
end

%%
% subplot(q+1,1,n+1); plot(Hidx(plotrange));
% figure(4);
% set(gcf,'Name','ICA results');
% l=[1,4,6];
% for n=1:3
%   subplot(3,1,n); plot(unmixed(plotrange,l(n)));  
%    hold on; subplot(3,1,n); plot(mod(plotrange,(params.lSegments/2))/(params.lSegments/2)*2*std(unmixed(plotrange,l(n)))+mean(unmixed(plotrange,l(n))),'r'); 
% 
% end

%figure;plot(icasig')
%% PCA
%% PCA/SVD decomposition of the data
%[U,S,V]=svd(input,'econ');
%show recovered components
% for iter=1:nComps,
%  C_PCA{iter}=reshape(U(:,iter),[iS iS]);
% end;

%% spectral analysis to find the cardiac and respiratory components
TR=params.lEchoSpacing;

fs = 1/2/TR; 
% n = length(x);          % number of samples
n = size(unmixed,1);
f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range

df=fs/size(unmixed,1);% sample frequency (Hz)
%figure;
q=size(unmixed,2);
%%%%% prepmodulation frequency
fsm=1/(TR*params.lSegments);
%wsm=3*fs;
%%% notch filter
%%
  unmixed_decoupl= unmixed;
for fn=1:floor(max(f)/fsm)
d(fn) = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',fsm*fn-df/0.5,'HalfPowerFrequency2',fsm*fn+df/0.5, ...
               'DesignMethod','butter','SampleRate',fs);
           for p=1:q
            unmixed_decoupl(:,p)= filtfilt(d(fn),unmixed_decoupl(:,p));    
            end
end
%%
figure(103);
set(gcf,'Name','ICA notch filtered results');
for n=1:q
  %figure;imagesc(temp)
  %subplot(q+1,1,n); plot(Nav_interval*plotrange,unmixed(plotrange,n));  
 figure(103); subplot(q+1,1,n); plot((1:length(unmixed_decoupl(plotrange,n))),unmixed_decoupl(plotrange,n));  hold on
 %hold on; subplot(q+1,1,n); plot(mod(plotrange,(params.lSegments/2))/(params.lSegments/2)*2*std(unmixed(plotrange,n))+mean(unmixed(plotrange,n)),'r'); 
end
%%
%%%signal decouple from the Prep pulse and trajectory modulation

figure;
for p=1:q
%subplot(q+1,1,n); plot(Nav_interval*plotrange,unmixed_decoupl(plotrange,n));  
x=abs(unmixed_decoupl(:,p));

y = fft(x);
n = length(x);          % number of samples
harmonicsm=(abs(mod(f,fsm))<df)|(abs(mod(f,fsm))>(fsm-df)); %harmocin frequency for prep modulation
powr(p,:) = abs(y).^2/n;    % power of the DFT
frange=(f>0.001)&(f<5);
figure(120);subplot(q,1,p); plot(f(frange),powr(p,frange));%/sum(power(:,p)));  


Cardiacrange=(f>(45/60))&(f<(200/60)).*~harmonicsm;
respiratoryrange=(f>(5/60))&(f<(1/2)).*~harmonicsm;
AcceptRange=(f>(5/60))&(f<(300/60)).*~harmonicsm;
prepmodscore(p)=sum(powr(p,frange).*harmonicsm(frange))/sum(powr(p,AcceptRange));

cscore(p)=sum(powr(p,frange).*Cardiacrange(frange))/sum(powr(p,frange).*AcceptRange(frange));
rscore(p)=sum(powr(p,frange).*respiratoryrange(frange))/sum(powr(p,frange).*AcceptRange(frange));

end

Cardiacp=find(cscore==max(cscore))
Respp=find(rscore==max(rscore))

figure(150);
subplot(ntemp,1,ctemp); plot(f(frange),powr(Cardiacp,frange));  
 figure(170);
subplot(ntemp,1,ctemp); plot(f(frange),powr(Cardiacp,frange).*Cardiacrange(frange));  
clear powr Spectditrib

%% Cardiac frequency filtering for animals
fs=1/(params.lEchoSpacing*2);
df=fs/size(unmixed,1);
n = size(unmixed,1);          % number of samples
clear f
f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
harmonicsm=((abs(mod(f,fsm))<(df/2))|(abs(mod(f,fsm))>(fsm-(df/2)))); %harmocin frequency for prep modulation
winlp=2*floor((180/60)/df);
windowlp=zeros(size(unmixed,1),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((60/60)/df);
windowhp=zeros(size(unmixed,1),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
%hwindow=hwindow.*~harmonicsm';%demodulate the prepsignal


hf = designfilt('highpassiir','StopbandFrequency',45/60, ...
               'PassbandFrequency',60/60, ...
                'SampleRate',fs);
            Z= filtfilt(hf,unmixed_decoupl(:,Cardiacp));    

figure(160);hold on;subplot(ntemp,1,ctemp); plot(Z);
figure(180);hold on;subplot(ntemp,1,ctemp); plot(unmixed_decoupl(:,Cardiacp));
signalSeg(ctemp,:)=unmixed_decoupl(:,Cardiacp);
Zseg(ctemp,:)=Z;

%%
tm=1:length(Z);
%tm=tm*TR*2;
%%

 hrcutoff=180;
 clear triggerpoint
for n=1:length(hrcutoff);
%%
locsp =[]; 
locsv =[];
locspd =[];
locsvd=[];
%detect peak section by section
peakdetectintT=60;%sec
peakdetectint=floor(peakdetectintT/TR/2);
pstepn=ceil(length(Z)/peakdetectint);
for pseg=1:pstepn
    prange=((pseg-1)*peakdetectint+1:(pseg-1)*peakdetectint+...
        min(peakdetectint,length(Z((pseg-1)*peakdetectint+1:end))));
locsp=[locsp,ampd(Z(prange))+(pseg-1)*peakdetectint];
locsv=[locsv,ampd(-Z(prange))+(pseg-1)*peakdetectint];
locspd=[locspd,ampd(diff(movmean(unmixed_decoupl(prange,Cardiacp),round((1/20)/(2*TR)))))+(pseg-1)*peakdetectint];
locsvd=[locsvd,ampd(diff(movmean(-unmixed_decoupl(prange,Cardiacp),round((1/20)/(2*TR)))))+(pseg-1)*peakdetectint];
%locspvd=ampd(-abs(diff(movmean(-unmixed_decoupl(:,Cardiacp),round((1/20)/(2*TR))))));
end
Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';
recon=abs(dispimICA(reshape(reshape(U,Ny*Nx*Nz,[])*Phi_rt_temp(:,locsp),Ny,Nx,Nz,[])));
implay(abs(recon(:,:,:))/cw);

%
recon=abs(dispimICA(reshape(reshape(U,Ny*Nx*Nz,[])*Phi_rt_temp(:,locspd),Ny,Nx,Nz,[])));
implay(abs(recon(:,:,:))/cw);
%

Peakmean=mean(mean(abs(masksig(locsp,:))))
Vallymean=mean(mean(abs(masksig(locsv,:))))

% determine diastolic for Z
%put it back to the original Philines by multiplying SGblack

if Peakmean>Vallymean

systolep=[systolep,locsv+(cstep-1)*segtemp];
diastolep=[diastolep,locsp+(cstep-1)*segtemp];



else
%peak is systole
systolep=[systolep,locsp+(cstep-1)*segtemp];
diastolep=[diastolep,locsv+(cstep-1)*segtemp];


end
%
% determine diastolic for diffZ
[locpvd,idx]=sort([locsvd,locspd]);
locpvddif=diff(locpvd);
intvv=locpvddif(idx(1:end-1)<=length(locsvd));
intvp=locpvddif(idx(1:end-1)>length(locsvd));

if mean(intvp)>mean(intvv)
%valley is systole

systolepd=[systolepd,locsvd+(cstep-1)*segtemp];
diastolepd=[diastolepd,locspd+(cstep-1)*segtemp];


else
%peak is systole

systolepd=[systolepd,locspd+(cstep-1)*segtemp];
diastolepd=[diastolepd,locsvd+(cstep-1)*segtemp];

end

end
%% Respiratory binning using unmixed signal directly

fs=1/(params.lEchoSpacing*2);
df=fs/size(unmixed,1);
n = size(unmixed,1);          % number of samples
clear f
f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
harmonicsm=((abs(mod(f,fsm))<(df/2))|(abs(mod(f,fsm))>(fsm-(df/2)))); %harmocin frequency for prep modulation
winlp=2*floor((3/60)/df);
windowlp=zeros(size(unmixed,1),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((40/60)/df);
windowhp=zeros(size(unmixed,1),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
%hwindow=hwindow.*~harmonicsm';%demodulate the prepsignal


lf = designfilt('lowpassiir','StopbandFrequency',40/60, ...
               'PassbandFrequency',30/60, ...
                'SampleRate',fs);
            ZR= filtfilt(hf,unmixed_decoupl(:,Respp));    

figure(161);hold on;subplot(ntemp,1,ctemp); plot(ZR);
figure(181);hold on;subplot(ntemp,1,ctemp); plot(unmixed_decoupl(:,Respp));
signalSeg(ctemp,:)=unmixed_decoupl(:,Respp);
ZRseg(ctemp,:)=ZR;

%% detect inspiration and expiration
locsp =[]; 
locsv =[];
locspd =[];
locsvd=[];

%detect peak section by section
peakdetectintT=240;%sec
peakdetectint=floor(peakdetectintT/TR/2);
pstepn=ceil(length(unmixed_decoupl)/peakdetectint);
for pseg=1:pstepn
    prange=((pseg-1)*peakdetectint+1:(pseg-1)*peakdetectint+...
        min(peakdetectint,length(unmixed_decoupl((pseg-1)*peakdetectint+1:end,Respp))));
locsp=[locsp,ampd(unmixed_decoupl(prange,Respp))+(pseg-1)*peakdetectint];
locsv=[locsv,ampd(-unmixed_decoupl(prange,Respp))+(pseg-1)*peakdetectint];
locspd=[locspd,ampd(diff(movmean(unmixed_decoupl(prange,Respp),...
    round((1/20)/(2*TR)))))+(pseg-1)*peakdetectint];
locsvd=[locsvd,ampd(diff(movmean(-unmixed_decoupl(prange,Respp),...
    round((1/20)/(2*TR)))))+(pseg-1)*peakdetectint];
%locspvd=ampd(-abs(diff(movmean(-unmixed_decoupl(:,Cardiacp),round((1/20)/(2*TR))))));
end
% locsp=ampd(unmixed_decoupl(:,Respp));
% locsv=ampd(-unmixed_decoupl(:,Respp));
% locspd=ampd(diff(movmean(unmixed_decoupl(:,Respp),round((1/20)/(2*TR)))));
% locsvd=ampd(diff(movmean(-unmixed_decoupl(:,Respp),round((1/20)/(2*TR)))));

Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';
recon=abs(dispimICA(reshape(reshape(U,Ny*Nx*Nz,[])*Phi_rt_temp(:,locsp),Ny,Nx,Nz,[])));
implay(abs(recon(:,:,:))/cw);

%
recon=abs(dispimICA(reshape(reshape(U,Ny*Nx*Nz,[])*Phi_rt_temp(:,locspd),Ny,Nx,Nz,[])));
implay(abs(recon(:,:,:))/cw);
%

% determine expiration for Respp
if abs(mean(unmixed_decoupl(locsp,Respp)-mean(unmixed_decoupl(:,Respp))))>abs(mean(unmixed_decoupl(locsv,Respp)-mean(unmixed_decoupl(:,Respp))))

expp=[expp,locsv+(cstep-1)*segtemp];
insp=[insp,locsp+(cstep-1)*segtemp];



else
%peak is systole
expp=[expp,locsp+(cstep-1)*segtemp];
insp=[insp,locsv+(cstep-1)*segtemp];


end
%
% determine diastolic for diffZ
[locpvd,idx]=sort([locsvd,locspd]);
locpvddif=diff(locpvd);
intvv=locpvddif(idx(1:end-1)<=length(locsvd));
intvp=locpvddif(idx(1:end-1)>length(locsvd));

if mean(intvp)>mean(intvv)
%valley is experiation

exppd=[exppd,locsvd+(cstep-1)*segtemp];
inspd=[inspd,locspd+(cstep-1)*segtemp];


else

exppd=[exppd,locspd+(cstep-1)*segtemp];
inspd=[inspd,locsvd+(cstep-1)*segtemp];

end


%%
end
%%
Triggeredsystag=systolep;
 x=[Triggeredsystag,Triggeredsystag-1];
 y=[ones(length(Triggeredsystag),1);ones(length(Triggeredsystag),1)*cbins];
xi=1:size(Phi_rt_small,2);
[xs,xsorder]=sort(x);
y=y(xsorder);
Hidx = interp1q(xs',y,xi');
Hidx=round(Hidx);
%clear masksig


%%


Hidx(isnan(Hidx))=2;
Hidx(Hidx==0)=2;


%% View

  HRinterv=1;

 HRintv=diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*1000;

figure;plot(HRintv(1:HRinterv:end))
ylabel('RR interval (ms)')

%% Ridx

Triggeredsystag=insp;
 x=[Triggeredsystag,Triggeredsystag-1];
 y=[ones(length(Triggeredsystag),1);ones(length(Triggeredsystag),1)*rbins];
xi=1:size(Phi_rt_small,2);
[xs,xsorder]=sort(x);
y=y(xsorder);
Ridx = interp1(xs',y,xi','linear');
Ridx=round(Ridx);
%clear masksig


%%


Ridx(isnan(Ridx))=2;
Ridx(Ridx==0)=2;


%% View

  Resinterv=1;

 Resintv=diff(find(diff(Ridx)==(1-rbins)))*params.lEchoSpacing;

figure;plot(Resintv(1:Resinterv:end))
ylabel('Res interval (s)')

%% expand to Phi_full

 %   Ridx = interp1(1:length(Phi_rt_small),Ridx,1:0.5:length(Phi_rt_small),'previous','extrap');
 %   Hidx = interp1(1:length(Phi_rt_small),Hidx,1:0.5:length(Phi_rt_small),'previous','extrap');
Ridx=Ridx';
Hidx=Hidx';
 %% other index
% comment to put out to binning_ICA
% Segidx=ones(size(Hidx));
% %wall clock
% switch ScanType
%   case 'T2IR' %replace T1T2 subspace with T1 subspace (we will learn the T2 subspace)
%     Nseg=Nseg/5;
%     params.lSegments=params.lSegments/5;
%     [curveU,curveS]=gen_curve_subspace(Nseg,params.lEchoSpacing,alpha_deg,alpha0_deg,(params.alTR_seconds-params.lEchoSpacing*Nseg)/2);
%     curvePhi = curveU(:,1:cL);
%     wall_clock=repmat(vec(repmat(1:5,[Nseg/2 1])),[numel(Hidx)/(Nseg/2*5) 1]);
%     wall_clock=wall_clock(1:numel(Hidx));
%   case 'Cine'
%     wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
%     wall_clock(end+1)=wall_clock(end);
%     wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
%     curvePhi=[];
%     Segidx(:)=1;
%   case {'IR','SR'}
%     wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
%     wall_clock(end+1)=wall_clock(end);
%     wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
%   case 'T2prep'
%     wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
%     wall_clock(end+1)=wall_clock(end);
%     wall_clock=ceil(wall_clock/5);
%     wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
%     %Randy added 2021
%     Segidx=mod((1:size(Phi_rt_small,2))-1,params.lSegments/SGblock*(SGblock-1))+1;%Original
%     [curveU,curveS]=gen_curve_subspace(Nseg,params.lEchoSpacing,alpha_deg,alpha0_deg,0);
%     curvePhi = curveU(1:SGblock:Nseg,1:cL);
%     %%Randy added for bloch end
%   case 'T2star'
% %     Ridx=repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco]);
% %     Ridx=Ridx(:);
% %     Segidx=vec(repmat(nav_indices(1:NnavsPerBlock),[params.NEco 1]));
% %     Segidx=repmat(Segidx,[numel(Ridx)/numel(Segidx), 1]);
% % %     wall_clock=ceil((1:numel(Ridx))/(2*params.NEco*NnavsPerBlock)).'; %every 2 SR periods (~1 sec)
% % %     wall_clock=ceil((1:numel(Ridx))/(params.NEco*NnavsPerBlock)).'; %every SR period (~0.5 sec)
%     Ridx = interp1(1:params.NEco:(params.NEco*length(Ridx)),Ridx,1:(params.NEco*length(Ridx)),'previous','extrap');
%     Hidx = interp1(1:params.NEco:(params.NEco*length(Hidx)),Hidx,1:(params.NEco*length(Hidx)),'previous','extrap');
%     
% %     Ridx = reshape(repmat(Ridx,1,8).',[],1);
%     Segidx = Ridx;
%     Segidx(:) = 1;
% %     Hidx = reshape(repmat(Hidx,1,8).',[],1);
%     wall_clock = vec(repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco])).';
% end
% 
