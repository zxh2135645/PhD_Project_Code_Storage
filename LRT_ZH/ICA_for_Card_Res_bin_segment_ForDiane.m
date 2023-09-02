% Data binningwith ICA analysis
% Cedars Sinai 
% Randy Yang 2018 Oct
% Hsin-Jung.Yang@cshs.org

%% setup variables
clear masksig signalSeg Zseg f
ntemp=ceil(total_time/120);
Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';
segtemp=floor(size(Phi_rt_temp,2)/ntemp);
cstep=0;
% use the adjustvolume as the cardiac ROI
mask=ApplyAdjvolumeMask_LRT(twix_obj);
systolep=[];
diastolep=[];
systolepd=[];
diastolepd=[];
systoleprom=[];
diastoleprom=[];
systolewidth=[];
diastolwidth=[];
Diaframxy=0;

for ctemp=1:ntemp
    
    %% segmentwised realtime recon to get the ROI of heart
recon=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,(ctemp-1)*segtemp+1:(ctemp)*segtemp),Ny,Nx,[])));
cstep=cstep+1;
    tempmask=repmat(mask,[1,1,size(recon,3)]);
    tempmasksig=recon(tempmask);
    masksig((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3),:)=reshape(tempmasksig,sum(mask(:)),[])';



clear Phi_rt_temp tempmask tempmasksig
%% ICA analysis

%input=masksig;
input= masksig((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recon,3),:);
%
plotrange=1:size(input,1);
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

 figure(102); subplot(q+1,1,n); plot(unmixed(plotrange,n));  hold on
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
n = length(x);          % number of samples
f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range

df=fs/size(unmixed,1);% sample frequency (Hz)
%figure;
q=size(unmixed,2);
%%%%% prepmodulation frequency
fsm=1/(TR*2*params.lSegments);
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
figure(103);
set(gcf,'Name','ICA notch filtered results');
for n=1:q
  %figure;imagesc(temp)
  %subplot(q+1,1,n); plot(Nav_interval*plotrange,unmixed(plotrange,n));  
 figure(103); subplot(q+1,1,n); plot(unmixed_decoupl(plotrange,n));  hold on
 %hold on; subplot(q+1,1,n); plot(mod(plotrange,(params.lSegments/2))/(params.lSegments/2)*2*std(unmixed(plotrange,n))+mean(unmixed(plotrange,n)),'r'); 
end
%%
%%%signal decouple from the Prep pulse and trajectory modulation


%%%

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


Cardiacrange=(f>(45/60))&(f<(140/60)).*~harmonicsm;
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
%% Cardiac frequency filtering
fs=1/(params.lEchoSpacing*2);
df=fs/size(unmixed,1);
n = size(unmixed,1);          % number of samples
clear f
f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
harmonicsm=((abs(mod(f,fsm))<(df/2))|(abs(mod(f,fsm))>(fsm-(df/2)))); %harmocin frequency for prep modulation
winlp=2*floor((140/60)/df);
windowlp=zeros(size(unmixed,1),1);
windowlp(1:winlp)=1;
windowlp=circshift(windowlp,[-winlp/2, 0]);

winhp=2*floor((45/60)/df);
windowhp=zeros(size(unmixed,1),1);
windowhp(1:winhp)=1;
hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));
%hwindow=hwindow.*~harmonicsm';%demodulate the prepsignal


hf = designfilt('highpassiir','StopbandFrequency',30/60, ...
               'PassbandFrequency',45/60, ...
                'SampleRate',fs);
            Z= filtfilt(hf,unmixed_decoupl(:,Cardiacp));    

figure(160);hold on;subplot(ntemp,1,ctemp); plot(Z);
figure(180);hold on;subplot(ntemp,1,ctemp); plot(unmixed_decoupl(:,Cardiacp));
signalSeg(ctemp,:)=unmixed_decoupl(:,Cardiacp);
Zseg(ctemp,:)=Z;

   %% respiratory
   

%%
tm=1:length(Z);
%tm=tm*TR*2;
%%

 hrcutoff=140;
 clear triggerpoint
for n=1:length(hrcutoff);
%%

locsp=ampd(Z);
locsv=ampd(-Z);
locspd=ampd(diff(movmean(unmixed_decoupl(:,Cardiacp),round((1/20)/(2*TR)))));
locsvd=ampd(diff(movmean(-unmixed_decoupl(:,Cardiacp),round((1/20)/(2*TR)))));

Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';
recon=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,systolep),Ny,Nx,[])));
implay(abs(recon(:,:,:))/cw);

Peakmean=mean(mean(masksig(locsp,:)))
Vallymean=mean(mean(masksig(locsv,:)))

if Peakmean>Vallymean

systolep=[systolep,locsv+(cstep-1)*segtemp];
diastolep=[diastolep,locsp+(cstep-1)*segtemp];
systolepd=[systolepd,locsvd+(cstep-1)*segtemp];
diastolepd=[diastolepd,locspd+(cstep-1)*segtemp];


else
%peak is systole
systolep=[systolep,locsp+(cstep-1)*segtemp];
diastolep=[diastolep,locsv+(cstep-1)*segtemp];
systolepd=[systolepd,locspd+(cstep-1)*segtemp];
diastolepd=[diastolepd,locsvd+(cstep-1)*segtemp];

end

end
end
 x=[systolep,systolep-1];
 y=[ones(length(systolep),1);ones(length(systolep),1)*cbins];
xi=1:size(Phi_rt_init,2);
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

 HRintv=diff(find(diff(Hidx)==(1-cbins)))*params.lEchoSpacing*2000;

figure;plot(HRintv(1:HRinterv:end))
ylabel('RR interval (ms)')


