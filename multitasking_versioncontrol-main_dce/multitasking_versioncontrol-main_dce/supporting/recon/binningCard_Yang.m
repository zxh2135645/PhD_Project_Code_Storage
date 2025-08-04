function [dataArray,temporalBasis] = binningCard_Yang(params,reconOptions,dataArray,temporalBasis,spatialCoeff,selectROI)

% load variables
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L_init = size(Phi_rt_small_init,1);
ccL = min(10,size(curvePhi,2)); 
Segidx = mod(navIndices-1,linesPerShot*moduleLength)+1;

vec = @(x) x(:);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end

%% ROI
h = figure;imagesc(abs(fbpComposite(:,:,floor(Nz/2)+1,1)));axis equal tight;colormap('gray');title('Draw cardiac motion ROI')
roiResp = imellipse;
roiPosition = roiResp.getPosition();
roiPosition = round(roiPosition);
close(h);

Nzshift = floor(Nz/2)- floor(Nzorig/2);
Nxshift = roiPosition(1);
Nyshift = roiPosition(2);
ROINx = roiPosition(3);
ROINy = roiPosition(4);

windowfun = zeros(Ny,Nx);
windowfun(Nyshift+(1:ROINy),Nxshift+(1:ROINx)) = 1;
roi_weighting = zeros(Ny,Nx,Nz);
roi_weighting(:,:,Nzshift+(1:Nzorig)) = repmat(windowfun,[1 1 Nzorig]);

if MBfactor == 1
    roi_weighting = ifftshift(roi_weighting,3);
end
if isCartesian
    roi_weighting = ifftshift(roi_weighting,1);
end
%% ICA
Wti = reshape(U_init,[],L_init);
Wti = bsxfun(@times,Wti,roi_weighting(:));
Wti = Wti'*Wti;
Wti = sqrtm(Wti); %actually inverse of Wt
U = vec(reshape(U_init,[],L_init)/Wti);

Phi_rt       = Wti*Phi_rt_init;
Phi_rt_full  = Wti*Phi_rt_full_init;
Phi_rt_small = Wti*Phi_rt_small_init;

temporalBasis.Phi_rt = Phi_rt;
temporalBasis.Phi_rt_full  = Phi_rt_full;
temporalBasis.Phi_rt_small = Phi_rt_small;

ntemp=5;
segtemp=floor(size(Phi_rt_init,2)/ntemp);
cstep=0;
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :);

navdata_temp=reshape(reshape(navData,[],Ncoils)*mixer(:,1:(Ncoils/3)),size(navData,1),[]);
clear Phi_new signalSeg Zseg recon_mask

for ctemp=1:ntemp
    recontemp=abs(reshape(reshape(U_init,Ny*Nx,[])*Phi_rt_init(:,(ctemp-1)*segtemp+1:(ctemp)*segtemp),Ny,Nx,[]));
    recontemp=double(recontemp);
    cw=prctile(abs(recontemp(:)),99);
    implay(dispim(double(recontemp))./cw);
    cstep=cstep+1;
    tempmask=repmat(roi_weighting,[1,1,size(recontemp,3)]);
    recon_mask_temp=recontemp(logical(tempmask));
    recon_mask((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recontemp,3),:)=reshape(recon_mask_temp,sum(roi_weighting(:)),[])';
    mask2=repmat(roi_weighting,[1,1,size(reshape(U_init,Ny,Nx,[]),3)]);
    U_temp=reshape(U_init,Ny,Nx,[]); U_mask=U_temp(logical(mask2));
    recon_mask_temp2=recon_mask((cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recontemp,3),:);
    recon_mask_temp2=recon_mask_temp2';
    proj_nav=pinv(reshape(U_mask,[],size(reshape(U_init,Ny*Nx,[]),2)))*recon_mask_temp2;
    Phi_new(:,(cstep-1)*segtemp+1:(cstep-1)*segtemp+size(recontemp,3))=proj_nav;
    proj_nav=[real(proj_nav);imag(proj_nav)];

    clear tempmask recon_mask_temp cscore rscore

    absnav=abs(navdata_temp)';
    figure; imagesc(absnav(:,1:segtemp)); colorbar; hold on;

    q=8; % number of features
    Mdl1=rica(proj_nav',q,'NonGaussianityIndicator',ones(q,1),'IterationLimit',5000,'Lambda',1.0);
    nav_ica=transform(Mdl1,proj_nav');
    plotrange=1:size(proj_nav',1);
    figure; hold on; set(gcf,'Name',['ICA results - ', num2str(ctemp)]);
    for i=1:q
        subplot(q,1,i);
        plot(nav_ica(plotrange,i));
    end
    %% Notch Filter
    TR=Params.lEchoSpacing;
    fs = 1/2/TR;
    n = size(nav_ica,1);          % number of samples
    f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
    f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
    df=fs/size(nav_ica,1);% sample frequency (Hz)
    q=size(nav_ica,2);
    fsm=1/(TR*2*Params.lSegments); % prepmodulation frequency
    temp=nav_ica;
    for fn=1:floor(max(f)/fsm)
        d(fn) = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',fsm*fn-df/0.5,'HalfPowerFrequency2',fsm*fn+df/0.5, ...
            'DesignMethod','butter','SampleRate',fs);
        for i=1:q
            temp(:,i)= filtfilt(d(fn),temp(:,i));
        end
    end
    filt_nav=temp;
    figure; set(gcf,'Name',['ICA notch filtered results - ', num2str(ctemp)]); hold on
    for i=1:q
        subplot(q,1,i); plot(filt_nav(plotrange,i));
    end
    clear temp
    unmixed_decoupl=filt_nav;
    %% Signal decouple
    for p=1:q
        x=abs(unmixed_decoupl(:,p));

        y = fft(x);
        n = length(x);          % number of samples
        % f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
        % f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
        harmonicsm=(abs(mod(f,fsm))<df)|(abs(mod(f,fsm))>(fsm-df)); %harmocin frequency for prep modulation
        powr(p,:) = abs(y).^2/n;    % power of the DFT
        frange=(f>0.001)&(f<5);
        figure(100);subplot(q,1,p); plot(f(frange),powr(p,frange)); hold on%/sum(power(:,p)));

        Cardiacrange=(f>(35/60))&(f<(120/60)).*~harmonicsm;
        respiratoryrange=(f>(4/60))&(f<(30/60)).*~harmonicsm;
        AcceptRange=(f>(4/60))&(f<(300/60)).*~harmonicsm;
        prepmodscore(p)=sum(powr(p,frange).*harmonicsm(frange))/sum(powr(p,AcceptRange));
        cscore(p)=sum(powr(p,frange).*Cardiacrange(frange))/sum(powr(p,frange).*AcceptRange(frange));
        rscore(p)=sum(powr(p,frange).*respiratoryrange(frange))/sum(powr(p,frange).*AcceptRange(frange));
    end

    Cardiacp=find(cscore==max(cscore)); disp(['Cardiacp = ',num2str(Cardiacp)]);
    Respp=find(rscore==max(rscore)); disp(['Respp = ',num2str(Respp)]);
    Carid(ctemp)=Cardiacp;
    Respid(ctemp)=Respp;

    figure(101);
    subplot(ntemp,1,ctemp); plot(f(frange),powr(Cardiacp,frange)); set(gcf,'Name','Cardiac phase spectrum');
    figure(102);
    subplot(ntemp,1,ctemp); plot(f(frange),powr(Cardiacp,frange).*Cardiacrange(frange)); set(gcf,'Name','Cardiac phase spectrum .* cardiac range');
    figure(103);
    subplot(ntemp,1,ctemp); plot(f(frange),powr(Respp,frange)); set(gcf,'Name','Respiratory phase spectrum');
    figure(104);
    subplot(ntemp,1,ctemp); plot(f(frange),powr(Respp,frange).*respiratoryrange(frange)); set(gcf,'Name','Respiratory phase spectrum .* respiratory range');
    %% Cardiac filtering
    fs=1/(Params.lEchoSpacing*2);
    df=fs/size(nav_ica,1);
    n = size(nav_ica,1);          % number of samples
    clear f
    f(1:n/2) = (0:n/2-1)*(fs/n);     % frequency range
    f(n/2+1:n) = (-(n)/2:1:-1)*(fs/n);     % frequency range
    harmonicsm=((abs(mod(f,fsm))<(df/2))|(abs(mod(f,fsm))>(fsm-(df/2)))); %harmocin frequency for prep modulation
    winlp=2*floor((140/60)/df);
    windowlp=zeros(size(nav_ica,1),1);
    windowlp(1:winlp)=1;
    windowlp=circshift(windowlp,[-winlp/2, 0]);

    winhp=2*floor((45/60)/df);
    windowhp=zeros(size(nav_ica,1),1);
    windowhp(1:winhp)=1;
    hwindow=windowlp.*(1-circshift(windowhp,[-winhp/2, 0]));

    hf = designfilt('highpassiir','StopbandFrequency',40/60, ...
        'PassbandFrequency',120/60, ...
        'SampleRate',fs);


    Z= filtfilt(hf,unmixed_decoupl(:,Cardiacp));

    figure(105);hold on; subplot(ntemp,1,ctemp); plot(unmixed_decoupl(:,Cardiacp)); set(gcf,'Name','Cardiac signal before filtering');
    figure(106);hold on; subplot(ntemp,1,ctemp); plot(Z); set(gcf,'Name','Cardiac signal after filtering');
    signalSeg(ctemp,:)=unmixed_decoupl(:,Cardiacp);
    Zseg(ctemp,:)=Z;
end
%% RR interval & save
systolep=[];
diastolep=[];
Zsegint=Zseg;
Zseg=Zseg';
Zseg=Zseg(:);
Zsegtemp=buffer(Zseg,segtemp+500,500,'nodelay');
signalSegint=signalSeg;
signalSeg=signalSeg';
for i=1:size(Zsegtemp,2)
    Z=Zsegtemp(:,i);
    sig=signalSeg(:,i);
    locsp=ampd(Z);
    locsv=ampd(-Z);
%     [~,locsp]=findpeaks(Z,'MinPeakDistance',45,'MinPeakHeight',0.0015);
%     locsp=locsp';
%     [~,locsv]=findpeaks(-Z,'MinPeakDistance',45,'MinPeakHeight',0.0015);
%     locsv=locsv';
    %locsv=ampd(-Z);
    locspd=ampd(diff(movmean(sig,round((1/20)/(2*TR)))));
    locsvd=ampd(diff(movmean(-sig,round((1/20)/(2*TR)))));
    
    tempv=locsv+(i-1)*segtemp;
    tempp=locsp+(i-1)*segtemp;
    Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';
    recon=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,tempv),Ny,Nx,[])));
    recon1=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,tempp),Ny,Nx,[])));
    %implay(recon./max(recon(:)));
    %mask=ApplyAdjvolumeMask_LRT(twix_obj);
    %[newmask]=gaussmask(mask);
    newmask=dispim(roi_weighting);
    temp_recon_v=recon.*repmat(newmask,[1,1,size(recon,3)]);
    temp_recon_p=recon1.*repmat(newmask,[1,1,size(recon1,3)]);
    sbloodp=0; sbloodv=0;
    for k=1:size(temp_recon_p,3)
        non0p=nonzeros(temp_recon_p(:,:,k));
        threp=max(non0p)*0.6;
        bloodp=non0p(non0p>threp);
        sbloodp=size(bloodp,1)+sbloodp;
    end
    sbloodp=sbloodp/size(temp_recon_p,3);
    for j=1:size(temp_recon_v,3)
        non0v=nonzeros(temp_recon_v(:,:,j));
        threv=max(non0v)*0.6;
        bloodv=non0v(non0v>threv);
        sbloodv=size(bloodv,1)+sbloodv;
    end
    sbloodv=sbloodv/size(temp_recon_v,3);
    
    locsv=locsv(locsv<segtemp+1);
    locsp=locsp(locsp<segtemp+1);
    
    if sbloodp>sbloodv
        % valley is systole
        % RR=diff(locsv);
        systolep=[systolep,locsv+(i-1)*segtemp];
        diastolep=[diastolep,locsp+(i-1)*segtemp];
        
    else
        %peak is systole mean_v>mean_p
        systolep=[systolep,locsp+(i-1)*segtemp];
        diastolep=[diastolep,locsv+(i-1)*segtemp];
    end
    peaknum(i)=sbloodp;
    valleynum(i)=sbloodv;
end
systolep=unique(systolep);
systolep=systolep';
diastolep=unique(diastolep);
diastolep=diastolep';

disp(['Peaknum = ', num2str(peaknum)]);
disp(['Vallynum = ', num2str(valleynum)]);

recon=abs(dispim(reshape(reshape(U_init,Ny*Nx,[])*Phi_rt_init(:,systolep),Ny,Nx,[])));
implay(abs(recon(:,:,:))./cw);
recon=abs(dispim(reshape(reshape(U_init,Ny*Nx,[])*Phi_rt_init(:,diastolep),Ny,Nx,[])));
implay(abs(recon(:,:,:))/cw);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(4,1,1); plot(Zseg(:),'-','Marker','square','MarkerIndices',systolep,'MarkerEdgeColor','r'); xlim([0 inf]); title('systolep');
subplot(4,1,2); plot(Zseg(:),'-','Marker','square','MarkerIndices',diastolep,'MarkerEdgeColor','r'); xlim([0 inf]); title('diastolep');

% Gating
cbins=ReconOptions.cbins;
x=[systolep;systolep-1];
y=[ones(length(systolep),1);ones(length(systolep),1)*cbins];
xi=1:size(Phi_rt_init,2);
%xi=1:size(Zsegtemp,1);
[xs,xsorder]=sort(x);
y=y(xsorder);
Hidx = interp1q(xs,y,xi');
Hidx=round(Hidx);
Hidx(isnan(Hidx))=2;
Hidx(Hidx==0)=2;
figure; plot(Hidx); hold on; yyaxis right; plot(Zseg(:));

HRinterv_temp=1;
HRintv=diff(find(diff(Hidx)==(1-cbins)))*Params.lEchoSpacing*2000;
%HRintv=movmean(HRintv,5);
figure;plot(HRintv(1:HRinterv_temp:end));
ylabel('RR interval (ms)');
%% Output

