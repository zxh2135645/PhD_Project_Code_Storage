if isunix
  sep='/';
else
  sep='\';
end

dicoms=dir('DICOM/*.IMA');
AIFcount=1;
perfcount=1;
for j=1:numel(dicoms)
  header=dicominfo(strcat('DICOM',sep,dicoms(j).name));
  if strcmp(header.SeriesDescription,'Rest1_wip713srBIR4_tPAT2_8mmSL4_112x176_Res2p7x2p2_AIF')
    AIF(:,:,AIFcount)=dicomread(strcat('DICOM',sep,dicoms(j).name));
    AIFcount=AIFcount+1;
  elseif strcmp(header.SeriesDescription,'Rest1_wip713srBIR4_tPAT2_8mmSL4_112x176_Res2p7x2p2_MOCO_NORM')
    perfusion(:,:,perfcount)=dicomread(strcat('DICOM',sep,dicoms(j).name));
    interval(perfcount)=header.NominalInterval;
    perfcount=perfcount+1;
  elseif strcmp(header.SeriesDescription,'Rest1_wip713srBIR4_tPAT2_8mmSL4_112x176_Res2p7x2p2_MAP_SLOPE')
    slope=dicomread(strcat('DICOM',sep,dicoms(j).name));
  end
end


%%
h=figure;
imshow(abs(mean(AIF,3)),[])
title('Draw LV')
LV_roi = impoly(gca,'Closed',1);
LV_mask = createMask(LV_roi);
close(h)

perfmean=max(perfusion,[],3);
if false %auto-detect myocardium
  h=figure;
  imshow(perfmean,[])
  title('Draw Around Myo')
  myo_roi = impoly(gca,'Closed',1);
  myo_mask = createMask(myo_roi);
  close(h)

  f1=perfmean(myo_mask);
  centroid=regionprops(myo_mask,'centroid');
  [ygrid,xgrid]=ndgrid((1:size(myo_mask,1))-0.5,(1:size(myo_mask,2))-0.5);
  dfromc=sqrt((ygrid-centroid.Centroid(1)).^2+(xgrid-centroid.Centroid(2)).^2);
  f2=dfromc(myo_mask);
  
  f1=f1(:)/std(f1(:));
  f2=f2(:)/std(f2(:));
  [class,c]=kmeans([f1, f2],4); %air, RVblood, LVblood, myo
  
  [~,dist]=sort(c(:,2));
  [~,cont]=sort(c(:,1));
  if dist(2)==cont(2) %second-closest to center, second-darkest
    myo_class=dist(2);
  end
  
  classim=double(myo_mask);
  classim(myo_mask)=class;
  myo_mask=(classim==myo_class);
  temp=reshape(perfusion,[],size(perfusion,3));
  maxes=double(max(temp(myo_mask,:).'));
  outliers=(maxes>(median(maxes)+1.96*std(maxes)))&(maxes<(median(maxes)-1.96*std(maxes)));
  outlier_inds=find(myo_mask(:));
  outlier_inds=outlier_inds(outliers);
  myo_mask(outlier_inds)=0;
else %draw myocardium
  h=figure;
  imshow(perfmean,[])
  title('Draw Outer Myo Border')
  myo_roi = impoly(gca,'Closed',1);
  myo_mask = createMask(myo_roi);
  close(h)
  
  h=figure;
  imshow(perfmean,[])
  title('Draw Inner Myo Border')
  myo_roi = impoly(gca,'Closed',1);
  myoi_mask = createMask(myo_roi);
  close(h)
  
  myo_mask=myo_mask&~myoi_mask;
  
  
end

perfmean=repmat(double(perfmean)/double(max(perfmean(:))),[1 1 3]);
perfmean=perfmean*.98+.01;
tempb=perfmean(:,:,3);
tempb(myo_mask)=1;
perfmean(:,:,3)=tempb;
h=figure;
imshow(perfmean)
title('Draw Anterior Half')
ant_roi = impoly(gca,'Closed',1);
ant_mask = createMask(ant_roi);
close(h)

tempb=perfmean(:,:,3);
tempb(myo_mask & ~ant_mask)=0;
perfmean(:,:,3)=tempb;
h=figure;
imshow(perfmean)
title('Draw Septal Third')
sept_roi = impoly(gca,'Closed',1);
sept_mask = createMask(sept_roi);
close(h)

tempr=perfmean(:,:,1);
tempr(myo_mask & sept_mask)=1;
tempr(myo_mask & ~sept_mask)=0;
perfmean(:,:,1)=tempr;
h=figure;
imshow(perfmean)
title('Draw Lateral Third')
lat_roi = impoly(gca,'Closed',1);
lat_mask = createMask(lat_roi);
close(h)

tempg=perfmean(:,:,2);
tempg(myo_mask & lat_mask)=1;
tempg(myo_mask & ~lat_mask)=0;
perfmean(:,:,2)=tempg;
figure,imshow(perfmean)

seg_mask(:,:,1)=prod(perfmean==repmat(reshape([0 0 1],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,2)=prod(perfmean==repmat(reshape([1 0 1],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,3)=prod(perfmean==repmat(reshape([1 0 0],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,4)=prod(perfmean==repmat(reshape([0 0 0],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,5)=prod(perfmean==repmat(reshape([0 1 0],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,6)=prod(perfmean==repmat(reshape([0 1 1],1,1,3),[size(perfmean,1),size(perfmean,2),1]),3);
seg_mask(:,:,7)=myo_mask;
seg_mask=logical(seg_mask);

AIF=reshape(AIF,[],size(AIF,3));
perfusion=reshape(perfusion,[],size(perfusion,3));
LV=squeeze(mean(double(AIF(LV_mask,:)),1));
LV_sat=squeeze(mean(double(perfusion(imresize(LV_mask,sqrt(size(perfusion,1))*[1 1])>.5,:)),1));
for j=1:7
  myo(:,j)=squeeze(mean(double(perfusion(seg_mask(:,:,j),:)),1));
end


%remove 2 garbage ims
LV=LV(3:end);
LV_sat=LV_sat(3:end);
myo=myo(3:end,:);

%detect baseline length
for j=4:size(LV_sat,2) %at least length 3
  z=(LV_sat(j)-mean(LV_sat(1:j-1)))/std(LV_sat(1:j-1));
  if tcdf(z,j-2) > 0.95 %if 95% confident is wash-in
    base_end = j-1;
    break;
  end
end
LV=LV-mean(LV(1:base_end)); LV(1:base_end)=0;
LV_sat=LV_sat-mean(LV_sat(1:base_end)); LV_sat(1:base_end)=0;
%myo=myo-repmat(mean(myo(1:base_end,:)),[size(myo,1) 1]); myo(1:base_end,:)=0;

%calculate scaling factor
[template,idx]=sort(LV_sat(base_end+1:end));
input=LV(idx+base_end);
for j=ceil(numel(input)/2):numel(input)
  [~,s,v]=svd([input(1:j).',template(1:j).'],'econ');
  if s(end)/norm(diag(s)) > 0.025
    scalar=v(2)/v(1)
    break;
  end
end

%%
dt=median(interval)/1000/60; % in minutes

%truncate at local minimum after peak (stop before recirculation)
spline=[1 1 1]/3.';
spline=conv(conv(spline,spline,'full'),spline,'full');
LV_sm=conv(LV,spline,'valid'); %smooth AIF
[~,peak]=max(LV_sm);
LV_sm=LV_sm(peak:end);
len=3+find(sign(diff(LV_sm))>=0,1)+peak-1; %find zero-crossing of finite differences, then correct for convolution and wash-in removal

for seg=1:7
  l=LV*scalar;
  m=myo(:,seg).';
  
  maxm=max(m);
  l=l(1:len)/maxm;
  m=m(1:len)/maxm;
  % l=l/maxm; l(len+1:end)=0;
  % m=m/maxm;
  
  fermi=@(f)[zeros(1,f(4)), f(1)./(exp(((0:(len-1))-f(2))*f(3))+1)];
  resmax=inf;
  for j=len+(0:10)
    [f,res]=lsqnonlin(@(f)conv(l,fermi([f(1:3) j])*dt,'same')-m+mean(m(1:base_end+j-len)),[1 10 1],[0 0 0]);
    if res < resmax
      fbest=[f(1:3) j];
      resmax = res;
    end
  end
  figure(seg),subplot(3,1,1),plot(fermi(fbest))
  figure(seg),subplot(3,1,2),plot(l*maxm)
  figure(seg),subplot(3,1,3),plot(conv(l,fermi(fbest)*dt*maxm,'same'))
  hold all
  plot((m-mean(m(1:base_end+fbest(4)-len)))*maxm)
  drawnow
  
  % max(diff(m)/dt)/max(l)/1.05
  ambf(seg)=max(fermi(fbest)/1.05);
end

ambf

ambf_im=zeros(size(myo_mask));
for seg=1:6
  ambf_im(seg_mask(:,:,seg))=ambf(seg);
end
figure,imagesc(ambf_im,[0 2]),colormap(jet(256));
