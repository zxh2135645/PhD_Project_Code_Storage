%% Select IR-FLASH recon
if ispc
  sep = '\';
else
  sep = '/';
end

folders=dir;
found=[];
k=0;
for j=1:numel(folders);
  if folders(j).isdir
    file=dir(sprintf('%s%sAC_recon.mat',folders(j).name,sep));
    if numel(file) > 0;
      found{k+1}=folders(j).name;
      k=k+1;
    end
  end
end

folder=listdlg('ListString',found,'SelectionMode','single','ListSize',[500 150]);
isfirstbolus = questdlg('Is this the first bolus?', ...
  'First bolus?', ...
  'Yes','No', 'Yes');
isfirstbolus = strcmp(isfirstbolus,'Yes');

load(sprintf('%s%sAC_recon.mat',found{folder},sep),'U','Phi','L','Gr',...
  'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec','wall_clock','Ridx')

%% Choose recovery images
Phi=reshape(Phi,[L sizes(2:end)]);

temp=Gr\reshape(Phi(:,end,:,1,:),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp=abs(dispim(temp))/max(abs(vec(dispim(temp))));
temp=mean(reshape(temp,N-ovs,N-ovs,sizes(3),[]),4);

implay(temp);

cphase=input('Which phase [1]? ');
if isempty(cphase)
  cphase=1;
end

recon = Gr\reshape(Phi(:,:,cphase,1,:),L,[]);
recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
recon=reshape(recon,N-ovs,N-ovs,Nseg/2,[]);

implay(abs(recon(:,:,end,:))/max(abs(recon(:))))

domoco=input('Do MOCO? (y/n) [y]: ','s');
if isempty(domoco)
  domoco=true;
else
  domoco=strcmpi(domoco,'y') || strcmpi(domoco,'yes');
end
if domoco
  mocostuff;
  usemoco=input('Use MOCO? (y/n) [y]: ','s');
  if isempty(usemoco)
    usemoco=true;
  else
    usemoco=strcmpi(usemoco,'y') || strcmpi(usemoco,'yes');
  end
  if usemoco
    recon=moco;
  end
end

%%
perfmean=mean(abs(recon(:,:,end,:)),4);
h=figure;
imshow(perfmean,[])
title('Draw LV')
LV_roi = impoly(gca,'Closed',1);
LV_mask = createMask(LV_roi);
close(h)

perfmean=max(max(abs(recon),[],3),[],4);
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

%%
Nt=sizes(5);
recon=reshape(recon,[],Nseg/2,Nt);

LV=squeeze(mean(recon(LV_mask,:,:)./repmat(sign(recon(LV_mask,end,end)),[1 Nseg/2 Nt])));
myo=[];
for seg=1:7
  myo(:,:,seg)=squeeze(mean(recon(seg_mask(:,:,seg),:,:)./repmat(sign(recon(seg_mask(:,:,seg),end,end)),[1 Nseg/2 Nt])));
end

alpha = params.adFlipAngleDegree*pi/180;

e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) ./ (1-cos(alpha)*e);

n = 1:2:Nseg;
Sint = @(A,e,alpha,B)A * repmat(Mss(e,alpha),[1 Nseg/2]) .* (1 + repmat(B-1,[1 Nseg/2]).*(repmat(e*cos(alpha),[1 Nseg/2]).^repmat(n-1,[Nt 1]))) * sin(alpha);
S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);

ppinv=@(x,y)(y(:)'*x(:))/norm(x(:))^2; %fast left-sided pseudoinverse function

minB = -.5;
maxB = 1;


opts=[];
opts.MaxFunEvals = 10000;

%%
if isfirstbolus
  prior=1.15;
else
  prior=0.6;
end
seg=7;
curve=myo(:,:,seg).';

normcurve = max(curve(end,:));

Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve/normcurve); %parameterize solution to A as function of R1,alpha,

x0=real(curve(:,ceil(end/2)).'/normcurve); x0=x0/prior/x0(1);
x0=[x0, alpha, 0];
xmin=[ones(1,Nt)/3, 0, minB];
xmax=[ones(1,Nt)*100, alpha, maxB];
cost=@(x)abs([vec(x(1:10)-1/prior); vec(S(Avp(vec(x(1:Nt)),x(Nt+1),x(Nt+2)*ones(Nt,1)),vec(x(1:Nt)),x(Nt+1),x(Nt+2)*ones(Nt,1))-curve/normcurve)]);
x=lsqnonlin(cost, x0, xmin, xmax, opts);

myofit(seg).A=Avp(vec(x(1:Nt)),x(Nt+1),x(Nt+2)*ones(Nt,1))*normcurve;
myofit(seg).R1=x(1:Nt);
myofit(seg).alpha=x(Nt+1);
myofit(seg).B=x(Nt+2);


for seg=1:6
  curve=myo(:,:,seg).';
  
  normcurve = max(curve(end,:));
  
  Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve/normcurve); %parameterize solution to A as function of R1,alpha,
  
  x0=real(curve(:,ceil(end/2)).'/normcurve); x0=x0/prior/x0(1);
  xmin=ones(1,Nt)/3;
  xmax=ones(1,Nt)*100;
  cost=@(x)abs([vec(x(1:10)-1/prior); vec(S(Avp(vec(x(1:Nt)),myofit(7).alpha,myofit(7).B*ones(Nt,1)),vec(x(1:Nt)),myofit(7).alpha,myofit(7).B*ones(Nt,1))-curve/normcurve)]);
  x=lsqnonlin(cost, x0, xmin, xmax, opts);
  
  myofit(seg).A=Avp(vec(x(1:Nt)),myofit(7).alpha,myofit(7).B*ones(Nt,1))*normcurve;
  myofit(seg).R1=x(1:Nt);
  myofit(seg).alpha=myofit(end).alpha;
  myofit(seg).B=myofit(7).B;
end

%%
if isfirstbolus
  prior=1.9;
else
  prior=0.6;
end

curve=LV.';

normcurve = max(curve(end,:));

Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve/normcurve); %parameterize solution to A as function of R1,alpha,B

x0=real(curve(:,ceil(end/2)).'/normcurve); x0=x0/prior/x0(1);
x0=[x0, 0];
xmin=[ones(1,Nt)/3, minB];
xmax=[ones(1,Nt)*100, maxB];
cost=@(x)abs([vec(x(1:5)-1/prior)/5; vec(S(Avp(vec(x(1:Nt)),myofit(7).alpha,x(Nt+1)*ones(Nt,1)),vec(x(1:Nt)),myofit(7).alpha,x(Nt+1)*ones(Nt,1))-curve/normcurve)]);
x=lsqnonlin(cost, x0, xmin, xmax, opts);

LVfit.A=Avp(vec(x(1:Nt)),myofit(7).alpha,x(Nt+1)*ones(Nt,1))*normcurve;
LVfit.R1=x(1:Nt);
LVfit.alpha=myofit(7).alpha;
LVfit.B=x(Nt+1);

%%
%detect baseline length
for j=6:size(LVfit.R1,2) %at least length 5
  z=(LVfit.R1(j)-mean(LVfit.R1(1:j-1)))/std(LVfit.R1(1:j-1));
  if tcdf(z,j-2) > 0.95 %if 95% confident is wash-in
    base_end = j-1;
    break;
  end
end

%truncate at local minimum after peak (stop before recirculation)
spline=[1 1 1]/3.';
spline=conv(conv(spline,spline,'full'),spline,'full');
LV_sm=conv(LVfit.R1,spline,'valid'); %smooth AIF
[~,peak]=max(LV_sm);
LV_sm=LV_sm(peak:end);
len=3+find(sign(diff(LV_sm))>=0,1)+peak-1; %find zero-crossing of finite differences, then correct for convolution and wash-in removal
if isempty(len)
  len=numel(LVfit.R1);
end


%%
dt=find(wall_clock==2,1)*params.alTR_seconds/(Nseg/2)/60; % in minutes
for seg=1:7
  l=LVfit.R1-mean(LVfit.R1(1:base_end)); l(1:base_end)=0;
  m=myofit(seg).R1;
  
  maxm=max(m);
  l=l(1:len)/maxm;
  m=m(1:len)/maxm;
  % l=l/maxm; l(len+1:end)=0;
  % m=m/maxm;
  
  fermi=@(f)[zeros(1,f(4)), f(1)./(exp(((0:(len-1))-f(2))*f(3))+1)];
  resmax=inf;
  for j=len+(0:10)
    [f,res]=lsqnonlin(@(f)conv(l,fermi([f(1:3) j])*dt,'same')-m+mean(m(1:base_end+j-len)),[1 10 1],[0 0 1]);
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
  plot((base_end+fbest(4)-len)*[1 1],ylim,'--')
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

%%
%save(sprintf('%s%sperf_analysis.mat',found{folder},sep))