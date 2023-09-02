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
load(sprintf('%s%sAC_recon.mat',found{folder},sep),'U','Phi','L','Gr',...
  'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec')

%% Choose recovery images
Phi=reshape(Phi,[L sizes(2:end)]);

temp=Gr\reshape(Phi(:,21,:,1,1),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp1=abs(dispim(temp))/max(abs(vec(dispim(temp))));

temp=Gr\reshape(Phi(:,71,:,1,1),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp2=abs(dispim(temp))/max(abs(vec(dispim(temp))));

temp=Gr\reshape(Phi(:,111,:,1,1),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp3=abs(dispim(temp))/max(abs(vec(dispim(temp))));

implay(cat(2,temp1,temp2,temp3));
clear temp*;

cphase=input('Which phase [1]? ');
if isempty(cphase)
  cphase=1;
end

recon = Gr\reshape(Phi(:,:,cphase,1,:),L,[]);
recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
recon=reshape(recon,N-ovs,N-ovs,Nseg/2,[]);

% final_val=false;
% while ~final_val
%   implay(abs(recon(:,:,:,1))/max(abs(recon(:))))
%   imstart=input('Choose starting image [21]? ');
%   if isempty(imstart)
    imstart=21;
%   end
%   implay(abs(recon(:,:,imstart,:))/max(abs(recon(:))))
%   final_val=input(sprintf('Confirm starting image %d (y/n) [y]: ',imstart),'s');
%   if isempty(final_val)
%     final_val=true;
%   else
%     final_val=strcmpi(final_val,'y') || strcmpi(final_val,'yes');
%   end
% end

recon=recon(:,:,imstart:end,:);
Nseg=2*size(recon,3);
recon=recon(:,:,:);

%% Fit
sliceprof=false; %true;
alpha0_deg=180;

planningim=abs(recon(:,:,1475-imstart*5));
%planningim=ind2rgb(uint8(255./fits(4).fits(:,:,2)/2.3),jet(256));

h=figure;
imshow(planningim,[])
title('Draw Outer Myo Border')
myo_roi = impoly(gca,'Closed',1);
myo_mask = createMask(myo_roi);
close(h)

h=figure;
imshow(planningim,[])
title('Draw Inner Myo Border')
myo_roi = impoly(gca,'Closed',1);
myoi_mask = createMask(myo_roi);
close(h)

im_mask=myo_mask&~myoi_mask;

planningim=repmat(double(planningim)/double(max(planningim(:))),[1 1 3/size(planningim,3)]);
planningim=planningim*.98+.01;
tempb=planningim(:,:,3);
tempb(im_mask)=1;
planningim(:,:,3)=tempb;
h=figure;
imshow(planningim)
title('Draw Septal Region')
sept_roi = impoly(gca,'Closed',1);
sept_mask = createMask(sept_roi) & im_mask;
close(h)

%%
alpha = params.adFlipAngleDegree*pi/180;

e1 = @(R1)exp(-params.lEchoSpacing*R1);
TEs=[12 20 30 40 50]*1e-3;
e2 = @(R2)exp(-TEs*R2);
Mss = @(e1,alpha)(1-e1) / (1-cos(alpha)*e1);
n = params.lSegments-Nseg+1:2:params.lSegments;

if true %if new
  %Calculate pct. of steady state reached
  syms R1 R2 B Eff_last alpha
  Eff_prev=Eff_last;
  for j=1:numel(TEs)
    Eff_sym(j)=(1 + ((exp(-params.lEchoSpacing*R1)*cos(alpha))^params.lSegments).'*(B*Eff_prev*exp(-TEs(j)*R2)-1));
    Eff_prev=Eff_sym(j);
  end
  Eff_last=solve(Eff_last==Eff_sym(end),Eff_last);
  Eff_sym(end)=Eff_last;
  evalstring='Eff=@(R1,R2,alpha,B)cat(2';
  for j=1:numel(TEs)
    evalstring=strcat(evalstring,sprintf(',%s',char(subs(Eff_sym(j)))));
  end
  evalstring=strcat(evalstring,');');
  eval(evalstring);
  
  Sint = @(A,e1,e2,alpha,B,Eff)vec(A * Mss(e1,alpha) * (1 + ((e1*cos(alpha)).^(n-1)).'*(B*e2.*Eff-1)) * sin(alpha)).';
  S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B,Eff(R1,R2,alpha,B));
else
  Sint = @(A,e1,e2,alpha,B)vec(A * Mss(e1,alpha) * (1 + ((e1*cos(alpha)).^(n-1)).'*(B*e2-1)) * sin(alpha)).';
  S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B);
end
ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)

minB = -1;
maxB = 1;

if ~exist('im_mask','var')
  h=figure;
  imshow(abs(recon(:,:,end)),[]);
  title('Select Fitting ROI')
  im_mask=createMask(imrect);
  close(h)
  drawnow
  
  im_mask=im_mask.*abs(recon(:,:,end));
  [im_mask,centroids]=kmeans(im_mask(:),2);
  im_mask=reshape(im_mask==1+(centroids(2)>centroids(1)),size(recon,1),size(recon,2));
end
figure,imshow(im_mask),drawnow;

[rows, cols]=ind2sub(size(im_mask),find(im_mask));
recontemp=reshape(recon,size(recon,1)*size(recon,2),[]);
recontemp=recontemp(im_mask(:),:);

opts=[];
opts.MaxFunEvals = 1000;

for type=4
  alpha = params.adFlipAngleDegree*pi/180;
  switch type
    case 1 %fixed alpha
      alpha=2.5*pi/180
    case 2 %fixed alpha from septum
      curve=reshape(recon,[],size(recon,3));
      [~,~,curve]=svd(curve(sept_mask,:),'econ');
      curve=curve(:,1).';
      normcurve = curve(end);
      curve = curve/normcurve;
      Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
      tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve),...
        [1/1.1,20,alpha/2,median([minB, real(curve(1)/curve(end)), maxB])],...
        [1/3,1,0,minB], [1/.1,1/.01,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
      alpha=tempfit(3);
    case 3 %fixed alpha from whole myocardium
      curve=reshape(recon,[],size(recon,3));
      [~,~,curve]=svd(curve(im_mask,:),'econ');
      curve=curve(:,1).';
      normcurve = curve(end);
      curve = curve/normcurve;
      Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
      tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve),...
        [1/1.1,20,alpha/2,median([minB, real(curve(1)/curve(end)), maxB])],...
        [1/3,1,0,minB], [1/.1,1/.01,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
      alpha=tempfit(3);
  end
  
  fitmat = zeros(numel(rows),6);
  switch type
    case {1, 2, 3}
      parfor j=1:numel(rows)
        curve = double(recontemp(j,:));
        normcurve = curve(end);
        curve = curve/normcurve;
        Avp = @(R1,R2,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
        tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha,x(3))-curve),...
          [1/1.1,20,median([minB, real(curve(1)/curve(end)), maxB])],...
          [1/3,1/.1,minB], [1/.1,1/.01,maxB], opts);
        tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit(1), tempfit(2), alpha, tempfit(3)];
        res=norm(curve*normcurve-S(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5))); %residual
        tempfit(6)=res;
        fitmat(j,:) = tempfit;
      end
    case 4 %fit alpha
      parfor j=1:numel(rows)
        curve = double(recontemp(j,:));
        normcurve = curve(end);
        curve = curve/normcurve;
        Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
        tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve),...
          [1/1.1,20,alpha/2,median([minB, real(curve(1)/curve(end)), maxB])],...
          [1/3,1/.1,alpha/4,minB], [1/.1,1/.01,alpha,maxB], opts);
        tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4))*normcurve, tempfit];
        res=norm(curve*normcurve-S(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5))); %residual
        tempfit(6)=res;
        fitmat(j,:) = tempfit;
      end
  end
  
  tempfit=fitmat;
  fits(type).fits=zeros(size(recon,1)*size(recon,2),6);
  fits(type).fits(im_mask(:),:)=tempfit;
  fits(type).fits=reshape(fits(type).fits,size(recon,1),size(recon,2),6);
  
  figure,imagesc(1./fits(type).fits(:,:,2),[0 2.3]),colormap(jet(256))
  figure,imagesc(1./fits(type).fits(:,:,3),[0 0.2]),colormap(jet(256))
  %   figure,imagesc(fits(type).fits(:,:,5)),colormap(jet(256))
  drawnow
  
  fitmat=reshape(fits(type).fits,[],size(fits(type).fits,3));
  myoT1(type,1)=mean(1./fitmat(im_mask,2));
  septT1(type,1)=mean(1./fitmat(sept_mask,2));
  myoT2(type,1)=mean(1./fitmat(im_mask,3));
  septT2(type,1)=mean(1./fitmat(sept_mask,3));
  
  myoT1_std(type,1)=std(1./fitmat(im_mask,2));
  septT1_std(type,1)=std(1./fitmat(sept_mask,2));
  myoT2_std(type,1)=std(1./fitmat(im_mask,3));
  septT2_std(type,1)=std(1./fitmat(sept_mask,3));
  
  temp=reshape(recon,size(recon,1)*size(recon,2),[]);
  w = 1-fitmat(im_mask,6)./sqrt(sum(abs(temp(im_mask,:)).^2,2)); %max(-fitmat(im_mask,5),0);
  myoT1(type,2)=sum(w./fitmat(im_mask,2))/sum(w);
  w = 1-fitmat(sept_mask,6)./sqrt(sum(abs(temp(sept_mask,:)).^2,2)); %max(-fitmat(sept_mask,5),0);
  septT1(type,2)=sum(w./fitmat(sept_mask,2))/sum(w);
  w = 1-fitmat(im_mask,6)./sqrt(sum(abs(temp(im_mask,:)).^2,2)); %max(-fitmat(im_mask,5),0);
  myoT2(type,2)=sum(w./fitmat(im_mask,3))/sum(w);
  w = 1-fitmat(sept_mask,6)./sqrt(sum(abs(temp(sept_mask,:)).^2,2)); %max(-fitmat(sept_mask,5),0);
  septT2(type,2)=sum(w./fitmat(sept_mask,3))/sum(w);
end
myoT1
septT1
myoT2
septT2

save(strcat(found{folder},sep,'fit_',datestr(now,30)),'fits','recon','im_mask','cphase','imstart','myo*','sept*')
