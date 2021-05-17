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
  'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec','cbins')

imstart=21;
Nsegnew=2*(Nseg/2-imstart+1);

%% Setup functions
sliceprof=false; %true;
alpha0_deg=180;
alpha = params.adFlipAngleDegree*pi/180;

e1 = @(R1)exp(-params.lEchoSpacing*R1);
TEs=[12 20 30 40 50]*1e-3;
e2 = @(R2)exp(-TEs*R2);
Mss = @(e1,alpha)(1-e1) / (1-cos(alpha)*e1);
n = params.lSegments-Nsegnew+1:2:params.lSegments;

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

%% Fit
Phi=reshape(Phi,[L sizes(2:end)]);

for cphase=1:cbins
  
  recon = Gr\reshape(Phi(:,:,cphase,1,:),L,[]);
  recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
  recon=reshape(recon,N-ovs,N-ovs,Nseg/2,[]);
  
  recon=recon(:,:,imstart:end,:);
  recon=recon(:,:,:);
  
  if ~exist('im_mask','var')
    h=figure;
    imshow(abs(recon(:,:,end)),[]);
    title('Select Fitting ROI')
    im_mask=createMask(imrect);
    close(h)
    drawnow
    
    %   im_mask=im_mask.*abs(recon(:,:,end));
    %   [im_mask,centroids]=kmeans(im_mask(:),2);
    %   im_mask=reshape(im_mask==1+(centroids(2)>centroids(1)),size(recon,1),size(recon,2));
    
    figure,imshow(im_mask),drawnow;
  end
  
  
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
    fits(cphase).fits=zeros(size(recon,1)*size(recon,2),6);
    fits(cphase).fits(im_mask(:),:)=tempfit;
    fits(cphase).fits=reshape(fits(cphase).fits,size(recon,1),size(recon,2),6);
    
    figure(21),imagesc(wmedfilt2(1./fits(cphase).fits(:,:,2)),[0 2.3]),colormap(jet(256))
    figure(22),imagesc(wmedfilt2(1./fits(cphase).fits(:,:,3)),[0 0.2]),colormap(jet(256))
    %   figure,imagesc(fits(type).fits(:,:,5)),colormap(jet(256))
    drawnow
    
    T1cine(:,:,cphase)=1./fits(cphase).fits(:,:,2);
    T2cine(:,:,cphase)=1./fits(cphase).fits(:,:,3);
  end
end
  
save(strcat(found{folder},sep,'2017cinefit_',datestr(now,30)),'fits',...
  'T1cine','T2cine','im_mask','cphase','imstart')
