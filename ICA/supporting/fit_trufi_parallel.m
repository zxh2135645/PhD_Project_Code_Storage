
alpha = params.adFlipAngleDegree*pi/180;

e1 = @(R1)exp(-params.lEchoSpacing*R1);
e2 = @(R2)exp(-params.lEchoSpacing*R2);
e1star=@(e1,e2,alpha)e1*cos(alpha/2)^2+e2*sin(alpha/2)^2;
Mss = @(e1,e2,alpha)(1-e1)/(1-(e1-e2)*cos(alpha)-e1*e2);
n = 1:2:Nseg;
Sint = @(A,e1,e2,alpha,B)A * Mss(e1,e2,alpha) * (1 + (B-1)*e1star(e1,e2,alpha).^n) * sin(alpha);
% Sint = @(A,e,alpha,B)A * (1 + (B-1)*(e.^(n-1)));

if sliceprof
    S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B)+Sint(A,e1(R1),e2(R2),alpha/2,B);
else
    S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B);
end

ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)

if alpha0_deg == 90
    minB = 0;
    maxB = (1 - e(1/.01)^((params.alparams.lEchoSpacing_seconds-params.lEchoSpacing*Nseg)/2/params.lEchoSpacing))/(Mss(e(1/.01),alpha)); %assume 2 params.lEchoSpacing to first imaging pulse, shortest T1
elseif alpha0_deg == 180
    minB = -1;
    maxB = 1;
end

%%
h=figure;
imshow(abs(recon(:,:,end)),[]);
title('Select Fitting ROI')
im_mask=createMask(imrect);
close(h)
drawnow

im_mask=im_mask.*abs(recon(:,:,end));
[im_mask,centroids]=kmeans(im_mask(:),2);
im_mask=reshape(im_mask==1+(centroids(2)>centroids(1)),size(recon,1),size(recon,2));
figure,imshow(im_mask),drawnow;

[rows, cols]=ind2sub(size(im_mask),find(im_mask));
recontemp=reshape(recon,size(recon,1)*size(recon,2),[]);
recontemp=recontemp(im_mask(:),:);

opts=[];
opts.MaxFunEvals = 1000;

%%
fits = zeros(numel(rows),5);
parfor j=1:numel(rows)
    curve = double(recontemp(j,:));
    
    normcurve = curve(end);
    curve = curve/normcurve;
    
    Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve)+max(x(1)-x(2),0),...
        [2/3,10,alpha,median([minB, real(curve(1)/curve(end)), maxB])],...
        [0,0,0,minB], [10,50,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
    tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4))*normcurve, tempfit];
    
    fits(j,:) = tempfit;
    
end

tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),5);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),5);

%%
figure,imagesc(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imagesc(1./fits(:,:,3),[0 0.4]),colormap(jet(256))
figure,imagesc(fits(:,:,4)*180/pi),colormap(jet(256))

% alpha=input('alpha');

%%
fits = zeros(numel(rows),5);
parfor j=1:numel(rows)
    curve = double(recontemp(j,:));
    
    normcurve = curve(end);
    curve = curve/normcurve;
    
    Avp = @(R1,R2,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha,x(3))-curve)+max(x(1)-x(2),0),...
        [2/3,10,median([minB, real(curve(1)/curve(end)), maxB])],...
        [0,0,minB], [10,50,maxB], opts); 
    tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit(1), tempfit(2), alpha, tempfit(3)];
    
    fits(j,:) = tempfit;
    
end

tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),5);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),5);

%%
figure,imagesc(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imagesc(1./fits(:,:,3),[0 0.4]),colormap(jet(256))
figure,imagesc(fits(:,:,4)*180/pi),colormap(jet(256))