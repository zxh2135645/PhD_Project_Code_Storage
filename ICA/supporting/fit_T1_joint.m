alpha = params.adFlipAngleDegree*pi/180;

e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);

n = 1:2:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);
% Sint = @(A,e,alpha,B)A * (1 + (B-1)*(e.^(n-1)));


S = @(A,R1pre,R1post,alpha,Bpre,Bpost)[Sint(A,e(R1pre),alpha,Bpre) Sint(A,e(R1post),alpha,Bpost)];


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
fits = zeros(numel(rows),6);
parfor j=1:numel(rows)
    curve = double(recontemp(j,:));
    normcurve = curve(end);
    curve = curve/normcurve;
    
    curve(1:Nseg/2)=curve(1:Nseg/2)*exp(1i*(median(angle(curve(3*Nseg/4+1:end)))-median(angle(curve(Nseg/4+1:Nseg/2)))));
        
    Avp = @(R1pre,R1post,alpha,Bpre,Bpost)ppinv(S(1,R1pre,R1post,alpha,Bpre,Bpost),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5)),x(1),x(2),x(3),x(4),x(5))-curve),...
        [2/3,3/2,alpha,median([minB, real(curve(1)/curve(Nseg/2)), maxB]),median([minB, real(curve(Nseg/2+1)/curve(end)), maxB])],...
        [1/3,1/3,0,minB,minB], [1/.1,1/.1,alpha,maxB,maxB], opts); %be careful of alpha bounds. 0.58?
    tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5))*normcurve, tempfit];
    
    fits(j,:) = tempfit;
    
end

%%
tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),6);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),6);

figure,imagesc(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imagesc(1./fits(:,:,3),[0 2.3]),colormap(jet(256))
figure,imagesc(fits(:,:,3)-fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imagesc(fits(:,:,4)*180/pi),colormap(jet(256))
