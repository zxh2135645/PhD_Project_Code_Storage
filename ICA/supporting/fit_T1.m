
alpha = params.adFlipAngleDegree*pi/180;

e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);

n = 1:2:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);

if sliceprof
  S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B)+Sint(A,e(R1),alpha/2,B);
else
  S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);
end

ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)

if alpha0_deg == 90
  minB = 0;
  maxB = (1 - e(1/.01)^((params.alparams.lEchoSpacing_seconds-params.lEchoSpacing*Nseg)/2/params.lEchoSpacing))/(Mss(e(1/.01),alpha)); %assume 2 params.lEchoSpacing to first imaging pulse, shortest T1
  % maxB = 1;
elseif alpha0_deg == 180
  %   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
  minB = -1;
  maxB = 1;
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^2)/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha));
%   maxB = max(maxB,-.5)
end

h=figure;
imshow(abs(recon(:,:,end,1)),[]);
title('Select Fitting ROI')
im_mask=createMask(imrect);
close(h)
drawnow

im_mask=im_mask.*abs(recon(:,:,end,1));
[im_mask,centroids]=kmeans(im_mask(:),2);
im_mask=reshape(im_mask==1+(centroids(2)>centroids(1)),size(recon,1),size(recon,2));
figure,imshow(im_mask),drawnow;
% im_mask = true(N,N);

% fat_mask=(real(recon(:,:,46,1)./recon(:,:,end,1))>0 & real(recon(:,:,1,1)./recon(:,:,end,1))<0 & im_mask);

[rows, cols]=ind2sub(size(im_mask),find(im_mask));
frames = size(recon,4);
fits = zeros(size(recon,1),size(recon,2),4,frames);
% fitsmodif=fits;

opts=[];
opts.MaxFunEvals = 1000;
for frame=1:frames
  for j=1:numel(rows)
    row = rows(j);
    col = cols(j);
    curve = double(squeeze(recon(row,col,:,frame))).';
    %   [u,s,v]=svd([real(curve(:)), imag(curve(:))],'econ');
    %   normcurve = s(1)*complex(v(1),v(2))*u(end,1);
    %   curve = u(:,1).'/u(end,1);
    normcurve = curve(end);
    curve = curve/normcurve;
    
    %   Avp = @(R1,alpha,B)curve*pinv(S(1,R1,alpha,B));
    Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),x(3))-curve),...
      [1.5,alpha,median([minB, real(curve(1)/curve(end)), maxB])],...
      [1/3,0,minB], [1/.1,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
    fits(row,col,1,frame) = Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve;
    fits(row,col,2:4,frame) = tempfit;
    
%     Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B)*modif,curve*modif);
%     tempfit=lsqnonlin(@(x)abs((S(Avp(x(1),x(2),x(3)),x(1),x(2),x(3))-curve)*modif),...
%       [1.5,alpha,median([minB, real(curve(1)/curve(end)), maxB])],...
%       [1/3,alpha/2,minB], [1/.1,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
%     fitsmodif(row,col,1,frame) = Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve;
%     fitsmodif(row,col,2:4,frame) = tempfit;
    
  end
end