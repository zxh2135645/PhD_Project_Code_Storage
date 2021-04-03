
alpha = params.adFlipAngleDegree*pi/180;

e1 = @(R1)exp(-params.lEchoSpacing*R1);
e2 = @(R2)exp(-[12 20 30 40 50]*1e-3*R2);
Mss = @(e1,alpha)(1-e1) / (1-cos(alpha)*e1);
n = params.lSegments-Nseg+1:2:params.lSegments;
Sint = @(A,e1,e2,alpha,B)vec(A * Mss(e1,alpha) * (1 + ((e1*cos(alpha)).^(n-1)).'*(B*e2-1)) * sin(alpha)).';
S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B);

ppinv=@(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)

minB = -1;
maxB = 1;


%%
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

%%
% fits = zeros(numel(rows),5);
% parfor j=1:numel(rows)
%     curve = double(recontemp(j,:));
%
%     normcurve = curve(end);
%     curve = curve/normcurve;
%
%     Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
%     tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve),...
%         [2/3,10,alpha,median([minB, real(curve(1)/curve(end)), maxB])],...
%         [1/3,1,0,minB], [1/.1,1/.01,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
%     tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4))*normcurve, tempfit];
%
%     fits(j,:) = tempfit;
%
% end
%
% tempfit=fits;
% fits=zeros(size(recon,1)*size(recon,2),5);
% fits(im_mask(:),:)=tempfit;
% fits=reshape(fits,size(recon,1),size(recon,2),5);
%
% %%
% figure,imagesc(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
% figure,imagesc(1./fits(:,:,3),[0 0.1]),colormap(jet(256))
% figure,imagesc(fits(:,:,4)*180/pi),colormap(jet(256))
%
% alpha=input('alpha? ');

% if exist('sept_mask','var') %if sept_mask is available
%   %Fit myocardial B1
%   curve=reshape(recon,[],size(recon,3));
%   [~,~,curve]=svd(curve(sept_mask,:),'econ');
%   curve=curve(:,1).';
%   normcurve = curve(end);
%   curve = curve/normcurve;
%   Avp = @(R1,R2,alpha,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
%   tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3),x(4)),x(1),x(2),x(3),x(4))-curve),...
%     [1/1.1,20,alpha/2,median([minB, real(curve(1)/curve(end)), maxB])],...
%     [1/3,1,0,minB], [1/.1,1/.01,alpha,maxB], opts); %be careful of alpha bounds. 0.58?
%   myoT1=1/tempfit(1)
%   myoT2=1/tempfit(2)
%   myoalpha=tempfit(3)*180/pi
%   alpha=tempfit(3);
% else
  alpha=2.5*pi/180
% end

%%
fits = zeros(numel(rows),5);
parfor j=1:numel(rows)
  curve = double(recontemp(j,:));
  
  normcurve = curve(end);
  curve = curve/normcurve;
  
  Avp = @(R1,R2,B)ppinv(S(1,R1,R2,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
  tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),alpha,x(3))-curve),...
    [2/3,10,median([minB, real(curve(1)/curve(end)), maxB])],...
    [1/3,1,minB], [1/.1,1/.01,maxB], opts);
  tempfit=[Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve, tempfit(1), tempfit(2), alpha, tempfit(3)];
  
  fits(j,:) = tempfit;
end

tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),5);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),5);

%%
figure,imagesc(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imagesc(1./fits(:,:,3),[0 0.1]),colormap(jet(256))
figure,imagesc(fits(:,:,5)),colormap(jet(256))