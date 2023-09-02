%% Select IR-FLASH recon

if ispc
    sep = '\';
else
    sep = '/';
end

folders=dir;
found=[];
for j=1:numel(folders);
    if folders(j).isdir
        file=dir(sprintf('%s%sAC_recon.mat',folders(j).name,sep));
        if numel(file) > 0;
        found{end+1}=folders(j).name;
        end
    end
end

folder=listdlg('ListString',found,'SelectionMode','single','ListSize',[500 150]);
load(sprintf('%s%sAC_recon.mat',found{folder},sep),'U','Phi','L','Gr',...
    'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec')

%% Choose recovery images
Phi=reshape(Phi,[L sizes(2:end)]);

temp=Gr\reshape(Phi(:,1,:,1,:),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp1=abs(dispim(temp))/max(abs(vec(dispim(temp))));

temp=Gr\reshape(Phi(:,71,:,1,:),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp2=abs(dispim(temp))/max(abs(vec(dispim(temp))));

temp=Gr\reshape(Phi(:,111,:,1,:),L,[]);
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
temp3=abs(dispim(temp))/max(abs(vec(dispim(temp))));

implay(cat(2,temp1,temp2,temp3));
clear temp*;

cphase=input('Which phase is diastole [1]? ');
if isempty(cphase)
    cphase=1;
end

recon = Gr\reshape(Phi(:,:,cphase,1,ceil(end/2)),L,[]);
recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
recon=reshape(recon,N-ovs,N-ovs,Nseg/2);

implay(abs(recon)/max(abs(recon(:))))

imstart=input('Choose starting image [1]? ');
if isempty(imstart)
    imstart=1;
end
recon=recon(:,:,imstart:end);
Nseg=2*size(recon,3);

%% Choose MOLLI images

MOLLIfolders=dir;
for j=1:numel(MOLLIfolders);
    if MOLLIfolders(j).isdir
        files=dir(sprintf('%s%s*.0001.*',MOLLIfolders(j).name,sep));
        if numel(files) > 0;
            DICOM_folder=j;
            break;
        end
    end
end

MOLLI_list=[];
MOLLIs=[];
for j=1:numel(files)
    info=dicominfo(strcat(MOLLIfolders(DICOM_folder).name,sep,files(j).name));
    isMOLLI_T1=~isempty(regexp(info.SeriesDescription,'MOCO_T1', 'once'));
    isMOLLI_T1S=~isempty(regexp(info.SeriesDescription,'MOCO_T1S', 'once'));
    if isMOLLI_T1 && ~isMOLLI_T1S
        MOLLI_list(end+1)=j;
        try
        MOLLIs(:,:,:,numel(MOLLI_list))=ind2rgb(uint8(255*double(dicomread(strcat(MOLLIfolders(DICOM_folder).name,sep,files(j).name)))/2300),jet(256));
        catch
            MOLLIs(:,:,:,numel(MOLLI_list))=ind2rgb(uint8(255*double(dicomread(strcat(MOLLIfolders(DICOM_folder).name,sep,files(j).name)))/2300).',jet(256));
        end
    end
end

implay(MOLLIs)
MOLLI=input('Choose MOLLI reference [1]? ');
if isempty(MOLLI)
    MOLLI=1;
end
T1_ref=double(dicomread(strcat(MOLLIfolders(DICOM_folder).name,sep,files(MOLLI_list(MOLLI)).name)))/1000;
crop=size(T1_ref)-N+ovs;
croprect=[flip(crop,2)/2 + 1, N - ovs - 1, N - ovs - 1];
T1_ref=imcrop(T1_ref,croprect);

h=figure;
imshow(T1_ref,[0 2.3]),colormap(jet(256))
title('Draw Septum')
myo_roi = impoly(gca,'Closed',1);
myo_mask = createMask(myo_roi);
close(h)
T1 = median(T1_ref(myo_mask))

%% Fit
sliceprof=false; %true;
alpha0_deg=180;

alpha = params.adFlipAngleDegree*pi/180;

e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);

n = 1:2:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);
% Sint = @(A,e,alpha,B)A * (1 + (B-1)*(e.^(n-1)));

if sliceprof
    S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B)+Sint(A,e(R1),alpha/2,B);
else
    S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);
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
fits = zeros(numel(rows),4);
parfor j=1:numel(rows)
    curve = double(recontemp(j,:));
    
    normcurve = curve(end);
    curve = curve/normcurve;
    
    Avp = @(alpha,B)ppinv(S(1,1/T1,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2)),1/T1,x(1),x(2))-curve),...
        [alpha/2,median([minB, real(curve(1)/curve(end)), maxB])],...
        [0,minB], [alpha,maxB], opts); %be careful of alpha bounds. 0.58?
    tempfit=[Avp(tempfit(1),tempfit(2))*normcurve, 1/T1, tempfit];
    
    fits(j,:) = tempfit;
    
end

tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),4);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),4);

%%
B1map=fits(:,:,3)/alpha;
h=figure;
imshow(B1map),colormap(jet(256))
title('Draw Septum')
myo_roi = impoly(gca,'Closed',1);
myo_mask = createMask(myo_roi);
close(h)
alpha = alpha*median(B1map(myo_mask));

fid=fopen('alpha','w')
fprintf(fid,'%0.2f degrees',alpha*180/pi);
fclose(fid);
printf('%0.2f degrees',alpha*180/pi)

%%
fits = zeros(numel(rows),4);
parfor j=1:numel(rows)
    curve = double(recontemp(j,:));
    
    normcurve = curve(end);
    curve = curve/normcurve;
    
    Avp = @(R1,B)ppinv(S(1,R1,alpha,B),curve); %parameterize solution to A as function of R1,alpha,B
    tempfit=lsqnonlin(@(x)abs(S(Avp(x(1),x(2)),x(1),alpha,x(2))-curve),...
        [2/3,median([minB, real(curve(1)/curve(end)), maxB])],...
        [1/3,minB], [1/.1,maxB], opts); 
    tempfit=[Avp(tempfit(1),tempfit(2))*normcurve, tempfit(1), alpha, tempfit(2)];
    
    fits(j,:) = tempfit;
    
end

tempfit=fits;
fits=zeros(size(recon,1)*size(recon,2),4);
fits(im_mask(:),:)=tempfit;
fits=reshape(fits,size(recon,1),size(recon,2),4);

%%
figure,imshow(1./fits(:,:,2),[0 2.3]),colormap(jet(256))
figure,imshow(T1_ref,[0 2.3]),colormap(jet(256))