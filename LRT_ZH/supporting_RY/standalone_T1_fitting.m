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
% 
% folders=dir;
% found=[];
% k=0;
% for j=1:numel(folders);
%     if folders(j).isdir
%         files=dir(sprintf('%s%s*.0001.*',folders(j).name,sep));
%         if numel(file) > 0;
%         break;
%         end
%     end
% end


%% Fit
sliceprof=false; %true;
alpha0_deg=180;
fit_T1_parallel
