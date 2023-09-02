function SEs=espirit_agc(composite_k,winlen,zstack_flag)

if nargin<2
  winlen=5
  zstack_flag=false;
elseif nargin <3
  zstack_flag=false;
end
winlen=floor(winlen/2)*2+1;

[Ny,Nx,Nz,coils]=size(composite_k);

halfwin=floor(winlen/2);
windex=-halfwin:halfwin;
if (Nz>1) && ~zstack_flag
  zindex=windex;
else
  zindex=0;
end

DC_fe=floor(Nx/2+1);
DC_ky=floor(Ny/2+1);
DC_kz=floor(Nz/2+1);

GP.mask = composite_k(:,DC_fe,:,1)~=0;
GP.Ny2 = find(ifftshift(composite_k(:,DC_fe,DC_kz,1)==0),1)-1;
GP.Ny1 = Ny-find(ifftshift(composite_k(:,DC_fe,DC_kz,1)==0),1,'last');
GP.Nyw = GP.Ny1+GP.Ny2;
if ~isempty(GP.Nyw)
  GP.Nyw;
  GP.Nyw_indices = (DC_ky-GP.Ny1):(DC_ky+GP.Ny2-1);
else
  GP.Nyw = Ny;
  GP.Nyw_indices = 1:Ny;
end

if zstack_flag
  kk_ind=1:Nz;
elseif Nz>1
  kk_ind=(halfwin+1):(Nz-halfwin);
else
  kk_ind=1;
end

if isunix
  [~,bytes]=unix('/usr/local/bin/free -b | sed -e "1d;2d;4d" | awk ''{print $4}''');
  bytes=str2num(bytes);
else
  [~,bytes]=memory;
  bytes=bytes.PhysicalMemory.Available;
end
double_entries=bytes/8;
entries_required = (Nx-2*halfwin)*(GP.Nyw-2*halfwin)*numel(kk_ind)*winlen^2*numel(zindex)*coils*2;
if ( entries_required > double_entries )
  kk_ind_w = floor(numel(kk_ind)*double_entries/entries_required/4);
  fprintf('Warning: reducing from %d kz-lines to %d.\n',numel(kk_ind),kk_ind_w*2+1);
  kk_ind = kk_ind(floor(numel(kk_ind)/2+1)) + (-kk_ind_w:kk_ind_w);
end

A=complex(zeros((Nx-2*halfwin)*(GP.Nyw-2*halfwin)*numel(kk_ind),winlen^2*numel(zindex)*coils));
ksize=winlen^2*numel(zindex);
index=0;
for kk=zindex
  for jj=windex
    for ii=windex
      A(:,index + (1:ksize:(ksize*coils)))=reshape(permute(composite_k(GP.Nyw_indices((halfwin+1):(GP.Nyw-halfwin))+ii,((halfwin+1):(Nx-halfwin))+jj,kk_ind+kk,:),[2 1 3 4]),[],coils);
      index=index+1;
    end
  end
end

disp('ESPiRIT SVD...')
% [~,S,V]=svd(A,'econ');
[~,S,V]=svde(A);

figure,plot(20*log10(diag(S)/S(1)))
% lasteig=find(diag(S)/S(1)>0.01,1,'last')
lasteig=input('Last eigenvalue: ');
Vp=V(:,1:lasteig);
clear A S V;

if zstack_flag
  Nz=1;
  DC_kz=1;
end

if isunix
  [~,bytes]=unix('/usr/local/bin/free -b | sed -e "1d;2d;4d" | awk ''{print $4}''');
  bytes=str2num(bytes);
else
  [~,bytes]=memory;
  bytes=bytes.PhysicalMemory.Available;
end
double_entries=bytes/8;
max_mem=floor(double_entries/(Ny*Nx*Nz*coils*2)); %COMPLEX no safety factor
if lasteig > max_mem
  disp('Not enough memory, reducing subspace dimension...')
  lasteig=max_mem
  Vp=Vp(:,1:lasteig);
end

G=complex(zeros(Ny,Nx,Nz,coils,size(Vp,2)));
G(DC_ky+windex,DC_fe+windex,DC_kz+zindex,:,:)=reshape(Vp,winlen,winlen,numel(zindex),coils,[]);
G=fft(fft(fft(G,[],1),[],2),[],3);

% eigval=zeros(Ny,Nx,Nz);
eigvec=zeros(Ny,Nx,Nz,coils);
% pool=parpool(6);
for qr=1:Ny
  for qc=1:Nx
    for qs=1:Nz
      Gq=squeeze(G(qr,qc,qs,:,:)).';
%       [~,~,v]=svd(Gq,'econ');
        [~,~,v]=svde(Gq);
%       eigval(qr,qc,qs,:)=s(1);
      eigvec(qr,qc,qs,:)=v(:,1);
    end
  end
end
% delete(pool)

SEs=ifft(ifft(ifft(...
  ifftshift(ifftshift(ifftshift(...
  fft(fft(fft(eigvec,[],1),[],2),[],3)...
  ,1),2),3)...
  ,[],1),[],2),[],3);

end
