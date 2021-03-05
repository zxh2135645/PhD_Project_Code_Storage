function [s_both, im_c_both, s_cbd, im_c_cbd, s_noreg]=sensemaps_cbd_ov(ims,kernel_size,varargin)
% Estimate smooth, normalized sensitivity maps (constrained blind deconvolution)
% v.0.3 by Anthony Christodoulou, 15 December 2016 - 20 December 2016
%
% sensemaps_sn(ims,kernel_size)
% sensemaps_sn(kspace,kernel_size,'kspace')
%
% ims or kspace:    Phase-encode X freq-encode X coils
%                or
%                   Phase-encode X freq-encode X partition-encode X coils
% kernel_size: convolution kernel size (default is [7 7] or [7 7 3])
%
% Options:
% 'kspace':         when inputting kspace data
% 'acs',kspace_mask:   when inputting incomplete kspace data with ACS lines
% 'figures':        to generate figures during calculation
% 'debug':          to output cost function information and more
%
% e.g., sensemaps_sn(ims,[8 8],'figures')
%
% Note: If using k-space data, data should preferably already be
%          fftshifted in k-space.
%       If using multislice 2D with large slice gaps, better to run this
%          one slice at a time, as the sensitivities may not be smooth in
%          the slice direction. Multislice mode may appear in a future
%          version.

if (numel(size(ims)) == 4) || (size(ims,3) == 1) %if 3D
  dims = 3;
  center_z = floor(size(ims,3)/2+1);
else %if 2D
  ims = reshape(ims,size(ims,1),size(ims,2),1,size(ims,3));
  dims = 2;
  center_z = 1;
end

if nargin < 2
  if dims==2
    kernel_size=[7 7];
  else
    kernel_size=[7 7 5];
  end
elseif isempty(kernel_size)
  if dims==2
    kernel_size=[7 7];
  else
    kernel_size=[7 7 5];
  end
elseif ischar(kernel_size) %if kernel_size omitted but options included
  varargin{end+1}=kernel_size;
  if dims==2
    kernel_size=[7 7];
  else
    kernel_size=[7 7 5];
  end
end

[Np, Nf, Ns, coils]=size(ims);
vec=@(x)x(:);
sos=@(s)sum(abs(reshape(s,Np,Nf,Ns,coils,[])).^2,4);

if dims == 2
  imshow3=@(x,caxis)montage(x,'DisplayRange',caxis);
  fftim = @(x)fft2(x);
  ifftim = @(x)ifft2(x);
else %3D
  imshow3=@(x,caxis)montage(x(:,:,center_z,:),'DisplayRange',caxis);
  fftim = @(x)fft(fft(fft(x,[],1),[],2),[],3);
  ifftim = @(x)ifft(ifft(ifft(x,[],1),[],2),[],3);
end

% ismember doesn't work well in newer MATLAB when using mixed-class cell
% arrays, so handle potential acs mask and k-space windows first
kspace_loc = find(strcmp(varargin,'acs')) + 1;
acs_flag = ~isempty(kspace_loc);
if acs_flag
  kspace_mask=varargin{kspace_loc};
  varargin(kspace_loc)=[];
end

kspace_flag=ismember('kspace',varargin); %all lowercase
figure_flag=ismember('figures',varargin); %all lowercase
debug_flag=ismember('debug',varargin); %all lowercase

if debug_flag
  kernel_size
  varargin
end

if kspace_flag
  kspace=ims;
else
  kspace=fftim(ifftshift(ifftshift(ifftshift(ims,1),2),3));
end

temp=sos(kspace);
temp_shift=ifftshift(temp);
if temp_shift(1) > temp(1) %if kspace isn't yet ifftshifted
  kspace = ifftshift(ifftshift(ifftshift(kspace,1),2),3);
  if acs_flag
    kspace_mask = ifftshift(kspace_mask);
  end
end
ims=fftshift(fftshift(fftshift(ifftim(kspace),1),2),3);
clear temp temp_shift

if acs_flag
  % Identify ACS region...assumes cube or rectangle at center of k-space
  rows=find(sum(sum(kspace_mask,2),3));
  if numel(rows)<Np
    rows=rows([1:find(diff(rows)>1,1), find(diff(rows)>1,1,'last')+1:end]);
  end
  
  cols=find(sum(sum(kspace_mask,1),3));
  if numel(cols)<Nf
    cols=cols([1:find(diff(cols)>1,1), find(diff(cols)>1,1,'last')+1:end]);
  end
  
  slices=find(sum(sum(kspace_mask,1),2));
  if numel(slices)<Ns
    slices=slices([1:find(diff(slices)>1,1), find(diff(slices)>1,1,'last')+1:end]);
  end
  
  kspace_mask=repmat(kspace_mask,[1 1 1 coils]);
  kspace=kspace_mask.*fftim(ims)/sqrt(Np*Nf*Ns);
  
  center_mask=zeros(Np,Nf,Ns);
  center_mask(rows,cols,slices)=1;
  center_mask=repmat(center_mask,[1 1 1 coils]);
  center_data=center_mask.*fftim(ims);
  
  window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ...
    .* repmat( reshape( hamming(numel(slices),'periodic'), 1,1,[]), [numel(rows), numel(cols), 1]);
  window=zeros(Np,Nf,Ns,coils);
  window(rows,cols,slices,:)=repmat(ifftshift(window_sm),[1 1 1 coils]);
  
  ims=ifftim(window.*center_data);
end

% Create coil combination image: sum-of-squares magnitude,
%                                first eigenimage phase
[~,~,V]=svd(reshape(ims,[],coils),0);
im_c=sqrt(sos(ims)).*sign(reshape(reshape(ims,[],coils)*V(:,1),Np,Nf,Ns));
s_noreg=bsxfun(@rdivide,ims,im_c);

if figure_flag
  if acs_flag
    im_orig=ifftim(kspace)*sqrt(Np*Nf*Ns);
    im_orig=sqrt(sos(im_orig)).*sign(reshape(reshape(im_orig,[],coils)*V(:,1),Np,Nf,Ns));
  else
    im_orig=im_c;
  end
  figure,subplot(1,2,1),imshow(abs(im_orig(:,:,center_z)),[]),title('SoS combination (Magnitude)')
  subplot(1,2,2),imshow(angle(im_orig(:,:,center_z)),[-pi pi]),title('1st eigenimage (Phase)')
  figure,imshow3(abs(s_noreg),[0 1]),title('Initial Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_noreg),[-pi pi]),title('Initial Sensitivity Maps (Phase)')
  drawnow
  clear im_orig
end

disp('Constrained blind deconvolution...')
hmask=false([Np Nf Ns coils]);
coords=floor(kernel_size/2);
if dims==2
  hmask([Np-coords(1)+1:Np, 1:coords(1)+1],[Nf-coords(2)+1:Nf, 1:coords(2)+1],1,:)=1;
  its=25;
else
  hmask([Np-coords(1)+1:Np, 1:coords(1)+1],[Nf-coords(2)+1:Nf, 1:coords(2)+1],[Ns-coords(3)+1:Ns, 1:coords(3)+1],:)=1;
  its=10;
end
s_cbd=s_noreg;

im_c_cbd=cat(5,im_c,zeros(size(im_c)));
s_cbd(:,:,:,:,2)=0; %temporarily

for it=1:its
  if acs_flag
    fwdh=@(h)bsxfun(@times,kspace_mask,fftim(sum(bsxfun(@times,im_c_cbd,ifftim(bsxfun(@times,reshape(h,Np,Nf,Ns,coils,2),hmask))),5)));
    adjh=@(x)vec(bsxfun(@times,fftim(bsxfun(@times,conj(im_c_cbd),ifftim(kspace_mask.*reshape(x,Np,Nf,Ns,coils)))),hmask));
    adjfwdh=@(h)adjh(fwdh(h));
    adjhb=adjh(kspace);
  else
    fwdh=@(h)sum(bsxfun(@times,im_c_cbd,ifftim(bsxfun(@times,reshape(h,Np,Nf,Ns,coils,2),hmask))),5)*sqrt(Np*Nf*Ns);
    adjh=@(x)vec(bsxfun(@times,fftim(bsxfun(@times,conj(im_c_cbd),reshape(x,Np,Nf,Ns,coils))),hmask))/sqrt(Np*Nf*Ns);
    adjfwdh=@(h)adjh(fwdh(h));
    adjhb=adjh(ims);
  end
  
  h=bsxfun(@times,fftim(s_cbd),hmask)/sqrt(Np*Nf*Ns);
  h=reshape(pcg(@(x)adjfwdh(x),adjhb,[],10,[],[],h(:)),Np,Nf,Ns,coils,2);
  
  s_cbd=ifftim(h)*sqrt(Np*Nf*Ns);
  s_cbd=bsxfun(@rdivide,s_cbd,sqrt(sos(s_cbd)));
  
  if (it==1)
    if acs_flag
      temp=ifftim(window.*(center_data-center_mask.*fftim(bsxfun(@times,reshape(im_c_cbd(:,:,:,:,1),Np,Nf,Ns),s_cbd(:,:,:,:,1)))));
    else
      temp=ims-bsxfun(@times,im_c_cbd(:,:,:,:,1),s_cbd(:,:,:,:,1));
    end
    s_cbd(:,:,:,:,2)=bsxfun(@rdivide,temp,sqrt(sos(temp)));
    s_cbd(:,:,:,:,2)=ifftim(bsxfun(@times,fftim(s_cbd(:,:,:,:,2)),hmask));
    s_cbd=bsxfun(@rdivide,s_cbd,sqrt(sos(s_cbd)));
  end
  
  if acs_flag
    A=@(x)bsxfun(@times,kspace_mask,fftim(sum(bsxfun(@times,reshape(x,Np,Nf,Ns,1,2)/sqrt(Np*Nf*Ns),s_cbd),5)));
    Ah=@(x)vec(sum(bsxfun(@times,conj(s_cbd),ifftim(sqrt(Np*Nf*Ns)*kspace_mask.*reshape(x,Np,Nf,Ns,coils))),4));
    AhA=@(x)Ah(A(x));
    im_c_cbd=reshape(pcg(AhA,Ah(kspace),[],10,[],[],im_c_cbd(:)),Np,Nf,Ns,1,2);
  else %is there a closed-form?
    A=@(x)sum(bsxfun(@times,reshape(x,Np,Nf,Ns,1,2),s_cbd),5);
    Ah=@(x)vec(sum(bsxfun(@times,conj(s_cbd),reshape(x,Np,Nf,Ns,coils)),4));
    AhA=@(x)Ah(A(x));
    im_c_cbd=reshape(pcg(AhA,Ah(ims),[],10,[],[],im_c_cbd(:)),Np,Nf,Ns,1,2);
  end
  
  if figure_flag && debug_flag
    ims_both=sum(bsxfun(@times,s_cbd,im_c_cbd),5);
    im_c_both=sqrt(sos(ims_both)).*sign(reshape(reshape(ims_both,[],coils)*V(:,1),Np,Nf,Ns));
    s_both=bsxfun(@rdivide,ims_both,im_c_both);
    
    try
      figure(fig(1))
    catch
      fig(1)=figure;
    end
    imshow3(abs(s_both),[0 1]),title('CBD Sensitivity Maps (Magnitude)')
    
    try
      figure(fig(2))
    catch
      fig(2)=figure;
    end
    imshow3(angle(s_both),[-pi pi]),title('CBD Sensitivity Maps (Phase)')
    
    try
      figure(fig(3))
    catch
      fig(3)=figure;
    end
    subplot(1,2,1),imshow(abs(im_c_both(:,:,center_z)),[0 max(abs(im_c(:)))]),title('Cartesian SENSE with CBD maps')
    subplot(1,2,2),imshow(angle(im_c_both(:,:,center_z)),[-pi pi])
    drawnow
  end
end

if figure_flag
  ims_both=sum(bsxfun(@times,s_cbd,im_c_cbd),5);
  im_c_both=sqrt(sos(ims_both)).*sign(reshape(reshape(ims_both,[],coils)*V(:,1),Np,Nf,Ns));
  s_both=bsxfun(@rdivide,ims_both,im_c_both);
  figure,imshow3(abs(s_both),[0 1]),title('CBD Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_both),[-pi pi]),title('CBD Sensitivity Maps (Phase)')
  figure,subplot(1,2,1),imshow(abs(im_c_both(:,:,center_z)),[0 max(abs(im_c(:)))]),title('Cartesian SENSE with CBD maps')
  subplot(1,2,2),imshow(angle(im_c_both(:,:,center_z)),[-pi pi])
  drawnow
end

return