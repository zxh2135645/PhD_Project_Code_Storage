function [s_walsh, s_noreg]=sensemaps_walsh(ims,window_length,varargin)
% Estimate sensitivity maps via Walsh method
% v1.0.0 by Anthony Christodoulou, 17 October 2016
%
% sensemaps_sn(ims,window_length)
% sensemaps_sn(kspace,'kspace')
%
% ims or kspace:    Phase-encode X freq-encode X coils
%                or 
%                   Phase-encode X freq-encode X partition-encode X coils
% window_length:    patch size
%
% Options:
% 'kspace':         when inputting kspace data
% 'acs',ACS_mask:   when inputting incomplete kspace data with ACS lines
% 'window',window:  optional k-space weightings (e.g., DCF for non-Cartesian
%                   acquisition). Currently incompatible with 'acs' mode.
% 'figures':        to generate figures during calculation
% 'debug':          to output cost function information and more
% 
% e.g., sensemaps_sn(ims,'figures')
%
% Note: If using k-space data, data should preferably already be
%          fftshifted in k-space, but not in image space.
%       If using 'acs' or 'window' with 'figures', coil combination will be
%          windowed rather than full-resolution
%       If using multislice 2D with large slice gaps, better to run this
%          one slice at a time, as the sensitivities may not be smooth in
%          the slice direction. Multislice mode may appear in a future
%          version.

if nargin < 2
  window_length=5
elseif isempty(window_length)
  window_length=5
elseif ischar(window_length) %if window_length omitted but options included
  varargin{end+1}=window_length;
  window_length=5
end

if (numel(size(ims)) == 4) || (size(ims,3) == 1) %if 3D
  dims = 3;
  center_z = floor(size(ims,3)/2+1);
else %if 2D
  ims = reshape(ims,size(ims,1),size(ims,2),1,size(ims,3));
  dims = 2;
  center_z = 1;
end

[Np, Nf, Ns, coils]=size(ims);
sizes = size(ims);
vec=@(x)x(:);
sos=@(s)sum(abs(reshape(s,sizes(1),sizes(2),[],sizes(4))).^2,4);

if dims == 2
  imshow3=@(x,caxis)montage(x,caxis);
  fftim = @(x)fft2(x);
  ifftim = @(x)ifft2(x);
else %3D
  imshow3=@(x,caxis)montage(x(:,:,center_z,:),caxis);
  fftim = @(x)fft(fft(fft(x,[],1),[],2),[],3);
  ifftim = @(x)ifft(ifft(ifft(x,[],1),[],2),[],3);
end

% ismember doesn't work well in newer MATLAB when using mixed-class cell
% arrays, so handle potential acs mask and k-space windows first
acs_loc = find(strcmp(varargin,'acs')) + 1;
acs_flag = ~isempty(acs_loc);
if acs_flag
  center_mask=varargin{acs_loc};
  varargin(acs_loc)=[];
end

window_loc = find(strcmp(varargin,'window')) + 1;
window_flag = ~isempty(window_loc);
if window_flag
  window=varargin{window_loc};
  window=window/max(abs(window(:))); %normalize
  varargin(window_loc)=[];
end

kspace_flag=ismember('kspace',varargin); %all lowercase
figure_flag=ismember('figures',varargin); %all lowercase
debug_flag=ismember('debug',varargin); %all lowercase

if debug_flag
  varargin
end
  
if kspace_flag
  kspace=ims;
else
  kspace=fftim(ifftshift(ims));
end

temp=sos(kspace);
temp_shift=ifftshift(temp);
if temp_shift(1) > temp(1) %if kspace isn't yet ifftshifted
  kspace = ifftshift(ifftshift(ifftshift(kspace,1),2),3);
  if acs_flag
    center_mask = ifftshift(center_mask);
  end
  if window_flag
    window = ifftshift(window);
  end
end
ims=fftshift(ifftim(kspace));
clear temp temp_shift

if window_flag
  window=repmat(window,[1 1 1 coils]);
  
  center_mask=abs(window)>0;

  center_data=window.*center_mask.*fftim(ims);

  window_comp=ones(size(window)); %no need to compensate for window;
  
  ims=sqrt(Np*Nf*Ns)*ifftim(center_data); %use normalized fft here
elseif acs_flag
  % Identify ACS region...assumes cube or rectangle at center of k-space  
  rows=find(sum(sum(center_mask,2),3));
  cols=find(sum(sum(center_mask,1),3));
  slices=find(sum(sum(center_mask,1),2));
  
  center_mask=repmat(center_mask,[1 1 1 coils]);
  center_data=center_mask.*fftim(ims);
  
  window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ...
    .* repmat( reshape( hamming(numel(slices),'periodic'), 1,1,[]), [numel(rows), numel(cols), 1]);
  window=zeros(Np,Nf,Ns,coils);
  window(rows,cols,slices,:)=repmat(ifftshift(window_sm),[1 1 1 coils]);

  window_comp=1./window; %window compensation
  window_comp(isinf(window_comp))=0;
  ims=sqrt(Np*Nf*Ns)*ifftim(window.*center_data); %use normalized fft here
end

%Normalize ims to similar scale everytime. Does not affect scaling of
% sensitivity maps.
if acs_flag || window_flag
  center_data=center_data/norm(ims(:))*coils;
end
ims=ims/norm(ims(:))*coils;

% Create coil combination image: sum-of-squares magnitude,
%                                first eigenimage phase
[U,~,mixer]=svd(reshape(ims,[],coils),0);
im_c=reshape(U(:,1),Np,Nf,Ns);
im_c=sqrt(sos(ims)).*exp(1i*angle(im_c));
if figure_flag
  figure,subplot(1,2,1),imshow(abs(im_c(:,:,center_z)),[]),title('SoS combination (Magnitude)')
  subplot(1,2,2),imshow(angle(im_c(:,:,center_z)),[-pi pi]),title('1st eigenimage (Phase)')
  drawnow
end

% Solve arg min_s || ims - A(s) ||
% noreg: no regularization
s_noreg=ims./repmat(im_c,[1 1 1 coils]); %closed-form
if figure_flag
  figure,imshow3(abs(s_noreg),[0 1]),title('Naive Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_noreg),[-pi pi]),title('Naive Sensitivity Maps (Phase)')
  drawnow
end


% do image-space blockwise SVD
s_walsh=s_noreg;
lb=floor(window_length/2);
win=@(x,xmax)mod(x+(-lb:lb)-1,xmax)+1;
if dims == 2
  winz=@(x,xmax)1;
else
  winz=win;
end
disp('SVDs...')
for j=1:Np
  for k=1:Nf
    for l=1:Ns
      [~,~,V]=svde(reshape(ims(win(j,Np),win(k,Nf),winz(l,Ns),:),[],coils));
      s_walsh(j,k,l,:)=conj(V(:,1))*exp(1i*angle(mean(V(:,1))./mean(mixer(:,1))));
    end
  end
end

if figure_flag
  figure,imshow3(abs(s_walsh),[0 1]),title('Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_walsh),[-pi pi]),title('Sensitivity Maps (Phase)')
  
%   figure,imshow3(abs(s_walsh-s_noreg),[0 1]),title('Change in Sensitivity Maps (Magnitude)')
%   figure,imshow3(angle(s_walsh./s_noreg),[-pi pi]),title('Change in Sensitivity Maps (Phase)')
  
  im_c_walsh=sum(ims.*conj(s_walsh),4)./sos(s_walsh);
  figure,subplot(1,2,1),imshow(abs(im_c_walsh(:,:,center_z)),[0 max(abs(im_c(:)))]),title('Coil combination with smooth maps')
  subplot(1,2,2),imshow(angle(im_c_walsh(:,:,center_z)),[-pi pi])
  figure,imshow(sqrt(sos(s_walsh(:,:,center_z,:))),[]),caxis(caxis.*[0 1]),title('Relative SNR (Sensitivity map combination)')
  drawnow
end

return