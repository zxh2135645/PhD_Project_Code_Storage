function [s_reg, s_reg_unc, s_noreg, s_reg_proj]=sensemaps_sn(ims,lambda,varargin)
% Estimate smooth, normalized sensitivity maps
% v2.0.1 by Anthony Christodoulou, 12 October 2015
%
% sensemaps_sn(ims,lambda)
% sensemaps_sn(kspace,lambda,'kspace')
%
% ims or kspace:    Phase-encode X freq-encode X coils
%                or 
%                   Phase-encode X freq-encode X partition-encode X coils
% lambda: Smoothness parameter
%
% Options:
% 'kspace':         when inputting kspace data
% 'acs',ACS_mask:   when inputting incomplete kspace data with ACS lines
% 'window',window:  optional k-space weightings (e.g., DCF for non-Cartesian
%                   acquisition). Currently incompatible with 'acs' mode.
% 'figures':        to generate figures during calculation
% 'debug':          to output cost function information and more
% 
% e.g., sensemaps_sn(ims,lambda,'figures')
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
  lambda=1
elseif isempty(lambda)
  lambda=1
elseif ischar(lambda) %if lambda omitted but options included
  varargin{end+1}=lambda;
  lambda=1
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
  lambda
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
[U,~,~]=svd(reshape(ims,[],coils),0);
im_c=reshape(U(:,1),Np,Nf,Ns);
im_c=sqrt(sos(ims)).*exp(1i*angle(im_c));
if figure_flag
  figure,subplot(1,2,1),imshow(abs(im_c(:,:,center_z)),[]),title('SoS combination (Magnitude)')
  subplot(1,2,2),imshow(angle(im_c(:,:,center_z)),[-pi pi]),title('1st eigenimage (Phase)')
  drawnow
end

if acs_flag || window_flag
  % "A" applies the coil combination image to the sense maps to generate
  %   truncated, sensitivity-encoded k-space data
  A=@(s)vec(window_comp.*center_mask.*(fftim(repmat(im_c,[1 1 1 coils]).*reshape(s,sizes))/sqrt(Np*Nf*Ns)));
  Ah=@(x)vec(repmat(conj(im_c),[1 1 1 coils]).*ifftim(conj(window_comp).*center_mask.*reshape(x,sizes))*sqrt(Np*Nf*Ns));
  AhA=@(s)vec(repmat(conj(im_c),[1 1 1 coils]).*ifftim((abs(window_comp).^2).*center_mask.*fftim(repmat(im_c,[1 1 1 coils]).*reshape(s,sizes))));
  Ahb=Ah(center_data);
else
  % "A" applies the coil combination image to the sense maps to generate
  %   sensitivity-encoded images
  A=@(s)vec(repmat(im_c,[1 1 1 coils]).*reshape(s,sizes));
  Ah=@(x)vec(repmat(conj(im_c),[1 1 1 coils]).*reshape(x,sizes));
  AhA=@(s)vec(repmat(abs(im_c).^2,[1 1 1 coils]).*reshape(s,sizes));
  Ahb=Ah(ims);
end

% s_noreg=reshape(pcg(AhA,Ah(center_data),[],500),size(im_0)); %iterative, from 0

% Solve arg min_s || ims - A(s) ||
% noreg: no regularization
s_noreg=ims./repmat(im_c,[1 1 1 coils]); %closed-form
if figure_flag
  figure,imshow3(abs(s_noreg),[0 1]),title('Unregularized Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_noreg),[-pi pi]),title('Unregularized Sensitivity Maps (Phase)')
  drawnow
end

% build gradient functions
padvec=@(dim)(sizes - (sizes(dim) - 1) * ((1:4)==dim));
pad=@(x,dim)cat(dim,zeros(padvec(dim)),x,zeros(padvec(dim)));

G1 = @(x,dim)diff(reshape(x,sizes),1,dim);
G1h = @(x,dim)-diff(pad(x,dim),1,dim);
if dims == 2
  G = @(x){G1(x,1), G1(x,2)};
  Gh = @(x)G1h(x{1},1) + G1h(x{2},2);
else
  G = @(x){G1(x,1), G1(x,2), G1(x,3)};
  Gh = @(x)G1h(x{1},1) + G1h(x{2},2) + G1h(x{3},3);
end
GhG = @(x)Gh(G(x));

% Solve arg min_s || ims - A(s) ||^2 + lambda * || G(s) ||^2
% reg_unc: regularized, unconstrained minimization
disp('Unconstrained minimization...')
s_reg_unc=zeros(size(s_noreg)); %initial guess
if acs_flag || window_flag
  s_reg_unc = s_noreg; % should already be pretty smooth
end
% s_reg_unc=imfilter(s_noreg,fspecial('gaussian',[10 10],5)); %other option

s_reg_unc=ifftim(fftim(s_noreg).*ifftshift(repmat(...
    reshape(kron(kron(hamming(Ns,'periodic'),hamming(Nf,'periodic')),hamming(Np,'periodic')),Np,Nf,Ns)...
    ,[1 1 1 coils])));
tic;
s_reg_unc=reshape(pcg(@(s)(AhA(s)+lambda*vec(GhG(s))),Ahb,[],50,[],[],s_reg_unc(:)),sizes); %iterative, from 0
toc;

if figure_flag
  figure,imshow3(abs(s_reg_unc),[0 1]),title('Smooth Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_reg_unc),[-pi pi]),title('Smooth Sensitivity Maps (Phase)')
  
%   figure,imshow3(abs(s_reg_unc-s_noreg),[0 1]),title('Change in Sensitivity Maps (Magnitude)')
%   figure,imshow3(angle(s_reg_unc./s_noreg),[-pi pi]),title('Change in Sensitivity Maps (Phase)')
  
  im_c_reg_unc=sum(ims.*conj(s_reg_unc),4)./sos(s_reg_unc);
  figure,subplot(1,2,1),imshow(abs(im_c_reg_unc(:,:,center_z)),[0 max(abs(im_c(:)))]),title('Coil combination with smooth maps')
  subplot(1,2,2),imshow(angle(im_c_reg_unc(:,:,center_z)),[-pi pi])
  figure,imshow(sqrt(sos(s_reg_unc(:,:,center_z,:))),[]),caxis(caxis.*[0 1]),title('Relative SNR (Sensitivity map combination)')
  drawnow
end

% Solve arg min_s || ims - A(s) ||^2 + lambda * || G(s) ||^2
%                                       s.t. sum_c |s_c(r)|^2 = 1 for all r
%
% use ADMM: arg min_s,z || ims - A(s) ||^2 + lambda * || G(s) ||^2
%                                       s.t. z in set of normalized maps, 
%                                            s = z
disp('Constrained minimization...')

%cost functions
if acs_flag || window_flag
  if dims == 2
    cost_quad=@(s,Gs)(norm(A(s)-center_data(:))^2+lambda*norm([Gs{1}(:); Gs{2}(:)])^2);
  else
    cost_quad=@(s,Gs)(norm(A(s)-center_data(:))^2+lambda*norm([Gs{1}(:); Gs{2}(:); Gs{3}(:)])^2);
  end
else
  if dims == 2
    cost_quad=@(s,Gs)(norm(A(s)-ims(:))^2+lambda*norm([Gs{1}(:); Gs{2}(:)])^2);
  else
    cost_quad=@(s,Gs)(norm(A(s)-ims(:))^2+lambda*norm([Gs{1}(:); Gs{2}(:); Gs{3}(:)])^2);
  end
end
cost_lagrange=@(s,y,z,rho)rho/2*norm(s(:)-z(:))^2+real(y(:)'*(s(:)-z(:)));
cost_equality=@(s)norm(vec(sos(s)-1))^2;

projs = @(s)s./sqrt(repmat(sos(s),[1 1 1 coils]));
s_reg_proj=projs(s_reg_unc);

s_reg = s_reg_unc;
s_reg_old = s_reg;
y=zeros(size(s_reg));
z=projs(s_reg);

costs=[cost_quad(s_reg,G(s_reg)), cost_equality(s_reg) cost_lagrange(s_reg,y,z,1)];
rho = costs(1)/costs(3);
costs(1,3)=rho*costs(1,3);

eps = inf;
i=0;
while (eps > 1e-3) && (i<10)
  i=i+1;
  z = projs(s_reg+y/rho);
  z(isnan(z))=0;
  y = y + rho*(s_reg-z);
  rho=rho*1.1;
  s_reg=reshape(pcg(@(s)(AhA(s)+lambda*vec(GhG(s))+rho/2*s(:)),Ahb+rho/2*(z(:)-y(:)/rho),[],min(i*10,50),[],[],z(:)),sizes);
  if debug_flag
    rho
    costs=[costs; [cost_quad(s_reg,G(s_reg)), cost_equality(s_reg) cost_lagrange(s_reg,y,z,rho)]];
    if figure_flag
      figure(99)
      subplot(4,1,1),plot(costs(:,1)),title('Quadratic')
      subplot(4,1,2),semilogy(costs(:,2)),title('Equality')
      subplot(4,1,3),semilogy(costs(:,3)),title('Lagrangian')
      subplot(4,1,4),semilogy(costs(:,1)+costs(:,3)),title('Quadratic + Lagrangian')
      drawnow
    end
  end
  eps = norm(s_reg(:)-s_reg_old(:))/norm(s_reg_old(:))
  s_reg_old = s_reg;
end

s_reg=projs(s_reg);
if debug_flag
  rho
  costs=[costs; [cost_quad(s_reg,G(s_reg)), cost_equality(s_reg) cost_lagrange(s_reg,y,projs(s_reg+y/rho),rho)]];
  if figure_flag
    figure(99)
    subplot(4,1,1),plot(costs(:,1)),title('Quadratic')
    subplot(4,1,2),semilogy(costs(:,2)),title('Equality')
    subplot(4,1,3),semilogy(costs(:,3)),title('Lagrangian')
    subplot(4,1,4),semilogy(costs(:,1)+costs(:,3)),title('Quadratic + Lagrangian')
    drawnow
  end
  
end

%Compare quadratic penalty portion of projected sense map and final estimated sense map
  [cost_quad(projs(s_reg_unc),G(projs(s_reg_unc))) cost_quad(s_reg,G(s_reg))]

if figure_flag
  figure,imshow3(abs(s_reg),[0 1]),title('Smooth, Normalized Sensitivity Maps (Magnitude)')
  figure,imshow3(angle(s_reg),[-pi pi]),title('Smooth, Normalized Sensitivity Maps (Phase)')
  
%   figure,imshow3(abs(s_reg-s_noreg),[0 1]),title('Change in Sensitivity Maps (Magnitude)')
%   figure,imshow3(angle(s_reg./s_noreg),[-pi pi]),title('Change in Sensitivity Maps (Phase)')
  
%   figure,imshow3(abs(s_reg-s_reg_proj),[0 1]),title('Change in Sensitivity Maps from Smooth Proj. (Magnitude)')
%   figure,imshow3(angle(s_reg./s_reg_proj),[-pi pi]),title('Change in Sensitivity Maps from Smooth Proj. (Phase)')
  
  im_c_reg=sum(ims.*conj(s_reg),4)./sos(s_reg);
  figure,subplot(1,2,1),imshow(abs(im_c_reg(:,:,center_z)),[0 max(abs(im_c(:)))]),title('Coil combination with smooth, normalized maps')
  subplot(1,2,2),imshow(angle(im_c_reg(:,:,center_z)),[-pi pi])
  figure,imshow(sqrt(sos(s_reg(:,:,center_z,:))),[0 1]),title('Relative SNR (Sensitivity map combination)')
  drawnow
end

return