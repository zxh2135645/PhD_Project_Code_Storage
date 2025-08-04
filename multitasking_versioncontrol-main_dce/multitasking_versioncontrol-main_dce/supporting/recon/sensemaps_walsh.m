function [s_walsh, im_c_walsh, s_noreg] = sensemaps_walsh(ims,window_length,varargin)
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
    window_length = 5;
elseif isempty(window_length)
    window_length = 5;
elseif ischar(window_length) %if window_length omitted but options included
    varargin{end+1} = window_length;
    window_length = 5;
end

if numel(window_length) == 1
    kernel = [window_length window_length min(window_length,size(ims,3))];
else
    window_length = [window_length 1];
    kernel = ones(1,3);
    kernel(1) = max(min(window_length(1),size(ims,1)),3);
    kernel(2) = max(min(window_length(2),size(ims,2)),3);
    kernel(3) = max(min(window_length(3),size(ims,3)),1);
end

if (numel(size(ims)) == 4) && (size(ims,3) > 1) %if 3D
    dims = 3;
    center_z = floor(size(ims,3)/2+1);
else    %if 2D
    ims  = reshape(ims,size(ims,1),size(ims,2),1,[]);
    dims = 2;
    center_z = 1;
end

[Np, Nf, Ns, coils] = size(ims);
sizes = size(ims);
vec   = @(x) x(:);
sos   = @(s) sum(abs(reshape(s,sizes(1),sizes(2),[],sizes(4))).^2,4);
prep  = @(x,sizes) reshape(x,sizes(1),sizes(2),sizes(3),[]);

if dims == 2
    imshow3 = @(x,caxis)montage(gather(x),'DisplayRange',caxis);
    fftim   = @(x)fft(fft(x,[],1),[],2)/sqrt(sizes(1))/sqrt(sizes(2));
    ifftim  = @(x)ifft(ifft(x,[],1),[],2)*sqrt(sizes(1))*sqrt(sizes(2));
    fftim3d   = @(x)ft1d(ft1d(x,1),2);
    ifftim3d  = @(x)ift1d(ift1d(x,1),2);
    fftshift3d  = @(x)fftshift(fftshift(x,1),2);
    ifftshift3d = @(x)ifftshift(ifftshift(x,1),2);
else %3D
    imshow3 = @(x,caxis)montage(gather(x(:,:,center_z,:)),'DisplayRange',caxis);
    fftim   = @(x)fft(fft(fft(x,[],1),[],2),[],3)/sqrt(sizes(1))/sqrt(sizes(2))/sqrt(sizes(3));
    ifftim  = @(x)ifft(ifft(ifft(x,[],1),[],2),[],3)*sqrt(sizes(1))*sqrt(sizes(2))*sqrt(sizes(3));
    fftim3d  = @(x)ft1d(ft1d(ft1d(x,1),2),3);
    ifftim3d = @(x)ift1d(ift1d(ift1d(x,1),2),3);
    fftshift3d  = @(x)fftshift(fftshift(fftshift(x,1),2),3);
    ifftshift3d = @(x)ifftshift(ifftshift(ifftshift(x,1),2),3);
end

% ismember doesn't work well in newer MATLAB when using mixed-class cell
% arrays, so handle potential acs mask and k-space windows first
kspace_loc = find(strcmp(varargin,'acs')) + 1;
acs_flag   = ~isempty(kspace_loc);
if acs_flag
    center_mask = varargin{kspace_loc};
    varargin(kspace_loc) = [];
end

window_loc  = find(strcmp(varargin,'window')) + 1;
window_flag = ~isempty(window_loc);
if window_flag
    window = varargin{window_loc};
    window = window/max(abs(window(:))); %normalize
    varargin(window_loc) = [];
end

figure_loc = find(strcmp(varargin,'figures')) + 1;
if isempty(figure_loc)
    figure_flag = false;
else
    figure_flag = varargin(figure_loc);
    figure_flag = figure_flag{1};
    varargin(figure_loc) = [];
end

voxelSpacing_loc  = find(strcmp(varargin,'voxelSpacing')) + 1;
if ~isempty(voxelSpacing_loc)
    voxelSpacing = varargin{voxelSpacing_loc};
    kernel(2) = min(floor(kernel(2)*voxelSpacing(1)/voxelSpacing(2)/2)*2 + 1,Nf);
    kernel(3) = min(floor(kernel(3)*voxelSpacing(1)/voxelSpacing(3)/2)*2 + 1,Ns);
    fprintf('Smoothing kernel size = [%d %d %d]. ',kernel(1),kernel(2),kernel(3));
    varargin(voxelSpacing_loc) = [];
end

kspace_flag = ismember('kspace',varargin);  % all lowercase
debug_flag  = ismember('debug',varargin);   % all lowercase

if debug_flag
    varargin
end
  
if kspace_flag
    kspace = ifftshift3d(ims);
    ims    = ifftim(kspace);
else
    ims = ifftshift3d(ims);
    kspace = fftim(ims);
end

temp = sos(kspace);
temp_shift = ifftshift3d(temp);
if temp_shift(1) > temp(1)  %if kspace isn't yet ifftshifted
    kspace = ifftshift3d(kspace);
    ims    = ifftim(kspace);
    if window_flag
        window = ifftshift(window);
    end
end
temp = sos(kspace);
temp_shift = ifftshift3d(temp);
if temp_shift(1) > temp(1)  %if kspace isn't yet ifftshifted
    kspace = ifftshift3d(kspace);
    ims    = ifftim(kspace);
end
if window_flag
    temp = window.^2;
    temp_shift = ifftshift3d(temp);
    if temp_shift(1) > temp(1)  %if kspace isn't yet ifftshifted
        window = ifftshift3d(window);
    end
end
if acs_flag
    temp = center_mask.^2;
    temp_shift = ifftshift3d(temp);
    if temp_shift(1) > temp(1)  %if kspace isn't yet ifftshifted
        center_mask = ifftshift3d(center_mask);
    end
end

clear temp temp_shift

if window_flag
    window = repmat(window,[1 1 1 coils]);
    center_mask = abs(window) > 0;
    center_data = window.*center_mask.*fftim(ims);
    window_comp = ones(size(window)); %no need to compensate for window;
    ims = ifftim(center_data); %use normalized fft here
elseif acs_flag
    % Identify ACS region...assumes cube or rectangle at center of k-space  
    rows = find(sum(sum(center_mask,2),3));
    if numel(rows) < Np
        rows = rows([1:find(diff(rows)>1,1) find(diff(rows)>1,1,'last')+1:end]);
    end
  
    cols = find(sum(sum(center_mask,1),3));
    if numel(cols) < Nf
        cols = cols([1:find(diff(cols)>1,1) find(diff(cols)>1,1,'last')+1:end]);
    end
  
    slices = find(sum(sum(center_mask,1),2));
    if numel(slices) < Ns
        slices = slices([1:find(diff(slices)>1,1) find(diff(slices)>1,1,'last')+1:end]);
    end
    
    center_mask = repmat(center_mask,[1 1 1 coils]);
    center_data = center_mask.*fftim(ims);

    if Ns > 5
        window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ...
                    .* repmat(reshape( hamming(numel(slices),'periodic'), 1,1,[]), [numel(rows), numel(cols), 1]);
    else
        window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ;
    end
    window = zeros(Np,Nf,Ns,coils);
    window(rows,cols,slices,:) = repmat(ifftshift(window_sm),[1 1 1 coils]);

    window_comp = 1./window; %window compensation
    window_comp(isinf(window_comp)) = 0;
    ims = ifftim(window.*center_data); % use normalized fft here
end

% Normalize ims to similar scale everytime. Does not affect scaling of
% sensitivity maps.
if acs_flag || window_flag
    center_data = center_data/norm(ims(:))*coils;
end
ims = ims/norm(ims(:))*coils;

ims = ifftshift3d(ims);

% Create coil combination image: sum-of-squares magnitude,
%                                first eigenimage phase
[U,~,mixer] = svd(reshape(ims,[],coils),0);
im_c = reshape(U(:,1),Np,Nf,Ns);
im_c = sqrt(sos(ims)).*exp(1i*angle(im_c));

if figure_flag
    figure,subplot(1,2,1),imshow(abs(im_c(:,:,center_z)),[]),title('SoS combination (Magnitude)');drawnow;
    subplot(1,2,2),imshow(angle(im_c(:,:,center_z)),[-pi pi]),title('1st eigenimage (Phase)');drawnow;
end

% Solve arg min_s || ims - A(s) ||
% noreg: no regularization
s_noreg = ims./repmat(im_c,[1 1 1 coils]); %closed-form
if figure_flag
    figure,imshow3(abs(s_noreg),[0 1]),title('Naive Sensitivity Maps (Magnitude)');drawnow;
    figure,imshow3(angle(s_noreg),[-pi pi]),title('Naive Sensitivity Maps (Phase)');colormap(hsv);drawnow;
end

% do image-space blockwise SVD
s_walsh = s_noreg;
lby  = floor(kernel(1)/2);
lbx  = floor(kernel(2)/2);
winy = @(x,xmax) mod(x+(-lby:lby)-1,xmax)+1;
winx = @(x,xmax) mod(x+(-lbx:lbx)-1,xmax)+1;
if dims == 2
    winz = @(x,xmax) 1;
else
    lbz  = floor(kernel(3)/2);
    winz = @(x,xmax) mod(x+(-lbz:lbz)-1,xmax)+1;
end
disp('SVDs...')
progress = waitbar(0,'','Name', sprintf('Progress'),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
setappdata(progress,'canceling',0)
waitbar(0, progress, sprintf('Estimating sensitivity maps...'));
for j = 1:Np
    for k = 1:Nf
        for l = 1:Ns
            [~,~,V] = svde(reshape(ims(winy(j,Np),winx(k,Nf),winz(l,Ns),:),[],coils));
            s_walsh(j,k,l,:) = conj(V(:,1))*exp(1i*angle(mean(V(:,1))./mean(mixer(:,1))));
        end
    end
    if getappdata(progress,'canceling')
        break
    else
        waitbar(j/Np, progress);
    end
end
delete(progress);
    
s_walsh = s_walsh./sqrt(sos(s_walsh));
ims = ifftshift3d(ifftim(kspace))/norm(ims(:))*coils;
im_c_walsh = sum(ims.*conj(s_walsh),4);
% if figure_flag
%   figure,imshow(sqrt(sos(s_walsh(:,:,center_z,:))),[]),caxis(caxis.*[0 1]),title('Relative SNR (Sensitivity map combination)'); drawnow;
%   figure,imshow3(abs(s_walsh-s_noreg),[0 1]),title('Change in Sensitivity Maps (Magnitude)')
%   figure,imshow3(angle(s_walsh./s_noreg),[-pi pi]),title('Change in Sensitivity Maps (Phase)')
    
%     figure,imshow3(abs(s_walsh),[0 1]),title('Smoothed sensitivity Maps (Magnitude)');drawnow;
%     figure,imshow3(angle(s_walsh),[-pi pi]),title('Smoothed sensitivity Maps (Phase)');colormap(hsv);drawnow;
    
%     figure,
%     subplot(1,2,1),imshow(abs(im_c_walsh(:,:,center_z)),[0 max(abs(im_c_walsh(:)))]);title('Coil combination with smooth maps');
%     subplot(1,2,2),imshow(angle(im_c_walsh(:,:,center_z)),[-pi pi]);drawnow;
% end

