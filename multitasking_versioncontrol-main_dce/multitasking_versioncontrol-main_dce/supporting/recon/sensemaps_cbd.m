function [s_cbd, im_c_cbd, ims, s_noreg] = sensemaps_cbd(ims,kernel_size,varargin)
% Estimate smooth, normalized sensitivity maps (constrained blind deconvolution)
% v0.5 by Anthony Christodoulou, 15 December 2016 - 18 April 2019
%
% sensemaps_sn(ims,kernel_size)
% sensemaps_sn(kspace,kernel_size,'kspace')
%
% ims or kspace:    Phase-encode X freq-encode X coils
%                or
%                   Phase-encode X freq-encode X partition-encode X coils
%
% kernel_size:      Convolution kernel size (default is [7 7] or [7 7 3])
%
% Options:
%
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

if nargin < 2
    kernel_size = [7 7 5];
elseif isempty(kernel_size)
    kernel_size = [7 7 5];
elseif ischar(kernel_size) %if kernel_size omitted but options included
    varargin = [{kernel_size}, varargin];
    kernel_size = [7 7 5];
end
if size(ims,3) < 5
    kernel_size(3) = size(ims,3);
end

try
    %ims = gpuArray(ims);
    
    [s_cbd, im_c_cbd, s_noreg] = sensemaps_cbd_calculation(ims,kernel_size,varargin);
    
    s_cbd    = gather(s_cbd);
    im_c_cbd = gather(im_c_cbd);
    s_noreg  = gather(s_noreg);
catch msg
    msg
    disp('Trying again on CPU');
    ims = gather(ims);
    [s_cbd, im_c_cbd, s_noreg] = sensemaps_cbd_calculation(ims,kernel_size,varargin);
end

function [s_cbd, im_c_cbd, s_noreg] = sensemaps_cbd_calculation(ims,kernel_size,varargin)

varargin = varargin{1}; %varargin is passed inside a cell

if (numel(size(ims)) == 4) && (size(ims,3) > 1) % if 3D
    dims = 3;
    center_z = floor(size(ims,3)/2 + 1);
else    % if 2D
    ims  = reshape(ims,size(ims,1),size(ims,2),1,[]);
    dims = 2;
    center_z = 1;
end

[Np, Nf, Ns, coils] = size(ims);
sizes = [Np, Nf, Ns, coils];
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
else % 3D
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
    kspace_mask = varargin{kspace_loc};
    varargin(kspace_loc) = [];
end

figure_loc = find(strcmp(varargin,'figures')) + 1;
if isempty(figure_loc)
    figure_flag = false;
else
    figure_flag = varargin(figure_loc);
    figure_flag = figure_flag{1};
    varargin(figure_loc) = [];
end

kspace_flag = ismember('kspace',varargin);  %all lowercase
debug_flag  = ismember('debug',varargin);   %all lowercase
realfirst_flag = ismember('realfirst',varargin);  %make first channel's map real during optimization. Experimental!

if debug_flag
    kernel_size
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
end
if acs_flag
    temp = kspace_mask.^2;
    temp_shift = ifftshift3d(temp);
    if temp_shift(1) > temp(1)  %if kspace isn't yet ifftshifted
        kspace_mask = ifftshift3d(kspace_mask);
    end
end

clear temp temp_shift

if acs_flag
    % Identify ACS region...assumes cube or rectangle at center of k-space
    rows = find(sum(sum(kspace_mask,2),3));
    if numel(rows) < Np
        rows = rows([1:find(diff(rows)>1,1) find(diff(rows)>1,1,'last')+1:end]);
    end
  
    cols = find(sum(sum(kspace_mask,1),3));
    if numel(cols) < Nf
        cols = cols([1:find(diff(cols)>1,1) find(diff(cols)>1,1,'last')+1:end]);
    end
  
    slices = find(sum(sum(kspace_mask,1),2));
    if numel(slices) < Ns
        slices = slices([1:find(diff(slices)>1,1) find(diff(slices)>1,1,'last')+1:end]);
    end
  
    % apply fbp_window
    kspace_mask = repmat(kspace_mask,[1 1 1 coils]);
    kspace_orig = fftim(ims);
    kspace      = kspace_mask .* kspace_orig;
  
    % apply acs mask with hamming window
    if Ns > 5
        window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ...
                    .* repmat(reshape( hamming(numel(slices),'periodic'), 1,1,[]), [numel(rows), numel(cols), 1]);
    else
        window_sm = repmat( hamming(numel(rows),'periodic') * hamming(numel(cols),'periodic').', [1 1 numel(slices) ] ) ;
    end
    window = zeros(Np,Nf,Ns,coils);
    window(rows,cols,slices,:) = repmat(ifftshift(window_sm),[1 1 1 coils]);
  
    % low-resolution image with only acs lines
    ims = ifftim(window.*kspace_orig);
else
    kspace_mask = ones(size(kspace));
end

% ceate coil combination image: sum-of-squares magnitude & first eigenimage phase
[~,~,V]  = svd(reshape(ims,[],coils),0);
im_c     = sqrt(sos(ims));%.*sign(reshape(reshape(ims,[],coils)*V(:,1),Np,Nf,Ns));
im_c_cbd = im_c;
s_noreg  = bsxfun(@rdivide, ims, im_c);
s_cbd    = s_noreg; 
s_noreg  = ifftshift3d(s_noreg);

if figure_flag
    if acs_flag
        im_orig = ifftim(kspace_orig);
        clear kspace_orig
        im_orig = sqrt(sos(im_orig)).*sign(reshape(reshape(im_orig,[],coils)*V(:,1),Np,Nf,Ns));
    else
        im_orig = im_c;
    end
    im_orig = ifftshift3d(im_orig);
    figure;
    subplot(1,2,1);imshow(abs(im_orig(:,:,center_z)),[]);title('SoS combination (Magnitude)');
    subplot(1,2,2);imshow(angle(im_orig(:,:,center_z)),[-pi pi]);title('1st eigenimage (Phase)');
    figure;imshow3(abs(s_noreg),[0 1]);title('Initial Sensitivity Maps (Magnitude)');drawnow;
    figure;imshow3(angle(s_noreg),[-pi pi]);title('Initial Sensitivity Maps (Phase)');colormap(hsv);drawnow;
    clear im_orig
end

fprintf('Constrained blind deconvolution...');
hmask  = zeros([Np Nf Ns coils]);
coords = floor(kernel_size/2);

if dims == 2
%      hmask([Np-coords(1)+1:Np, 1:coords(1)+1],[Nf-coords(2)+1:Nf, 1:coords(2)+1],1,:) = 1;
    hmask(1:kernel_size(1),1:kernel_size(2),1,:) = repmat(hanning(kernel_size(1))*hanning(kernel_size(2))',[1 1 1 coils]);
    hmask = circshift(hmask,-coords);
    hmask = (sqrt(sqrt(hmask)));
    its = 10;
else
%    hmask([Np-coords(1)+1:Np, 1:coords(1)+1],[Nf-coords(2)+1:Nf, 1:coords(2)+1],[Ns-coords(3)+1:Ns, 1:coords(3)+1],:) = 1;
    temp = reshape(hanning(kernel_size(3)),1,1,[]);
    hmask(1:kernel_size(1),1:kernel_size(2),1:kernel_size(3),:) = repmat(bsxfun(@times,hanning(kernel_size(1))*hanning(kernel_size(2))', temp),[1 1 1 coils]);
    hmask = circshift(hmask,-coords);
    its = 5;
end
% hmask = sqrt(sqrt(sqrt(hmask)));
hmask = sqrt(hmask);

% try 
%     [~,bias] = N4BiasCorrection(abs(im_c_cbd(:,:,center_z)));
% catch
%     hmask = sqrt(sqrt(sqrt(hmask)));
% end

s_cbd = ifftim(fftim(s_cbd).*hmask);
%s_cbd = bsxfun(@rdivide,s_cbd,sqrt(sos(s_cbd)));
if realfirst_flag
    s_cbd = bsxfun(@rdivide,s_cbd,sign(s_cbd(:,:,:,1))); % make first map real
end

% Progress bar
progress = waitbar(0,'','Name', sprintf('Progress'),...
                   'CreateCancelBtn',...
                   'setappdata(gcbf,''canceling'',1)');
setappdata(progress,'canceling',0)
waitbar(0, progress, sprintf('Estimating sensitivity maps...(1/%d)',its));
    
try
    for it = 1:its
        
        % update coil combination image
        if acs_flag
            Ah  = @(x)vec(sum(bsxfun(@times,conj(s_cbd),ifftim(kspace_mask.*prep(x,sizes))),4));
            AhA = @(x)vec(sum(bsxfun(@times,conj(s_cbd),ifftim(kspace_mask.^2.*fftim(bsxfun(@times,prep(x,sizes),s_cbd)))),4));
            [temp,flag] = pcg(AhA,Ah(kspace),[],10,[],[],im_c_cbd(:));
            im_c_cbd = reshape(temp,Np,Nf,Ns);
        else
            im_c_cbd = sum(ims.*conj(s_cbd),4)./sos(s_cbd);
        end
            
        % update sensitivity maps
        if acs_flag
            adjh    = @(x)vec(fftim(bsxfun(@times,conj(im_c_cbd),ifftim(kspace_mask.*prep(x,sizes).*hmask))).*hmask);
            adjfwdh = @(h)vec(fftim(bsxfun(@times,conj(im_c_cbd),ifftim(kspace_mask.^2.*fftim(bsxfun(@times,im_c_cbd,ifftim(prep(h,sizes)))).*hmask))).*hmask);
            adjhb   = adjh(kspace);
        else
            adjh    = @(x)vec(fftim(bsxfun(@times,conj(im_c_cbd),prep(x,sizes))).*hmask);
            adjfwdh = @(h)vec(fftim(bsxfun(@times,abs(im_c_cbd).^2,ifftim(prep(h,sizes).*hmask))).*hmask);
            adjhb   = adjh(ims);
        end
        
        h = fftim(s_cbd);
        [temp,flag] = pcg(@(x)adjfwdh(x),adjhb,[],10,[],[],h(:));
        h = reshape(temp,Np,Nf,Ns,coils);
        s_cbd = ifftim(h.*hmask);        
        %s_cbd = bsxfun(@rdivide,s_cbd,sqrt(sos(s_cbd)));
        if realfirst_flag
            s_cbd = bsxfun(@rdivide,s_cbd,sign(s_cbd(:,:,:,1))); % make first map real
        end
        s_cbd(~isfinite(s_cbd)) = 0;

        % update figure
        if figure_flag && debug_flag
            try
                figure(fig(1))
            catch
                fig(1) = figure;
            end
            subplot(1,2,1),imshow3(abs(fftshift3d(s_cbd)),[0 1]),title('CBD Sensitivity Maps (Magnitude)')
            subplot(1,2,2),imshow3(angle(fftshift3d(s_cbd)),[-pi pi]),title('CBD Sensitivity Maps (Phase)')

            try
                figure(fig(3))
            catch
                fig(3) = figure;
            end
            subplot(1,2,1),imshow(abs(fftshift3d(im_c_cbd(:,:,1))),[0 max(abs(im_c(:)))]),title('Cartesian SENSE with CBD maps');
            subplot(1,2,2),imshow(angle(fftshift3d(im_c_cbd(:,:,1))),[-pi pi]);title('Phase');
            drawnow;
        end

        % update progress bar
        drawnow;
        if getappdata(progress,'canceling')
            break
        elseif it < its
            waitbar(it/its, progress, sprintf('Estimating sensitivity maps... (%d/%d)',it+1,its));
        end
        
        fprintf('.');
    end
    fprintf(' done.\n');
catch errormsg
    fprintf(2,'Sensitivity map estimation failed.\n');
    fprintf(2,'%s\n', errormsg.message);
end
delete(progress);

im_c_cbd = sum(ims.*conj(s_cbd),4)./sos(s_cbd);

if realfirst_flag
    phase_to_restore = sign(reshape(reshape(s_cbd,[],coils)*V(:,1),[Np Nf Ns]));
    im_c_cbd = im_c_cbd.*phase_to_restore;
    s_cbd    = bsxfun(@rdivide,s_cbd,phase_to_restore);
end

s_cbd    = ifftshift3d(s_cbd);
im_c_cbd = ifftshift3d(im_c_cbd);
ims      = ifftshift3d(ims);
if figure_flag
    figure;imshow3(abs(s_cbd),[0 1]);title('CBD Sensitivity Maps (Magnitude)');
    figure;imshow3(angle(s_cbd),[-pi pi]);title('CBD Sensitivity Maps (Phase)');
    figure;
    subplot(1,2,1),imshow(abs(im_c_cbd(:,:,center_z)),[0 max(abs(im_c(:)))]),title('SENSE with CBD maps (Magnitude)');
    subplot(1,2,2),imshow(angle(im_c_cbd(:,:,center_z)),[-pi pi]),title('SENSE with CBD maps (Phase)');
end

return
