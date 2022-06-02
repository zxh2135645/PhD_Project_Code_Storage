clear all; 
close all;
%% 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,9,:),1);

temp = Gr\reshape(Phi(:,161,1,1,:), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));


% ax1 = implay(abs(temp/cw));
figure();
ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
% for i = 1:size(temp,3)
%     subplot(1,3,i) 
%     ax2 = imagesc(abs(temp(:,:,i)/cw)); axis image; colormap gray;
% end
%% whole-stack
figure();
for i = 1:Nz
    dispim = @(x)fftshift(x(:,:,i,:),1);
    
    temp = Gr\reshape(Phi(:,161,1,end,:), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    cw = max(vec(abs(temp)));
    
    
    % ax1 = implay(abs(temp/cw));
    subplot(4,4,i);
    ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
end

%% Pull up images - draw heart masks
addpath('../function/');
soi = [6, 5, 4, 3, 2, 1, 14, 13, 12, 11];
% soi = [5, 4, 3, 2, 1, 16, 15, 14, 13, 12];
t1oi = [21, 41, 61, 81, 101, 121, 141, 161, 181];
wx = 64;
wy = 64;
mask = zeros(Ny, Nx, Nz);
figure('Position',[100, 100, 720, 800]);
% tight_subplot(Nh, Nw, gap, marg_h, marg_w)
for i = 1:length(soi)
    slc = soi(i);
    dispim = @(x) fftshift(x(:,:,slc,:), 1);
    for j = 2:2
        % Seg#41
        t1_idx = t1oi(j);
        temp = Gr\reshape(Phi(:,t1_idx,1,end,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
        cw = max(vec(abs(temp)));
        imagesc(abs(temp(:,:,1)/cw)); axis image; 
        colormap gray; axis off;
        
        roi = drawpolygon;
        mask(:,:,slc) = createMask(roi);
    end
end

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    slc = soi(i);
    subplot(4,4,i)
    imagesc(mask(:,:,slc))
    s = regionprops(mask(:,:,slc),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

%% pull-up - T1 weighted
% tight_subplot(Nh, Nw, gap, marg_h, marg_w)
Nseg = size(Phi,2);
Nsect = size(Phi,5);
img_save = cat(2, fid_path, 'img/');
if ~exist(img_save, 'dir')
   mkdir(img_save); 
end
for nsec = 1:Nsect
    %nsec = 1;
    figure('Position',[100, 100, 720, 800]);
    [ha, pos] = tight_subplot(length(soi),length(t1oi),[.01 -.1],[.01 .01],[.01 .01]);
    for i = 1:length(soi)
        slc = soi(i);
        dispim = @(x) fftshift(x(:,:,slc,:), 1);
        temp = Gr\reshape(Phi(:,:,1,1,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
        %cw = 0.5*max(vec(abs(temp)));
        centroid = centroids{i};
        
        for j = 1:length(t1oi)
            t1_idx = t1oi(j);
            temp_crop = imcrop(abs(temp(:,:,t1_idx,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
            cw = 0.5*max(vec(temp_crop));
            
            axes(ha((i-1)*length(t1oi)+j));
            imagesc(temp_crop/cw); axis image;
            colormap gray; axis off;
        end
    end
    fname = cat(2, fid_file(1:17), fid_file(end-15:end-6), num2str(nsec), '.png');
    saveas(gcf,cat(2, img_save, fname));
    close all;
end

%% Pull-up T2* weighted - single slice
Nseg = size(Phi,2);
Nsect = size(Phi,5);
img2_save = cat(2, fid_path, 'img_t2star/');
if ~exist(img2_save, 'dir')
   mkdir(img2_save); 
end
NEco_old = params.NEco_old;
echo_f_glob = glob(cat(2, fid_path, '*Echo*.mat'));

soi = [6, 5, 4, 3, 2, 1, 14, 13, 12, 11];
slc = 3;
idx = find(soi == slc);
dispim = @(x) fftshift(x(:,:,slc,:), 1);
centroid = centroids{idx};

SpecialNote = '_HalfLambda';
for nsec = 1:Nsect
    %nsec = 1;
    for neco = 1:NEco_old
       if neco == 1
          temp_4D = zeros([Ny, Nx, Nseg, NEco_old]); 
       end
       load(echo_f_glob{neco}, 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nz', 'params', 'sizes');
       temp = Gr\reshape(Phi(:,:,1,1,nsec), L, []);
       temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, []);
       % cw = 0.5*max(vec(abs(temp)));
       temp_4D(:,:,:,neco) = temp;
    end
    
    figure('Position',[100, 100, 900, 600]);
    [ha, pos] = tight_subplot(NEco_old,length(t1oi),[.01 -.1],[.01 .01],[.01 .01]);
    
    for i = 1:NEco_old
        temp = temp_4D(:,:,:,i);
        for j = 1:length(t1oi)
            t1_idx = t1oi(j);
            temp_crop = imcrop(abs(temp(:,:,t1_idx)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
            cw = 0.5*max(vec(temp_crop));
            
            axes(ha((i-1)*length(t1oi)+j));
            imagesc(temp_crop/cw); axis image;
            colormap gray; axis off;
        end
    end
    fname = cat(2, fid_file(1:17), fid_file(end-15:end-6), num2str(nsec), '_slc', num2str(slc), SpecialNote, '.png');
    saveas(gcf,cat(2, img2_save, fname));
    close all;
end

%% pull-up - T1 weighted - PSIR recon
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');

Nseg = size(Phi,2);
Nsect = size(Phi,5);
img_save = cat(2, fid_path, 'img_psir/');
if ~exist(img_save, 'dir')
   mkdir(img_save); 
end

for nsec = 1:Nsect
    %nsec = 1;
    figure('Position',[100, 100, 720, 800]);
    [ha, pos] = tight_subplot(length(soi),length(t1oi),[.01 -.1],[.01 .01],[.01 .01]);
    for i = 1:length(soi)
        slc = soi(i);
        dispim = @(x) fftshift(x(:,:,slc,:), 1);
        temp = Gr\reshape(Phi(:,:,1,end,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
        
        phase_temp = angle(temp);
        phase_diff = phase_temp - phase_temp(:,:,end);
        cos_phase_diff = cos(phase_diff);
        %cos_phase_diff(cos_phase_diff >= 0) = 1;
        %cos_phase_diff(cos_phase_diff <  0) = -1;
        
        %cw = 0.5*max(vec(abs(temp)));
        centroid = centroids{i};
        
        for j = 1:length(t1oi)
            t1_idx = t1oi(j);
            % temp_crop = imcrop(abs(temp(:,:,t1_idx,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
            temp_psir = abs(temp(:,:,t1_idx,nsec)) .* cos_phase_diff(:,:,t1_idx,nsec);
            psir = (temp_psir - min(temp_psir(:))) ./ (max(temp_psir(:)) - min(temp_psir(:)));
            psir_crop = imcrop(psir, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
            % cw = 0.5*max(vec(psir_crop));

            axes(ha((i-1)*length(t1oi)+j));
            imagesc(psir_crop); axis image;
            colormap gray; axis off;
        end
    end
    fname = cat(2, fid_file(1:17), fid_file(end-15:end-6), num2str(nsec), '.png');
    saveas(gcf,cat(2, img_save, fname));
    close all;
end

% pull-up - montage of mag vs psir
Nseg = size(Phi,2);
Nsect = size(Phi,5);
img_save = cat(2, fid_path, 'img_montage/');
if ~exist(img_save, 'dir')
   mkdir(img_save); 
end

odd_array = (1:Nseg/2) * 2 - 1;
caxis_cell = {[0.3 0.9],[0.2 0.75]};
Nsect_array = [1, 15];
for n = 1:length(Nsect_array)
    nsec = Nsect_array(n);
    figure('Position',[0, 100, 1800, 800]);
    [ha, pos] = tight_subplot(1,3,[.01 -.1],[.01 .01],[.01 .01]);
    
        for slc = 3:3
            % slc = soi(i);
            dispim = @(x) fftshift(x(:,:,slc,:), 1);
            temp = Gr\reshape(Phi(:,:,1,end,:), L, []);
            temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
            
            phase_temp = angle(temp);
            phase_diff = phase_temp - phase_temp(:,:,end);
            cos_phase_diff = cos(phase_diff);
            cos_phase_diff(cos_phase_diff >= 0) = 1;
            cos_phase_diff(cos_phase_diff <  0) = -1;
            
            idx = find(soi == slc);
            centroid = centroids{idx};
            
            temp_crop = imcrop3(abs(temp(:,:,odd_array,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            temp_psir = abs(temp(:,:,odd_array,nsec)) .* cos_phase_diff(:,:,odd_array,nsec);
            psir = (temp_psir - min(temp_psir(:))) ./ (max(temp_psir(:)) - min(temp_psir(:)));
            psir_crop = imcrop3(psir, [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            % cos_phase_diff_crop = imcrop3(cos_phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            phase_diff_crop = imcrop3(phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
            cw1 = 0.5*max(vec(abs(temp(:,:,odd_array,nsec))));
            cw2 = 0.5*max(vec(phase_diff(:,:,odd_array,nsec)));
            %cw3 = 0.5*max(vec(psir));
            
            axes(ha(1));
            montage(temp_crop./cw1, 'Size', [12 8]);
            axes(ha(2));
            montage(phase_diff_crop./cw2, 'Size', [12 8]);
            caxis([-pi pi]);
            axes(ha(3));
            montage(psir_crop, 'Size', [12 8]);
            caxis(caxis_cell{n});
        end
    fname = cat(2, fid_file(1:17), fid_file(end-15:end-6), num2str(nsec), '.png');
    saveas(gcf,cat(2, img_save, fname));
    close all;
end