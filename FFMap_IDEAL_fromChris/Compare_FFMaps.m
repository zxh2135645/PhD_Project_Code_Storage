clear all;
close all;

addpath('../function/');
Scan_Type='Cardiac';

[HBfile, HBpath] = uigetfile('*.mat');
load(cat(2, HBpath, HBfile));

name = 'Sahara';
time_point = '1YR'
base_dir = uigetdir;
ff_map = cell(1, length(glob_names));
for f = 1:length(ff_map)
    ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
end

% convert ff_map to matrix
ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
for f = 1:length(ff_map)
    ff(:,:,f) = ff_map{f}.fwmc_ff;
end


%% Comparison
ff_hb = abs(fat) ./ (abs(fat) + abs(water));
slc = 1;

figure();
for slc = 1:size(fat, 3)
    subplot(1,2,1);
    imagesc(ff(:,:,slc)); caxis([0 20]); colormap jet; colorbar; axis image;
    subplot(1,2,2);
    imagesc(abs(ff_hb(:,:,slc))*100); caxis([0 20]); colormap jet; colorbar; axis image;
    pause;
end

%%
figure();
for slc = 1:size(fat, 3)
    subplot(1,2,1);
    imagesc(ff(:,:,slc).*AllPhasemap(slc).Mask); caxis([0 20]); colormap jet; colorbar; axis image;
    subplot(1,2,2);
    imagesc(abs(ff_hb(:,:,slc))*100.*AllPhasemap(slc).Mask); caxis([0 20]); colormap jet; colorbar; axis image;
    pause;
end