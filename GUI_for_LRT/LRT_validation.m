clear all;
close all;

base_dir = uigetdir;

load(cat(2, base_dir, '/Results/FID26795_21D05_LRT_Mappings.mat'));

t2star_map = squeeze(map_to_save.t2star_map(:,:,:,1,1,end));
reorder_slice = [8,9,10,11,12,1,2,3,4,5,6,7];
t2star_map = t2star_map(:,:,reorder_slice);

figure();
for i = 1:size(t2star_map, 3)
    subplot(3,4,i);
    imagesc(t2star_map(:,:,i)); axis image; caxis([0 100]);
    
end

% t1_map doesn't look right
t1_map = squeeze(map_to_save.t1_map(:,:,:,1,1,40));
reorder_slice = [8,9,10,11,12,1,2,3,4,5,6,7];
t1_map = t1_map(:,:,reorder_slice);
figure();
for i = 1:size(t1_map, 3)
    subplot(3,4,i);
    imagesc(t1_map(:,:,i)); axis image; caxis([0 500]); 
end

%% Draw Contours

roi_save = cat(2, base_dir, '/roi.mat');

if ~exist(roi_save, 'file')
    for i = 1:(size(t2star_map, 3))
        disp('Epicardium: ');
        figure('Position', [100 0 1600 1600]); imagesc(t2star_map(:,:,i)); axis image;
        caxis([0 100])
        epi = drawpolygon(gca);
        epi_coords = epi.Position;
        
        disp('Endocardium: ');
        figure('Position', [100 0 1600 1600]); imagesc(t2star_map(:,:,i)); axis image;
        caxis([0 100])
        endo = drawpolygon(gca);
        endo_coords = endo.Position;
        
        
        myo_coords_cell{i, 1} = epi.Position;
        myo_coords_cell{i, 2} = endo.Position;
        epi_mask = createMask(epi);
        endo_mask = createMask(endo);
        
        close all;
    end
    
    roi.myo_coords_cell = myo_coords_cell;
    save(roi_save, 'roi');
else
    load(roi_save);
    myo_coords_cell = roi.myo_coords_cell;
end

mask_save = cat(2, base_dir, '/mask.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    epi_mask = zeros(size(t2star_map));
    endo_mask = zeros(size(t2star_map));
    for i = 1:size(t2star_map, 3)
        imagesc(t2star_map(:,:,i)); caxis([0 100]);
        epi = drawpolygon(gca,'Position', [myo_coords_cell{i,1}(:,1), myo_coords_cell{i,1}(:,2)]);
        endo = drawpolygon(gca,'Position', [myo_coords_cell{i,2}(:,1), myo_coords_cell{i,2}(:,2)]);
        epi_mask(:,:,i) = createMask(epi);
        endo_mask(:,:,i) = createMask(endo);
    end
    
    myo_mask = epi_mask - endo_mask;
    mask_struct.epi_mask = epi_mask;
    mask_struct.endo_mask = endo_mask;
    mask_struct.myo_mask = myo_mask;
    
    save(mask_save, 'mask_struct');
    
else
    load(mask_save);
    epi_mask = mask_struct.epi_mask;
    endo_mask = mask_struct.endo_mask;
    myo_mask = mask_struct.myo_mask;
end

%% Draw insertion point
ref_pts_save = cat(2, base_dir, '/ref_pts.mat');
if ~exist(ref_pts_save, 'file')
   figure();
   x = zeros(size(t2star_map, 3),1);
   y = zeros(size(t2star_map, 3),1);
   for i = 1:size(t2star_map, 3)
      imagesc(t2star_map(:,:,i)); caxis([0 100]);
      [x(i),y(i)] = getpts(gca); 
   end
   ref_pts.x = x;
   ref_pts.y = y;
   save(ref_pts_save, 'ref_pts');
   
else
    load(ref_pts_save);
    x = ref_pts.x;
    y = ref_pts.y;
end

%% AHA 16 segment
BaseGroove = zeros(size(t2star_map, 3), 1);
for i = 1:size(t2star_map, 3)
    C = regionprops(myo_mask(:,:,i));
    x_centroid = C.Centroid(2);
    y_centroid = C.Centroid(1);
    BaseGroove(i) = atan2(x(i) - x_centroid, y(i) - y_centroid) * 180 / pi;
end

% tease out margin cases
m = 1;
BaseGroove_truc = BaseGroove((1+m):(end-m));
t2star_map_truc = t2star_map(:,:,(1+m):(end-m));
edg = zeros(size(t2star_map_truc));
myo_mask_truc = myo_mask(:,:,(1+m):(end-m));
for i = 1:size(t2star_map_truc, 3)
    edg(:,:,i)  = edge(squeeze(myo_mask_truc(:,:,i)),'Canny');
end

myo_mask_peel = myo_mask_truc - edg;


n = size(t2star_map_truc, 3);
mode = mod(n,3);
integ = fix(n/3);
if n >= 3
    switch mode
        case {0}
            aha_slice = cat(2, repmat([1], [1, integ]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {1}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {2}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ+1]), repmat([3], [1, integ]));
    end
else
    error("Available slice numbers are smaller than 3.");
end
%%
addpath('./function/');
% Basal
LocPixCount1 = zeros(6, 1);
SegTotalPixCount1 = zeros(6, 1);
idx_array = 1:size(t2star_map_truc,3);

basal_idx = idx_array(aha_slice == 1);
Groove = BaseGroove(basal_idx(1)) + 60;


[Segmentpix, stats, Mask_Segn] = AHASegmentation(t2star_map_truc(:,:,basal_idx), myo_mask_peel(:,:,basal_idx), 6, Groove);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount1(i) = LocPixCount1(i) + sum(Segmentpix{i,j});
        SegTotalPixCount1(i) = SegTotalPixCount1(i) + length(Segmentpix{i,j});
    end
end

% Mid-ventricular
LocPixCount2 = zeros(6, 1);
SegTotalPixCount2 = zeros(6, 1);

mid_idx = idx_array(aha_slice == 2);
Groove = BaseGroove(mid_idx(1)) + 60;

[Segmentpix, stats, Mask_Segn] = AHASegmentation(t2star_map_truc(:,:,mid_idx), myo_mask_peel(:,:,mid_idx), 6, Groove);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount2(i) = LocPixCount2(i) + sum(Segmentpix{i,j});
        SegTotalPixCount2(i) = SegTotalPixCount2(i) + length(Segmentpix{i,j});
    end
end

% Apical
LocPixCount3 = zeros(4, 1);
SegTotalPixCount3 = zeros(4, 1);

apical_idx = idx_array(aha_slice == 2);
Groove = BaseGroove(apical_idx(1)) + 75;

[Segmentpix, stats, Mask_Segn] = AHASegmentation(t2star_map_truc(:,:,apical_idx), myo_mask_peel(:,:,apical_idx), 4, Groove);

for i = 1:4
    for j = 1:size(Segmentpix, 2)
        LocPixCount3(i) = LocPixCount3(i) + sum(Segmentpix{i,j});
        SegTotalPixCount3(i) = SegTotalPixCount3(i) + length(Segmentpix{i,j});
    end
end

SegPixCount = [LocPixCount1; LocPixCount2; LocPixCount3];
SegTotalPixCount = [SegTotalPixCount1; SegTotalPixCount2; SegTotalPixCount3];
SegPixPerc = SegPixCount ./ SegTotalPixCount;


% =================== Create the AHA 17-segment bullseye =================
% ========================================================================
figure('Position', [400 400 800 800]);
PlotBullsEye(SegPixPerc);

