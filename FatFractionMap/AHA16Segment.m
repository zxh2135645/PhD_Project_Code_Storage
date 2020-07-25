clear all;
close all;
clc
current_dir = pwd;
%% 
addpath('..\function\');
addpath('..\AHA16Segment\');

base_dir = GetFullPath(cat(2, current_dir, '\..\..\Data\Diane\ContourData\'));
img_dir = GetFullPath(cat(2, base_dir, '../../Diane/ResultsFolder_180718/'));

sequence_label = {'MAG', 'T1Map', 'MultiEcho'};
output_label = {'LGE', 'T1', 'MultiEcho'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'freeROI'};
overwrite_label = 1;

name_glob = glob(cat(2, base_dir, '/*_*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end

ff_glob = glob(cat(2, img_dir, name, '\', 'T2_multiecho_sax*_bright*.mat'));
idx_array = [4 5 6 7 8 9 10 1];

for name_idx = 1:length(Names)
    name = Names{name_idx};
    
    CoordsFileName_T1 = [base_dir, name, '\', sequence_label{2}, '_coords.mat'];
    T1_coords = load(CoordsFileName_T1);
    
    CoordsFileName_ME = [base_dir, name, '\', sequence_label{3}, '_coords.mat'];
    ME_coords = load(CoordsFileName_ME);
    
    mat_glob = glob([base_dir, name, '\', output_label{2}, '\', sequence_label{2}, '_vol_img_*D.mat']);
    img = load(mat_glob{1});
    fdnames = fieldnames(img);
    img_3D = getfield(img, fdnames{1});
    img_size = size(img_3D);
    
    myo_glob = glob(cat(2, base_dir, name, '/', output_label{2}, '/Myocardium/mask*.mat'));
    myo = load(myo_glob{1});
    fdnames = fieldnames(myo);
    myo_3D = getfield(myo, fdnames{1});
    
    mi_glob = glob(cat(2, base_dir, name, '/', output_label{2}, '/freeROI/freeROI.mat'));
    mi = load(mi_glob{1});
    fdnames = fieldnames(mi);
    mi_3D = getfield(mi, fdnames{1});
    
    Segmentpix_T1 = cell(img_size(3),1);
    Segmentpix_ff = cell(img_size(3),1);
    Mask_Segn_t1 = zeros(img_size);
    
    n = length(idx_array);
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
    % But how do I know where is apex, where is base
    
    x = T1_coords.coords.x;
    y = T1_coords.coords.y;
    x_centroid = T1_coords.coords.x_centroid;
    y_centroid = T1_coords.coords.y_centroid;
    BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
    
    figure();
    for slc = 1:img_size(3)
        aha = aha_slice(slc);
        switch aha
            case {1}
                Groove = BaseGroove(slc) + 75;
                Segn = 4;
            case {2}
                Groove = BaseGroove(slc) + 60;
                Segn = 6;
            case {3}
                Groove = BaseGroove(slc) + 60;
                Segn = 6;
        end
        
        
        [Segmentpix_T1{slc}, stats, Mask_Segn_t1(:,:,slc)] = AHASegmentation(img_3D(:,:,slc), myo_3D(:,:,slc), Segn, Groove);
        sq = ceil(sqrt(img_size(3)));
        subplot(sq, sq, slc)
        img_2D = img_3D(:,:,slc);
        img_2D_scale = img_2D / max(img_2D(:));
        imagesc(Mask_Segn_t1(:,:,slc) + 2*img_2D_scale); axis image;
    end
    
    % For MultiEcho FF map  
    x = ME_coords.coords.x;
    y = ME_coords.coords.y;
    x_centroid = ME_coords.coords.x_centroid;
    y_centroid = ME_coords.coords.y_centroid;
    BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
    
    
    myo_ff_glob = glob(cat(2, base_dir, name, '/', output_label{3}, '/Myocardium/mask*.mat'));
    myo_ff = load(myo_ff_glob{1});
    fdnames = fieldnames(myo_ff);
    myo_ff_4D = getfield(myo_ff, fdnames{1});
    myo_ff_3D = squeeze(myo_ff_4D(:,:,1,:));
    
    mi_ff_glob = glob(cat(2, base_dir, name, '/', output_label{3}, '/freeROI/freeROI.mat'));
    mi_ff = load(mi_ff_glob{1});
    fdnames = fieldnames(mi_ff);
    mi_ff_4D = getfield(mi_ff, fdnames{1});
    mi_ff_3D = squeeze(mi_ff_4D(:,:,1,:));
    
    Mask_Segn_ff = zeros(size(myo_ff_3D));
    figure();
    for slc = 1:img_size(3)
        aha = aha_slice(slc);
        ff = load(ff_glob{idx_array(slc)});
        ff_map = ff.fwmc_ff;

        switch aha
            case {1}
                Groove = BaseGroove(slc) + 75;
                Segn = 4;
            case {2}
                Groove = BaseGroove(slc) + 60;
                Segn = 6;
            case {3}
                Groove = BaseGroove(slc) + 60;
                Segn = 6;
        end
        [Segmentpix_ff{slc}, stats, Mask_Segn_ff(:,:,slc)] = AHASegmentation(ff_map, myo_ff_3D(:,:,slc), Segn, Groove);
        sq = ceil(sqrt(img_size(3)));
        subplot(sq, sq, slc)
        ff_2D_scale = ff_map / 100;
        imagesc(Mask_Segn_ff(:,:,slc) + 2*ff_2D_scale); axis image; caxis([0 6])
    end
    
    
end

%% 
clear seg_mean_t1 seg_mean_ff
count_t1 = 1;
count_ff = 1;
for slc = 1:(img_size(3)-4)
    seg_size = length(Segmentpix_T1{slc});
    for seg = 1:seg_size
        seg_mean_t1(count_t1) = mean(Segmentpix_T1{slc}{seg});
        seg_mean_temp = Segmentpix_ff{slc}{seg};
        seg_mean_temp(seg_mean_temp < 0) = 0;
        seg_mean_temp(isnan(seg_mean_temp)) = 0;
        seg_mean_temp(seg_mean_temp > 100) = 100;
        seg_mean_ff(count_ff) = mean(seg_mean_temp);
        
        count_t1 = count_t1 + 1;
        count_ff = count_ff + 1;
    end
    
end


figure(); scatter(seg_mean_t1, seg_mean_ff); grid on;
figure(); yyaxis left; plot(seg_mean_t1); 
hold on; yyaxis right; plot(seg_mean_ff); 
grid on;

%% Apply mi to myo mask
clear ff_mean_mod t1_mean_mod
mi_3D_new = mi_3D .* myo_3D;
mi_ff_3D_new = mi_ff_3D .* myo_ff_3D;
Mask_Segn_ff_cell = cell(size(Mask_Segn_ff,3), 1);
Mask_Segn_t1_cell = cell(size(Mask_Segn_t1,3), 1);
count = 1;
count_t1 = 1;

count_array_t1 = zeros(1, size(Mask_Segn_ff, 3));
count_array_ff = zeros(1, size(Mask_Segn_ff, 3));

for slc = 1:size(Mask_Segn_ff, 3)
    ff = load(ff_glob{idx_array(slc)});
    ff_map = ff.fwmc_ff;
    Mask_Segn_ff_3D = SpreadLabels(Mask_Segn_ff(:,:,slc));
    Mask_Segn_t1_3D = SpreadLabels(Mask_Segn_t1(:,:,slc));
    Mask_Segn_ff_cell{slc} = Mask_Segn_ff_3D;
    Mask_Segn_t1_cell{slc} = Mask_Segn_t1_3D;
    inner_count_t1 = 0;
    inner_count_ff = 0;
    fig1 = figure();
    imagesc(Mask_Segn_ff(:,:,slc) .* mi_ff_3D_new(:,:,slc));
    fig2 = figure();
    imagesc(Mask_Segn_t1(:,:,slc) .* mi_3D_new(:,:,slc));
    for i = 1:size(Mask_Segn_ff_3D, 3)
        t1_map_masked = Mask_Segn_t1_3D(:,:,i) .* mi_3D_new(:,:,slc) .* img_3D(:,:,slc);
        ff_map_masked = Mask_Segn_ff_3D(:,:,i) .* mi_ff_3D_new(:,:,slc) .* ff_map;
        ff_map_masked(isnan(ff_map_masked)) = 0;
        t1_map_masked(isnan(t1_map_masked)) = 0;
        ff_map_masked(ff_map_masked<0) = 0;
        ff_map_masked(ff_map_masked>100) = 100;
        if ~isempty(nonzeros(ff_map_masked))
            
            ff_mean_mod(count) = mean(nonzeros(ff_map_masked));
            inner_count_ff = inner_count_ff + 1;
            count = count + 1;
        end
        if ~isempty(nonzeros(t1_map_masked))
            t1_mean_mod(count_t1) = mean(nonzeros(t1_map_masked));
            inner_count_t1 = inner_count_t1 + 1;
            count_t1 = count_t1 + 1;
            
        end
    end
    count_array_t1(slc) = inner_count_t1;
    count_array_ff(slc) = inner_count_ff;
end

figure(); scatter(ff_mean_mod, t1_mean_mod, 64); grid on;
hold on;
plot(ff_mean_mod, k*ff_mean_mod+interc);
figure(); yyaxis left; plot(t1_mean_mod); 
hold on; yyaxis right; plot(ff_mean_mod); 
grid on;

