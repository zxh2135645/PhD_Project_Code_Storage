clear all;
close all;

%% Read RYN T1 maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = 'D:\Data\London_Canine\';
read_dir = cat(2, base_dir, 'RYN_WK8_29NOV18\HEART_CEDARS_20181129_090537_179000\');
save_dir = cat(2, base_dir, 'RYN_WK8_29NOV18\Figures\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

%% Read files
t1_glob = glob(cat(2, read_dir, 'T1MAP_*'));
idx_array = zeros(length(t1_glob), 1);
for i = 1:length(t1_glob)
    strings = split(t1_glob{i}, '\');
    fname = strings{6};
    strings = split(fname, '_');
    idx_array(i) = str2num(strings{end});
end

manual_pick = [53]; % sax 6
for i = 1:length(manual_pick)
    idx = find(manual_pick(i) == idx_array);
    
    [images_apex, headers_apex] = dicomfolder(t1_glob{idx});
    if i == 1
       img_flow_apex = zeros(size(images_apex, 1), size(images_apex, 2), size(images_apex, 3), length(manual_pick)); 
       hdr_flow_apex = cell(size(images_apex, 3), length(manual_pick));
    end
    
    img_flow_apex(:,:,:,i) = images_apex;
    for j = 1:length(headers_apex)
        hdr_flow_apex{j,i} = headers_apex{j};
    end
end

%% Draw ROI
figure();
imagesc(img_flow_apex); axis equal;
colormap gray;

%% MT from ex-vivo heart in 10% formalin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Degree 45 and 15 (Heart 104 slice)
base_dir = 'D:\Data\Exvivo_Phantom\';
read_dir = cat(2, base_dir, 'JAMES_PHANTOM_EXVIVOHEART104_SLICE\BIRI_RESEARCH_JAMES_20190801_214319_715000\');
save_dir = cat(2, base_dir, 'JAMES_PHANTOM_EXVIVOHEART104_SLICE\Figures\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
glob_str = 'A_CV_STARS_CINE_FA*_BW1002_*';
[img_flow, hdr_flow] = ImgPick(read_dir, glob_str);
% [17 19]; % Fast44, LowSAR44

avg_label = 0; % Average over cardiac CINE phases
fname = 'roi_averaged.mat'; % The file name to save
slice_len = 1;
fsax = '';
[img_mean, roi_epi, roi_endo, roi_myo] = ...
    DrawROI(img_flow, save_dir, fname, avg_label, slice_len);
[img_flow_masked, img_flow_masked_blood] = ...
    MultiplyROI(img_mean, roi_epi, roi_myo, save_dir, fsax);

%% We don't see infarct in formalin
%% Show image
% figure();
% img_mean = zeros(size(img_flow,1), size(img_flow,2), size(img_flow,4));
% for i = 1:size(img_flow,4)
%     img_mean(:,:,i) = mean(img_flow(:,:,:,i), 3);
%     subplot(1,2,i)
%     imagesc(img_mean(:,:,i))
%     axis equal; axis off;
% end 

% Double check masks (Looks pretty OK)
% figure();
% for i = 1:size(img_mean, 3)
%     imagesc(img_flow_masked_apex_blood(:,:,i))
%     axis equal;
%     title(num2str(i))
%     pause(0.5)
% end
%% MTR from London dogs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gobi - Degree 45
base_dir = 'D:\Data\London_Canine\';
% read_glob = glob(cat(2, base_dir, '*\*\A_CV_STARS_FA*_BW1002_FAST_SA*'));
read_dir = cat(2, base_dir, 'GOBI_TEST_22JUL2019\HEART_DIANE_20190722_143948_306000\');
save_dir = cat(2, base_dir, 'GOBI_TEST_22JUL2019\Figures\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
glob_str = 'CINE_FA45_BW1002_*';
[img_flow, hdr_flow] = ImgPick(read_dir, glob_str);
% [17 19]; % LowSAR1, LowSAR2, ..., Fast1, Fast2, ...

avg_label = 33; % Average over cardiac CINE phases
fname = 'roi_cine.mat'; % The file name to save
slice_len = 1;

[img_mean, roi_epi, roi_endo, roi_myo] = ...
    DrawROI(img_flow, save_dir, fname, avg_label, slice_len);
fsax = 1:1:(size(img_mean,3)/2);
[img_flow_masked, img_flow_masked_blood] = ...
    MultiplyROI(img_mean, roi_epi, roi_myo, save_dir, fsax);
%% The rest of 6 dogs - Degree 45
base_dir = 'D:\Data\London_Canine\';
glob_str = 'A_CV_STARS_FA*_BW1002_*_SA*';
read_glob = glob(cat(2, base_dir, '*\*\', glob_str));
dog_names = {'MERRY', 'SUNNY', 'RYN', 'MOJAVE', 'SAHARA', 'HOPE'};
for i = 1:length(dog_names)
    dog = dog_names{i};
    index = find(contains(read_glob, dog));
    dog_glob = read_glob(index);
    last_pos = find(dog_glob{1} == '\', 2, 'last');
    read_dir = dog_glob{1}(1:last_pos(1));
    strings = split(read_dir, '\');
    save_dir = cat(2, base_dir, strings{4}, '\Figures\');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
    index = find(contains(dog_glob, 'ADJ'));
    dog_glob(index, :) = [];
    disp("=========================================================================================")
    disp(["The dog you are processing:", dog_names{i}])
    [img_flow, hdr_flow] = ImgPick(read_dir, glob_str);
    % % LowSAR1, LowSAR2, ..., Fast1, Fast2, ...
    strings_split = split(dog_glob, '_');
    sax_array = unique(strings_split(:,end-1));
    avg_label = 33; % Average over cardiac CINE phases
    fname = 'roi_cine.mat'; % The file name to save
    slice_len = 1;
    [img_mean, roi_epi, roi_endo, roi_myo] = ...
        DrawROI(img_flow, save_dir, fname, avg_label, slice_len);
    % fsax = 1:1:(size(img_mean,3)/2);
    fsax = sax_array;
    [img_flow_masked, img_flow_masked_blood] = ...
        MultiplyROI(img_mean, roi_epi, roi_myo, save_dir, fsax);
end

%% Image registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Register images for HOPE
dog = 'HOPE';
index = find(contains(read_glob, dog));
dog_glob = read_glob(index);
last_pos = find(dog_glob{1} == '\', 2, 'last');
read_dir = dog_glob{1}(1:last_pos(1));
strings = split(read_dir, '\');
[img_flow, hdr_flow] = ImgPick(read_dir, glob_str);


% Remove name with ADJ to have consistent split
index = find(contains(dog_glob, 'ADJ'));
dog_glob(index, :) = [];
strings_split = split(dog_glob, '_');
sax_array = unique(strings_split(:,end-1));
avg_label = 33; % Average over cardiac CINE phases
fname = 'roi_cine_ForReg.mat'; % The file name to save
slice_len = 1;
save_dir = cat(2, base_dir, strings{4}, '\Figures\');

% Draw Epi rois
[img_mean, roi_flow] = ...
        DrawROI_Generic(img_flow, save_dir, fname, avg_label, slice_len);

%%  Registration (rigid)
pairs = size(img_mean, 3) / 2;
img_reg_flow = zeros(size(img_mean, 1), size(img_mean, 2), pairs);
mtr_flow = zeros(size(img_mean, 1), size(img_mean, 2), pairs);
neo_img_mean = img_mean;
for i = 1:pairs
    FIXED = roi_flow(:,:,i);
    MOVING = roi_flow(:,:,i+pairs);
    [MOVINGREG, fixedRefObj, movingRefObj] = registerImages(MOVING,FIXED);
    img_fixed = img_mean(:,:,i);
    img_moving = img_mean(:,:,i+pairs);
    tform = MOVINGREG.Transformation;
    img_reg_flow(:,:,i) = imwarp(img_moving, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
    mtr_flow(:,:,i) = (img_fixed - img_reg_flow(:,:,i)) ./ img_fixed;
    
    figure();
    imagesc(mtr_flow(:,:,i).*epi_flow(:,:,i)); axis equal; axis off; caxis([0 0.5]);
    
    neo_img_mean(:,:,i+pairs) = img_reg_flow(:,:,i);
end

fname = 'roi_cine_epi.mat';
% Draw Epi rois
[~, epi_flow] = ...
    DrawROI_Generic(img_flow, save_dir, fname, avg_label, slice_len);
roi_epi = epi_flow(:,:,1:pairs);
roi_endo = roi_flow(:,:,1:pairs);
roi_myo = (roi_epi + roi_endo == 1);
strings_split = split(dog_glob, '_');
sax_array = unique(strings_split(:,end-1));
fsax = sax_array;
[img_flow_masked, img_flow_masked_blood] = ...
    MultiplyROI(neo_img_mean, roi_epi, roi_myo, save_dir, fsax);

%% See CNR on MTR map for Merry
%% 09/30/2019
dog = 'MERRY';
index = find(contains(read_glob, dog));
dog_glob = read_glob(index);
last_pos = find(dog_glob{1} == '\', 2, 'last');
read_dir = dog_glob{1}(1:last_pos(1));
strings = split(read_dir, '\');
[img_flow, hdr_flow] = ImgPick(read_dir, glob_str);
%%
save_dir = cat(2, base_dir, strings{4}, '\Figures\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
% LowSAR1, LowSAR2, ..., Fast1, Fast2, ...
strings_split = split(dog_glob, '_');
sax_array = unique(strings_split(:,end-1));
avg_label = 33; % Average over cardiac CINE phases
fname = 'roi_cine.mat'; % The file name to save
slice_len = 1;
[img_mean, roi_epi, roi_endo, roi_myo] = ...
    DrawROI(img_flow, save_dir, fname, avg_label, slice_len);
% fsax = 1:1:(size(img_mean,3)/2);
fsax = sax_array;
[img_flow_masked, img_flow_masked_blood] = ...
    MultiplyROI(img_mean, roi_epi, roi_myo, save_dir, fsax);

%%
figure();
for i = 1:size(img_flow_masked_blood, 3)
    subplot(3,4,i)
    imagesc(img_flow_masked_blood(:,:,i))
    axis equal
end

%% Get MTR map
pairs = size(img_flow_masked_blood, 3) / 2;
mi_roi = zeros(size(img_flow_masked_blood, 1), size(img_flow_masked_blood, 2), pairs);
remote_roi = zeros(size(img_flow_masked_blood, 1), size(img_flow_masked_blood, 2), pairs);
sMT45_blood_flow = zeros(size(img_flow_masked_blood, 1), size(img_flow_masked_blood, 2), pairs);
for i = 1:pairs
    s0 = img_flow_masked_blood(:,:,i);
    sMT45_blood = (s0 - img_flow_masked_blood(:,:,i+pairs)) ./ s0;
    sMT45_blood(sMT45_blood < 0) = 0;
    sMT45_blood(isnan(sMT45_blood)) = 0;
    sMT45_blood_flow(:,:,i) = sMT45_blood;
    
    for_mask = sMT45_blood; % at the LowSAR
    rescaled_for_mask = for_mask / (max(for_mask(:)) - min(for_mask(:)));
    
    figure();
    mi_roi(:,:,i) = roipoly(rescaled_for_mask);
    remote_roi(:,:,i) = roipoly(rescaled_for_mask);
end

%% If not smoothed
mi = mi_roi .* sMT45_blood_flow;
remote = remote_roi .* sMT45_blood_flow;
cnr_array = zeros(1, pairs);
mtr_mi = zeros(1, pairs);
mtr_remote = zeros(1, pairs);
std_remote = zeros(1, pairs);

for i = 1:length(cnr_array)
    cnr_array(i) = (mean(nonzeros(mi(:,:,i)))-mean(nonzeros(remote(:,:,i)))) / std(nonzeros(remote(:,:,i)));
    mtr_mi(i) = mean(nonzeros(mi(:,:,i)));
    mtr_remote(i) = mean(nonzeros(remote(:,:,i)));
    std_remote(i) = std(nonzeros(remote(:,:,i)));
end

% The original CNR is not big enough

%% If smoothed
sMT45_blood_medsmth = zeros(size(sMT45_blood_flow));
for i = 1:pairs
    sMT45_blood_medsmth(:,:,i) = medfilt2(sMT45_blood_flow(:,:,i), [3 3]);
end

mi = mi_roi .* sMT45_blood_medsmth;
remote = remote_roi .* sMT45_blood_medsmth;
cnr_array = zeros(1, pairs);
mtr_mi = zeros(1, pairs);
mtr_remote = zeros(1, pairs);
std_remote = zeros(1, pairs);

for i = 1:length(cnr_array)
    cnr_array(i) = (mean(nonzeros(mi(:,:,i)))-mean(nonzeros(remote(:,:,i)))) / std(nonzeros(remote(:,:,i)));
    mtr_mi(i) = mean(nonzeros(mi(:,:,i)));
    mtr_remote(i) = mean(nonzeros(remote(:,:,i)));
    std_remote(i) = std(nonzeros(remote(:,:,i)));
end
