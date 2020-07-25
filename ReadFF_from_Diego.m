clear all;
close all;
%%
addpath('D:\src\function');
base_dir = uigetdir;

folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

ll = 7;
label = labels{ll};
idx_array = contains(folder_glob, label);

if any(idx_array)
    num = numel(nonzeros(idx_array));
    ind_array = find(idx_array == 1);
    dst_files = cell(num, 1);
    for i = 1:num
        dst_files{i} = folder_glob{ind_array(i)};
    end
    strings = strsplit(dst_files{1});
    dst_name = strings{end};
    disp(dst_name);
    
    %sel_array = input('Please add an array here:  ');
    
    %char_array = num2str(sel_array', '%04.f');
    
    %ind_array2 = zeros(size(dst_name, 1), 1);
    %for i = 1:size(char_array, 1)
    %    cha = char_array(i, :);
    %    ind_array2 = ind_array2 + contains(dst_name, cha);
    %end
    
    %ind_array3 = find(ind_array2 == 1);
    %list_to_read = dst_files(ind_array3);
    
    %name_to_compare = ExtractNames(list_to_read);
    
    %order_to_read = zeros(length(list_to_read), 1);
    %for i = 1:length(list_to_read)
    %    order_to_read(i) = find(contains(name_to_compare, char_array(i, :)) == 1);
    %end
    
    save_dir = GetFullPath(cat(2, base_dir, '\..\img\'));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
    
    load(dst_name, 'fwmc_ff');
end

imshow3D(fwmc_ff, [0 20])
%% Show mGRE
dicom_dir = uigetdir('D:\Data\');
labels_to_read = [1];
folder_glob = glob(cat(2, dicom_dir, '\*'));

ll = labels_to_read(1);
label = labels{ll};
idx_array = contains(folder_glob, label);

num = numel(nonzeros(idx_array));
ind_array = find(idx_array == 1);
dst_files = cell(num, 1);


for i = 1:num
    dst_files{i} = folder_glob{ind_array(i)};
end


dst_name = ExtractNames(dst_files);
disp(dst_name);

sel_array = input('Please add an array here:  ');

char_array = num2str(sel_array', '%04.f');

ind_array2 = zeros(size(dst_name, 1), 1);
for i = 1:size(char_array, 1)
    cha = char_array(i, :);
    ind_array2 = ind_array2 + contains(dst_name, cha);
end

ind_array3 = find(ind_array2 == 1);
list_to_read = dst_files(ind_array3);

name_to_compare = ExtractNames(list_to_read);

order_to_read = zeros(length(list_to_read), 1);
for i = 1:length(list_to_read)
    order_to_read(i) = find(contains(name_to_compare, char_array(i, :)) == 1);
end

save_dir = GetFullPath(cat(2, dicom_dir, '\..\img\'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
%%
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit = dicom23D(f);
    loc = 6;
    
    figure(); imagesc(whatsinit(:,:,loc));
    epi_t1 = drawpolygon(gca);
    endo_t1 = drawpolygon(gca);
    
    %
    mask_epi_t1 = createMask(epi_t1);
    mask_endo_t1 = createMask(endo_t1);
    mask_myocardium_t1 = mask_epi_t1 - mask_endo_t1;
    figure(); 
    imagesc(mask_myocardium_t1.* whatsinit(:,:,loc))
end
%%


mGRE = reshape(volume_image, size(volume_image, 1), size(volume_image, 2), 8, []);

figure(); imshow3D(mGRE(:,:,:,20));
%%
% % draw on slice 1
% figure(); imagesc(mGRE(:,:,1,20));
% heart = drawpolygon(gca);
% heart_coords1 = heart.Position;
% mask1 = createMask(heart);

% draw on slice 8
figure(); imagesc(mGRE(:,:,8,20));
epi_ff8 = drawpolygon(gca);
epi_coords8 = epi_ff8.Position;
endo_ff8 = drawpolygon(gca);
endo_coords8 = endo_ff8.Position;
mask_epi_ff8 = createMask(epi_ff8);
mask_endo_ff8 = createMask(endo_ff8);


%%
mask_myocardium_ff8 = mask_epi_ff8 - mask_endo_ff8;
figure(); imagesc(mask_myocardium_ff8 .* mGRE(:,:,8,20));
figure(); imagesc(mask_myocardium_ff8 .* fwmc_ff(:,:,20)); caxis([0 20])
remote_ff8 = drawpolygon(gca);
mask_remote_ff8 = createMask(remote_ff8);
%%
figure(); imagesc(mask_myocardium_ff8 .* fwmc_ff(:,:,20)); caxis([0 20])
hemo_ff8 = drawpolygon(gca);
mask_hemo_ff8 = createMask(hemo_ff8);

%mean(nonzeros(mask_remote_ff8 .* ff)) 2.1622
%mean(nonzeros(mask_hemo_ff8.*mask_myocardium_ff8 .* ff)) 6.6504
%%
figure(); imagesc(mGRE(:,:,1,20).*mask1);
figure(); imagesc(mGRE(:,:,8,20).*mask8);
%%
figure(); imagesc(fwmc_ff(:,:,20).*mask1); caxis([0 20]);
figure(); imagesc(fwmc_ff(:,:,20).*mask8); caxis([0 20]);