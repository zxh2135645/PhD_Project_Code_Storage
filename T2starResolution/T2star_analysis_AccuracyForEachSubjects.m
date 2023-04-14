clear all;
close all;

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% aha_analysis
% T2star_meanSD_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob);

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name = input('Please type subject name here:  ', 's');
subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
if ~exist(subject_dir, 'dir')
    mkdir(subject_dir)
end
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

disp('Avg 0016 starts here: ');
avg_num = input('Please type average number here:  ');
avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));

%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100]);
end

%% Draw contours @ epi, endo, MI, remote, fluid
img = whatsinit{1};
myo_coords_cell = cell(size(img, 3), 2);
roi_save = cat(2, subject_data_dir, 'roi.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(img, 3))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    air = drawpolygon(gca);
    air_coords = air.Position;
    
    myo_coords_cell{i, 1} = epi.Position;
    myo_coords_cell{i, 2} = endo.Position;
    epi_mask = createMask(epi);
    endo_mask = createMask(endo);
    center_mask = createMask(center_line);
    
    close all;
end

roi.myo_coords_cell = myo_coords_cell;
roi.mi_coords = mi_coords;
roi.remote_coords = remote_coords;
roi.fluid_coords = fluid_coords;
roi.center_coords = center_coords;
roi.air_coords = air_coords;

save(roi_save, 'roi');

else
    load(roi_save);
    myo_coords_cell = roi.myo_coords_cell;
    mi_coords = roi.mi_coords;
    remote_coords = roi.remote_coords;
    fluid_coords = roi.fluid_coords;
    center_coords = roi.center_coords;
    air_coords = roi.air_coords;
end

% 28 different set of parameters
% Convert coords to masks for 28 images
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    for i = 1:length(whatsinit)
        img2 = whatsinit{i};
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); caxis([0 100]);
        epi = drawpolygon(gca,'Position', [myo_coords_cell{1}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{1}(:,2)/ratio + (ratio-1)/ratio]);
        endo = drawpolygon(gca,'Position', [myo_coords_cell{2}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{2}(:,2)/ratio + (ratio-1)/ratio]);
        mi = drawpolygon(gca,'Position', [mi_coords(:,1)/ratio + (ratio-1)/ratio, mi_coords(:,2)/ratio + (ratio-1)/ratio]);
        remote = drawpolygon(gca,'Position', [remote_coords(:,1)/ratio + (ratio-1)/ratio, remote_coords(:,2)/ratio + (ratio-1)/ratio]);
        fluid = drawpolygon(gca,'Position', [fluid_coords(:,1)/ratio + (ratio-1)/ratio, fluid_coords(:,2)/ratio + (ratio-1)/ratio]);
        air = drawpolygon(gca,'Position', [air_coords(:,1)/ratio + (ratio-1)/ratio, air_coords(:,2)/ratio + (ratio-1)/ratio]);
        center_line = drawpolygon(gca,'Position', [center_coords(:,1)/ratio + (ratio-1)/ratio, center_coords(:,2)/ratio + (ratio-1)/ratio]);
        
        epi_mask = createMask(epi);
        endo_mask = createMask(endo);
        myo_mask = epi_mask - endo_mask;
        
        mi_mask = createMask(mi);
        remote_mask = createMask(remote);
        fluid_mask = createMask(fluid);
        air_mask = createMask(air);
        center_mask = createMask(center_line);
        
        mask_struct(i).myo_mask = myo_mask;
        mask_struct(i).mi_mask = mi_mask;
        mask_struct(i).remote_mask = remote_mask;
        mask_struct(i).fluid_mask = fluid_mask;
        mask_struct(i).air_mask = air_mask;
        
        mask_struct(i).epi_mask = epi_mask;
        mask_struct(i).endo_mask = endo_mask;
        
        myo_mask_endo = myo_mask .* center_mask;
        myo_mask_epi = myo_mask - myo_mask_endo;
        mask_struct(i).myo_mask_endo = myo_mask_endo;
        mask_struct(i).myo_mask_epi = myo_mask_epi;
    end
    
    save(mask_save, 'mask_struct');
else
    load(mask_save);
end

%% Display image overlay with Mean-2SD mask
save_array = 1:1:length(whatsinit);
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    img2 = whatsinit{save_idx};
    thresh = mean(nonzeros(img2 .* mask_struct(save_idx).remote_mask)) - 2*std(nonzeros(img2 .* mask_struct(save_idx).remote_mask));
    hemo_mask = img2 < thresh;
    hemo_in_mi = hemo_mask.*mask_struct(save_idx).mi_mask.*mask_struct(save_idx).myo_mask;
    mi_perc = numel(nonzeros(hemo_in_mi)) ./ numel(nonzeros(mask_struct(save_idx).myo_mask))
end

%% Analysis
hemo_thickness = [2.40, 0.55, 0.69, 0.56, 1.21, 0.59, 0.53, 0.92, 0.59, 1.10];
accuracy = [0.79, 0.79, 0.64, 0.89, 0.73, 0.63, 0.94, 0.71, 0.58, 0.74];
snr = [13.11, 13.09, 18.09, 14.69, 20.51, 18.28, 17.23, 15.26, 17.45, 16.56];
bias_invivo = [0.25, 0.166666667, 1.155688623, 1, 0.145038168, 0.888888889, 1, 0.090909091, 3.482758621, 0.568788501];
snr_invivo = [4.46, 4.76, 11.02, 7.28, 5.21, 8.53, 4.99, 7.41, 8.81, 11.50];
figure();
scatter3(hemo_thickness, accuracy, snr);
xlabel('Hemo Thickness (mm)');
ylabel('Accuracy');
zlabel('SNR');

figure();
scatter3(hemo_thickness, bias_invivo, snr_invivo);
xlabel('Hemo Thickness (mm)');
ylabel('Bias');
zlabel('SNR Invivo');
%%
addpath('../function/')
figure();
% plotclr(hemo_thickness, accuracy, snr, 'd')
x = hemo_thickness;
y = accuracy;
y = bias_invivo;
v = snr;
v = snr_invivo;
map=colormap;
miv=min(v);
mav=max(v);
clrstep = (mav-miv)/size(map,1) ;
% Plot the points
hold on
for nc=1:size(map,1)
    iv = find(v>miv+(nc-1)*clrstep & v<=miv+nc*clrstep) ;
    plot3(x(iv),y(iv),v(iv),'o', 'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',12)
end
hold off
% Re-format the colorbar
h=colorbar;
set(h,'ylim',[1 length(map)]);
yal=linspace(1,length(map),10);
set(h,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,10);
s=char(10,4);
for i=1:10
    if min(abs(ytl)) >= 0.001
        B=sprintf('%-4.3f',ytl(i));
    else
        B=sprintf('%-3.1E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(h,'yticklabel',s);
grid on