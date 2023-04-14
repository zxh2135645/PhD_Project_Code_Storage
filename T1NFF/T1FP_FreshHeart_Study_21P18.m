close all; clear all;
% Specifically studying 21P18 Exvivo fresh heart scan
addpath('../function/')
dicom_dir = uigetdir;
dicom_glob = glob(cat(2, dicom_dir, '/T1_TSE_SAG_320*'));
%% MTR
dicom_glob_reorder = dicom_glob([1,2]);
dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    'TriggerTime',...
    'RepetitionTime',...
    'EchoTime',...
    'EchoTrainLength',...
    'InversionTime'
    };

whatsinit = cell(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end

img_mt = whatsinit{1}(:,:,:);
img_nomt = whatsinit{2}(:,:,2:end);
img_mtr = (img_nomt - img_mt) ./ img_nomt;
figure();
for i = 1:size(img_mtr, 3)
    subplot(4,4,i)
    imagesc(abs(img_mtr(:,:,i))); caxis([0 0.5]); 
end

% img_mtr2 = img_mt ./ img_nomt; % From K Haliot
% figure();
% for i = 1:size(img_mtr, 3)
%     subplot(4,4,i)
%     imagesc(img_mtr2(:,:,i)); caxis([0 1]); 
% end
% MTR is lowered in MI for this case

%% MRS
addpath('../function/');
mrs_glob = glob(cat(2, dicom_dir, '/SVS_*'));

water_ppm = 4.7;
ppm_ref = water_ppm;

disp(mrs_glob);
idx_array = input('Please refer to dir_glob for the index you want to import \n Enter cases with [ ] around them (Before | After):  '); % Before, After water suppression

InstanceNum = 1;
ss = cell(length(idx_array), 1);

for i = 1:length(idx_array)
    dst_dir = mrs_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    SeriesNum = str2num(f_parts{end});
    concent = f_parts{end-1};
    
    % Need to run startup.m in  
    ss{i} = loadSVSData(dst_dir, SeriesNum, InstanceNum); % Doesn't work for current version
    
    if i == 1
        s = ss{i};
        dwt = s.dwellTime; % in the unit of seconds
        n = s.samples;  % number of samples
        w = s.freqAxis; % frequency axis (Hz)
        t = s.timeAxis; % time axis (s)
        ppm_axis = s.ppmAxis + ppm_ref;
    end
    
    spectrum = ss{i}.spectra{1};
end

figure();
plot(abs(ss{1}.spectra{1}));
hold on;
plot(abs(ss{2}.spectra{1}));
legend({'MI', 'Remote'});

% Get the signal percentage
N = linspace(1,n,n);
perc = zeros(size(ss,1), 1);
cellToSave = cell(size(ss,1), 2);

for i = 1:size(ss, 1)
    dst_dir = mrs_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    
    spectrum = ss{i}.spectra{1};
    Segment_Area = cumtrapz(abs(spectrum));
    
    % Assuming the bandwidth of fatsat pulse is 200Hz thus 1.6 ppm width
    % Thus the fatsat pulse is at 0.4 - 2.0 ppm
    ppm_temp = abs(ppm_axis - 0.4);
    idx_low = find(min(ppm_temp) == ppm_temp);
    
    ppm_temp = abs(ppm_axis - 2.0);
    idx_high = find(min(ppm_temp) == ppm_temp);
    
    Fatsat_Area = Segment_Area(idx_high) - Segment_Area(idx_low);
    Total_Area = Segment_Area(end);
    
    perc(i) = Fatsat_Area / Total_Area;
    cellToSave{i, 1} = f_parts{end-1};
    cellToSave{i, 2} = perc(i);
end

% Plot spectrum
% No water suppression on the left
% Water suppression on the right
lim_mat = zeros(size(ss, 1), 2);
for i = 1:size(ss, 1)
    spectrum = ss{i}.spectra{1};
    lim_mat(i, 1) = max(abs(spectrum)); % max
    lim_mat(i, 2) = min(abs(spectrum)); % min
end

[mmax, mx] = max(max(lim_mat, [], 2));
[mmin, mm] = min(max(lim_mat, [], 2));
fixed_yl = ceil((mmax - 0) * 1.05);

% Original YLim
figure();
for i = 1:size(ss, 1)
    dst_dir = mrs_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    concent = f_parts{end-1};
    
    spectrum = ss{i}.spectra{1};
    
    subplot(1,2,i)
    plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity');
    hold on;
    rectangle('Position', [ppm_axis(idx_low), 0, ppm_axis(idx_high)-ppm_axis(idx_low), max(abs(spectrum))], 'FaceColor', [1 0 0 0.1])
    %plot(ppm_axis, medfilt1(abs(spectrum), 10));
    H = area(ppm_axis(idx_low:idx_high), abs(spectrum(idx_low:idx_high)));
    set(H(1), 'FaceColor', [0.8500    0.3250    0.0980])
    set(gca, 'Xdir', 'reverse')
    title(cat(2, 'MR Spectroscopy ', concent))
    grid on;
end

% Fixed YLim
figure();
for i = 1:size(ss, 1)
    dst_dir = mrs_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    concent = f_parts{end-1};
    
    spectrum = ss{i}.spectra{1};
    
    subplot(1,2,i)
    plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity');
    hold on;
    rectangle('Position', [ppm_axis(idx_low), 0, ppm_axis(idx_high)-ppm_axis(idx_low), max(abs(spectrum))], 'FaceColor', [1 0 0 0.1])
    %plot(ppm_axis, medfilt1(abs(spectrum), 10));
    H = area(ppm_axis(idx_low:idx_high), abs(spectrum(idx_low:idx_high)));
    set(H(1), 'FaceColor', [0.8500    0.3250    0.0980])
    set(gca, 'Xdir', 'reverse')
    title(cat(2, 'MR Spectroscopy ', concent))
    ylim([0 fixed_yl])
    grid on;
end


yl_cuttingwater(mx) = ceil(0.02*max(max(lim_mat, [], 2)));
yl_cuttingwater(mm) = ceil(0.1*min(max(lim_mat, [], 2)));
% Cutting water peak YLim
figure();
for i = 1:size(ss, 1)
    dst_dir = mrs_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    concent = f_parts{end-1};
    
    spectrum = ss{i}.spectra{1};
    
    subplot(1,2,i)
    plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity');
    hold on;
    rectangle('Position', [ppm_axis(idx_low), 0, ppm_axis(idx_high)-ppm_axis(idx_low), max(abs(spectrum))], 'FaceColor', [1 0 0 0.1])
    %plot(ppm_axis, medfilt1(abs(spectrum), 10));
    H = area(ppm_axis(idx_low:idx_high), abs(spectrum(idx_low:idx_high)));
    set(H(1), 'FaceColor', [0.8500    0.3250    0.0980])
    set(gca, 'Xdir', 'reverse')
    title(cat(2, 'MR Spectroscopy ', concent))
    ylim([0 yl_cuttingwater(i)])
    grid on;
end

% MRS of remote is much lower, and no fat peak is shown
%% T1-MAPIT, T2-MAPIT, T2STAR
dicom_glob1 = glob(cat(2, dicom_dir, '/T1_IMAGES_0011'));
dicom_glob2 = glob(cat(2, dicom_dir, '/T2_IMAGES_0013'));
dicom_glob3 = glob(cat(2, dicom_dir, '/T2STAR_IMAGES_0015'));

whatsinit = cell(length(dicom_glob1), 1);
for i = 1:length(dicom_glob1)
    f = dicom_glob1{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t1_mapit = whatsinit{1};

whatsinit = cell(length(dicom_glob2), 1);
for i = 1:length(dicom_glob2)
    f = dicom_glob2{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t2_mapit = whatsinit{1};

whatsinit = cell(length(dicom_glob3), 1);
for i = 1:length(dicom_glob3)
    f = dicom_glob3{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t2star_mapit = whatsinit{1};

%%
figure();
subplot(2,2,1);
imagesc(img_t1_mapit(:,:,10)); caxis([0 2000]); axis image; title('T1 (ms)');
subplot(2,2,2);
imagesc(img_t2_mapit(:,:,10)); caxis([0 200]); axis image; title('T2 (ms)');
subplot(2,2,3);
imagesc(img_t2star_mapit(:,:,10)); caxis([0 100]); axis image; title('T2star (ms)');

mask_f = GetFullPath(cat(2, dicom_dir, '/../mask.mat'));
mask_roi = zeros(size(whatsinit{1},1), size(whatsinit{1},2));
mask_remote = zeros(size(whatsinit{1},1), size(whatsinit{1},2));
mask = struct;
if ~exist(mask_f)
    for i = 1:1
        % for i = 1:size(whatsinit{1},3)

        figure();
        imagesc(img_t2_mapit(:,:,10)); axis image; caxis([0 200]);
        roi = drawpolygon;
        mask_roi(:,:,i) = createMask(roi);

        imagesc(img_t2_mapit(:,:,10)); axis image; caxis([0 200]);
        remote = drawpolygon;
        mask_remote(:,:,i) = createMask(remote);
    end
    mask.mask_roi = mask_roi;
    mask.mask_remote = mask_remote;
    save(mask_f, '-struct', 'mask');
else
    load(mask_f);
end

mask_hemo = (img_t2star_mapit(:,:,10) < thresh).*mask_roi;
thresh = mean(nonzeros(img_t2star_mapit(:,:,10).*mask_remote)) - 2*std(nonzeros(img_t2star_mapit(:,:,10).*mask_remote));
figure(); imagesc(mask_hemo); %caxis([0 50]);

t1_array_hemo = nonzeros(img_t1_mapit(:,:,10) .* mask_hemo);
t2_array_hemo = nonzeros(img_t2_mapit(:,:,10) .* mask_hemo);
t2star_array_hemo = nonzeros(img_t2star_mapit(:,:,10) .* mask_hemo);

t1_array_nonhemo = nonzeros(img_t1_mapit(:,:,10) .* (mask_roi -mask_hemo));
t2_array_nonhemo = nonzeros(img_t2_mapit(:,:,10) .* (mask_roi -mask_hemo));
t2star_array_nonhemo = nonzeros(img_t2star_mapit(:,:,10) .* (mask_roi -mask_hemo));

figure();
p = plot3(t1_array_hemo, t2_array_hemo, t2star_array_hemo); 
p.LineStyle = "none";
p.Color = "red";
p.Marker = "o";
hold on;
p1 = plot3(t1_array_nonhemo, t2_array_nonhemo, t2star_array_nonhemo); 
grid on;
p1.LineStyle = "none";
p1.Color = "blue";
p1.Marker = "x";
xlabel('T1 (ms)'); ylabel('T2 (ms)'); zlabel('T2star (ms)');
%%
dicom_glob1 = glob(cat(2, dicom_dir, '/T1MAP_PRECON_SAX7_MOCO_T1_0036'));
dicom_glob2 = glob(cat(2, dicom_dir, '/T2MAP_FLASH_SAX9_MOCO_T2_0078'));
dicom_glob3 = glob(cat(2, dicom_dir, '/T2_MAP_FLASH_SAX9_T2STAR_0115'));

whatsinit = cell(length(dicom_glob1), 1);
for i = 1:length(dicom_glob1)
    f = dicom_glob1{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t1_molli = whatsinit{1};

whatsinit = cell(length(dicom_glob2), 1);
for i = 1:length(dicom_glob2)
    f = dicom_glob2{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t2_cmr = whatsinit{1};

whatsinit = cell(length(dicom_glob3), 1);
for i = 1:length(dicom_glob3)
    f = dicom_glob3{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
end
img_t2star_cmr = whatsinit{1};

figure();
subplot(2,2,1);
imagesc(img_t1_molli); caxis([0 2000]); axis image; title('T1 (ms)');
subplot(2,2,2);
imagesc(img_t2_cmr*0.1); caxis([0 200]); axis image; title('T2 (ms)');
subplot(2,2,3);
imagesc(img_t2star_cmr*0.1); caxis([0 100]); axis image; title('T2star (ms)');

%%
mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_cmr.mat'));
mask_roi = zeros(size(whatsinit{1},1), size(whatsinit{1},2));
mask_remote = zeros(size(whatsinit{1},1), size(whatsinit{1},2));
mask_cmr = struct;
if ~exist(mask_f)
    for i = 1:1
        % for i = 1:size(whatsinit{1},3)

        figure();
        imagesc(img_t2_cmr*0.1); axis image; caxis([0 200]);
        roi = drawpolygon;
        mask_roi(:,:,i) = createMask(roi);

        imagesc(img_t2_cmr*0.1); axis image; caxis([0 200]);
        remote = drawpolygon;
        mask_remote(:,:,i) = createMask(remote);
    end
    mask_cmr.mask_roi = mask_roi;
    mask_cmr.mask_remote = mask_remote;
    save(mask_f, '-struct', 'mask_cmr');
else
    load(mask_f);
end
%%
thresh = mean(nonzeros(img_t2star_cmr .* mask_remote)) - std(nonzeros(img_t2star_cmr .* mask_remote));
mask_hemo = (img_t2star_cmr < thresh) .* mask_roi;
figure(); imagesc(mask_hemo); %caxis([0 50]);

% t1_array_hemo = nonzeros(img_t1_molli .* mask_hemo);
% img_t2star_cmr(img_t2star_cmr == 0) = nan;
t2star_array_hemo = nonzeros(img_t2star_cmr .* mask_hemo);
img_t2star_cmr_label = img_t2star_cmr .* mask_hemo > 0;
t2_array_hemo = img_t2_cmr(img_t2star_cmr_label);

%t1_array_nonhemo = nonzeros(img_t1_molli .* (mask_roi -mask_hemo));
t2star_array_nonhemo = nonzeros(img_t2star_cmr .* (mask_roi -mask_hemo));
img_t2star_cmr_label2 = img_t2star_cmr .* (mask_roi -mask_hemo) > 0;
t2_array_nonhemo = img_t2_cmr(img_t2star_cmr_label2);

figure();
p = plot(t2_array_hemo*0.1, t2star_array_hemo*0.1); 
p.LineStyle = "none";
p.Color = "red";
p.Marker = "o";
hold on;
p1 = plot(t2_array_nonhemo*0.1, t2star_array_nonhemo*0.1); 
grid on;
p1.LineStyle = "none";
p1.Color = "blue";
p1.Marker = "x";
xlabel('T2 (ms)'); ylabel('T2star (ms)');