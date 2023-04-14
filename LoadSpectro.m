close all;
clear all;
%% load SVS Data
% Run startup.m in D:\\OXSA\
base_dir = uigetdir;
%base_dir = 'D:\Data\Cedars_2020_pig\DHARMAKUMAR_5257_20P10_05222020_WEEK4\';
% addpath(cat(2, base_dir, 'src'));
% addpath('D:\Data\Exvivo_Phantom'); % To add glob function...
% addpath('D:\src\function\');
addpath('function/');

dir_glob = glob(cat(2, base_dir, '\SVS_*'));

water_ppm = 4.7;
ppm_ref = water_ppm;


disp(dir_glob);
idx_array = input('Please refer to dir_glob for the index you want to import \n Enter cases with [ ] around them (Before | After):  '); % Before, After water suppression

InstanceNum = 1;
ss = cell(length(idx_array), 1);

%%
for i = 1:length(idx_array)
    dst_dir = dir_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    SeriesNum = str2num(f_parts{end});
    concent = f_parts{end-1};
    
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
%% Get the signal percentage
N = linspace(1,n,n);
perc = zeros(size(ss,1), 1);
cellToSave = cell(size(ss,1), 2);

for i = 1:size(ss, 1)
    dst_dir = dir_glob{idx_array(i)};
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

% save(GetFullPath(cat(2, base_dir, '\..\perc_DiffParameters.mat')), 'cellToSave');
%% Plot spectrum
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
    dst_dir = dir_glob{idx_array(i)};
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
    dst_dir = dir_glob{idx_array(i)};
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
    dst_dir = dir_glob{idx_array(i)};
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
%% Add median filter
[mmin, mm] = min(max(lim_mat, [], 2));
yl_cuttingwater(mx) = ceil(0.02*max(max(lim_mat, [], 2)));
yl_cuttingwater(mm) = ceil(0.1*min(max(lim_mat, [], 2)));
% Cutting water peak YLim
figure();
for i = 1:size(ss, 1)
    dst_dir = dir_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    concent = f_parts{end-1};
    
    spectrum = ss{i}.spectra{1};
    spectrum_filt = medfilt1(abs(spectrum), 30);
    
    subplot(1,2,i)
    plot(ppm_axis, spectrum_filt); xlabel('ppm'); ylabel('intensity');
    hold on;
    rectangle('Position', [ppm_axis(idx_low), 0, ppm_axis(idx_high)-ppm_axis(idx_low), max(abs(spectrum))], 'FaceColor', [1 0 0 0.1])
    %plot(ppm_axis, medfilt1(abs(spectrum), 10));
    H = area(ppm_axis(idx_low:idx_high), spectrum_filt(idx_low:idx_high));
    set(H(1), 'FaceColor', [0.8500    0.3250    0.0980])
    set(gca, 'Xdir', 'reverse')
    title(cat(2, 'MR Spectroscopy ', concent))
    ylim([0 yl_cuttingwater(i)])
    grid on;
end


%% Show moive while moving Ylim
clear movieVector


yl_frame = (ceil(mmin)+1):1:fixed_yl;
figure();
for i = 1:length(yl_frame)
    spectrum = ss{mm}.spectra{1};
    plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity');
    set(gca, 'YLim', [0 yl_frame(i)]);
    grid on;
    %hold on;
    % drawnow limitrate
    % pause(0.002)
    movieVector(i) = getframe;
end

fname = GetFullPath(cat(2, base_dir, '\..\Spectro_Ylim.avi'));
RecordVideo(movieVector, fname, 200);