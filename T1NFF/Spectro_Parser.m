close all;
clear all;
% load SVS Data
% This doesn't work for MATLAB 2020a
% But works for 2017b
base_dir = uigetdir;
%addpath(cat(2, base_dir, 'src'));
addpath('../function/'); % To add glob function...
%oxsa_dir = 'D:\OXSA\main\';
%addpath(oxsa_dir);

dir_glob = glob(cat(2, base_dir, '/SVS_*'));

water_ppm = 4.6;
ppm_ref = water_ppm;

idx_array = [1, 2, 3, 4]; % Water Suppression


InstanceNum = 1;
ss = cell(length(idx_array), 1);

for i = 1:length(idx_array)
    dst_dir = dir_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    SeriesNum = str2num(f_parts{end});
    concent = f_parts{end-1};
    
    ss{i} = loadSVSData(dst_dir, SeriesNum, InstanceNum);
    
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
%% Plot
figure();
for i = 1:size(ss, 1)
    dst_dir = dir_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    concent = f_parts{end-1};
    
    spectrum = ss{i}.spectra{1};
    subplot(2,2,i)
    
    plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity')
    hold on;
    rectangle('Position', [ppm_axis(idx_low), 0, ppm_axis(idx_high)-ppm_axis(idx_low), max(abs(spectrum))], 'FaceColor', [1 0 0 0.1])
    H = area(ppm_axis(idx_low:idx_high), abs(spectrum(idx_low:idx_high)));
    set(H(1), 'FaceColor', [0.8500    0.3250    0.0980])
    set(gca, 'Xdir', 'reverse')
    title(cat(2, 'MR Spectroscopy ', concent))
    grid on;
end