close all;
clear all;
% load SVS Data
% This doesn't work for MATLAB 2020a
% But works for 2017b
base_dir = uigetdir;
%addpath(cat(2, base_dir, 'src'));
addpath('./function/'); % To add glob function...
%oxsa_dir = 'D:\OXSA\main\';
%addpath(oxsa_dir);

dir_glob = glob(cat(2, base_dir, '/svs_*'));

%%
addpath('../OXSA/');
water_ppm = 4.7;
ppm_ref = water_ppm;
% idx_array = [1, 2, 3, 4]; % Water Suppression
idx_array = [2]; % Water Suppression
InstanceNum = 1;

ss = cell(length(idx_array), 1);
for i = 1:length(idx_array)
    dst_dir = dir_glob{idx_array(i)};
    f_parts = strsplit(fileparts(dst_dir), '_');
    SeriesNum = str2num(f_parts{end});
    concent = f_parts{end-1};
    
    % ss{i} = loadSVSData(dst_dir, SeriesNum, InstanceNum);
    dt = Spectro.dicomTree('dir',dst_dir);
    matched = dt.searchForSeriesInstanceNumber(SeriesNum, InstanceNum);
    ss{i} = Spectro.dicom(matched);

    if i == 1
        s = ss{i};
        info = s.info{1};
        Num_avg = info.SharedFunctionalGroupsSequence.Item_1.MRAveragesSequence.Item_1.NumberOfAverages;
        dwt = 1/info.Private_0021_1401;
        % dwt = s.dwellTime; % in the unit of seconds

        n = size(info.SpectroscopyData, 1);
        % n = s.samples;  % number of samples

        tmp = size(info.SpectroscopyData, 1)/2;
        w = ((-tmp):(tmp-1)).'/(dwt*tmp*2);
        % w = s.freqAxis; % frequency axis (Hz)
        
        t = dwt*(0:n-1).';
        % t = s.timeAxis; % time axis (s)

        imagingFrequency = info.TransmitterFrequency;
        ppm_axis = w / imagingFrequency;
        % ppm_axis = s.ppmAxis + ppm_ref;
    end
    spectrum = ss{i}.info{1}.SpectroscopyData;
end

%% Plot
figure();
subplot(3,1,1)
dst_dir = dir_glob{idx_array(i)};
f_parts = strsplit(fileparts(dst_dir), '_');
concent = f_parts{end-1};
% plot(ppm_axis, abs(spectrum)); xlabel('ppm'); ylabel('intensity')
plot(abs(spectrum));  ylabel('intensity'); %xlabel('ppm');
title('Original Data');

subplot(3,1,2);
plot(spectrum(1:2:end)); ylabel('intensity'); %xlabel('ppm');
title('Odd Number');

subplot(3,1,3);
plot(spectrum(2:2:end)); ylabel('intensity'); %xlabel('ppm');
title('Even Number');

