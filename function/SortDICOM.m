%%% Sort "lumped" DICOM files into folders for each scan
% Modified by XZ, 06/25/2022
% Use 'SeriesNumber' to get the ordering of scans and sort out cases
% 'SeriesDecription' are the same.
% Save MAG and Phase in different folder if there is
% Minor changes to the waitbar

close all;
clear all;
clc;

% base_dir = 'C:\Users\ganthon\OneDrive - Indiana University\FromBox\MATLAB\T2starMappingMyocardium\GregCode';
% cd(base_dir)
addpath('../function/');
% addpath('../mfiles/');
% rmpath('/Users/jameszhang/Documents/MATLAB/BI265/mfiles/');
% This line needs to be modified

% Read files
data_dir = uigetdir;
% data_name = data_dir(strfind(data_dir, '/Data/')+6:end);
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
    'EchoTime',...
    };

listing = dir(data_dir);

N = numel(listing); % How many entries in the directory listing
    if (N<3)
        error('Empty folder');
        return
    end
    
h = waitbar(0,'Reading DICOM Files...','WindowStyle','modal');

true_index = 0; % a sequential index of dicom files, that is ignoring
% files of other types.

steps = length(listing) - 2;
for i = 3:length(listing) % loop through directory listing, but skip '.' and '..'
    isPhase = '';

    waitbar((i-2)/steps,h,'Reading DICOM Files...');
    filename = listing(i).name;
    [dummy_path, just_the_name, extension] = fileparts(filename);
    full_path = fullfile(data_dir, filename);

    goodfile = false;

    % Check for good dicom file
    if isdicom(full_path)
        true_index = true_index + 1;
        header = dicominfo(full_path);
        if isfield(header, 'SeriesDescription')
            scan_name = header.SeriesDescription;
            series_number = num2str(header.SeriesNumber);
            if length(series_number) == 4
                series_name = num2str(str2num(series_number(1)), '%.4d');
            elseif length(series_number) == 6
                series_name = num2str(str2num(series_number(1:3)), '%.4d');
            else
                series_name = num2str(str2num(series_number(1:2)), '%.4d');
            end

            image_type_name = strsplit(header.ImageType, '\');
            isPhase = image_type_name{3};

            if strcmp(isPhase, 'P')
                series_dir = cat(2, data_dir, '/', scan_name, '_PHASE_', series_name);
            else
                series_dir = cat(2, data_dir, '/', scan_name, '_', series_name);
            end
            if ~exist(series_dir, 'dir')
                mkdir(series_dir);
            end
            movefile(full_path, series_dir)

        end
    end
end

%%%
close(h);