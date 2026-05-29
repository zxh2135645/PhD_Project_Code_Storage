function dicomData = readAllDicom(folderPath)
% READALLDICOM Reads all DICOM files (.dcm or .IMA) in a folder recursively.
%
%   dicomData = READALLDICOM(folderPath)
%   Input:
%       folderPath - Path to the folder containing DICOM files
%
%   Output:
%       dicomData  - Struct array with fields:
%                       .Filename   : Full file path
%                       .Info       : DICOM metadata (from dicominfo)
%                       .Image      : DICOM image data (from dicomread)
%
%   Example:
%       data = readAllDicom('/Users/jameszhang/Documents/MRI/Scan01');
%       imshow(data(1).Image, []);

    % Validate input
    if nargin < 1 || ~isfolder(folderPath)
        error('Input must be a valid folder path.');
    end

    % Find all DICOM files recursively
    fileList = dir(fullfile(folderPath, '**', '*'));
    dicomFiles = fileList(~[fileList.isdir]); % exclude folders
    dicomFiles = dicomFiles(endsWith({dicomFiles.name}, {'.dcm', '.IMA'}, 'IgnoreCase', true));

    if isempty(dicomFiles)
        warning('No DICOM files found in the specified folder.');
        dicomData = [];
        return;
    end

    % Preallocate struct array
    dicomData = struct('Filename', {}, 'Info', {}, 'Image', {});

    % Read DICOM data
    for i = 1:numel(dicomFiles)
        fname = fullfile(dicomFiles(i).folder, dicomFiles(i).name);
        try
            info = dicominfo(fname);
            img = dicomread(info);
        catch ME
            warning('Skipping file %s: %s', fname, ME.message);
            continue;
        end
        dicomData(end+1) = struct('Filename', fname, 'Info', info, 'Image', img); %#ok<AGROW>
    end

    fprintf('Successfully loaded %d DICOM files from %s\n', numel(dicomData), folderPath);
end