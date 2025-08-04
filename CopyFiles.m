% Define the source and destination folders
addpath('../function/');
source_folder = uigetdir;
%%
strings = strsplit(source_folder, '/');
names = strings{end-1};
strings2 = strsplit(names, '_');
name = strings2{1};
time_point = strings2{2};

switch time_point
    case {'D5', 'D6', 'D7', 'D8'}
        label = 'acute';
    case {'WK8', 'WK8+2'}
        label = 'chronic';
end


destination_folder = GetFullPath(cat(2, source_folder, '/../../../ROC_Analysis/DICOM/', name, '_', label, '/' ,time_point, '/LRT/T2_MULTIECHO/'));
%%
% Create the destination folder if it doesn't exist
if ~exist(destination_folder, 'dir')
    mkdir(destination_folder);
end

% Get a list of all files in the source folder
file_list = dir(source_folder);

% Loop through each file in the source folder
for i = 1:length(file_list)
    % Get the current file name
    file_name = file_list(i).name;
    
    % Check if "(2)" is in the file name
    if contains(file_name, '(1)')
        % Construct full file paths
        source_file_path = fullfile(source_folder, file_name);
        destination_file_path = fullfile(destination_folder, file_name);
        
        % Copy the file to the destination folder
        copyfile(source_file_path, destination_file_path);
        fprintf('Copied: %s\n', file_name);
    end
end

disp('All matching files have been copied.');