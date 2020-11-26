clear all;
close all;
%% Delete duplicate files
% base directory
base_dir = uigetdir;
%% 
dir_glob = glob(cat(2, base_dir, '/*'));

for dd = 1:length(dir_glob)
    D = dir(dir_glob{dd});
    
    for i = 1:length(D)
        if ~strcmp(D(i).name, '.') || strcmp(D(i).name, '..')
            strings = strsplit(D(i).date, ' ');
            times = strings{end};
            hhmmss = strsplit(times, ':');
            hh = hhmmss{1};
            if strcmp(hh, '08')
                folder = D(i).folder;
                filename = D(i).name;
                delete(cat(2, folder, '/', filename));
            end
        end
    end
end