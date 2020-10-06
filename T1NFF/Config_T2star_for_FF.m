% Copy and paste T2star weighted images acquired by mGRE sequence for
% generating fat fraction map
clear all;
close all;

addpath('../function/');
base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, 'Data/*'));
Names = ExtractNames(folder_glob);

%labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
time_points = {'0D_baseline', '0D_occl', '0D', '1D', '2D', '3D', '5D', '7D', '14D', '21D', '28D', '6WK', '8WK', '12WK', '14WK', '6MO', '9MO', '1YR', '15YR'};

OutputPath = GetFullPath(cat(2, base_dir, 'Data_for_FF/'));
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end


for i = 1:length(Names)
    name = Names{i};
    for t = 1:length(time_points)
        tp = time_points{t};
        f_to_copy = cat(2, base_dir, name, '/', name, '_', tp, '/T2star/');
        if exist(f_to_copy, 'dir')
            f_to_out = cat(2, OutputPath, name, '/', name, '_', tp, '/T2star/');
            copyfile(f_to_copy, f_to_out); 
        end
    end
end

