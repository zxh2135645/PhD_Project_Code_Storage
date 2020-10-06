% Cropping cvi image with no centroids
clear all;
close all;

addpath('D:\src\function');
base_dir = 'D:\T1_Fat_Project\';
folder_glob = glob(cat(2, base_dir, 'Data\*'));

Names = ExtractNames(folder_glob);

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

time_points = {'0D', '0D_occl', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};
out_dir = cat(2, base_dir, 'cropped_cvi\');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%%

for n = 1:length(Names)
    name = Names{n};
    
    for t = 1:length(time_points)
        tp = time_points{t};
        cvi_glob = glob(cat(2, base_dir, 'Data\', name, '\', name, '_', tp, '\Screenshot\*'));
        for i = 1:length(cvi_glob)
            I = imread(cvi_glob{i});
            img_size = size(I);
            x_dis = round(img_size(2)/6);
            y_dis = round(img_size(1)/6);
            img_cropped = imcrop(I, [x_dis, y_dis, 4*x_dis, 4*y_dis]);
            strings = strsplit(cvi_glob{i},'\');
            fname = strings{end};
            strings = strsplit(fname, '-');
            fname_new = cat(2, name, '_', tp, '-', strings{2}, '-', strings{3}, '-', strings{4});
            f_save = cat(2, out_dir, fname_new, '.png');
            imwrite(img_cropped, f_save);
            
            %copyfile(cvi_glob{i}, out_dir);
        end
    end
    
end