% check if downsize for LGE works ok

clear all;
close all;

addpath('../function/');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '/*'));

labels = {'PSIR'};
target_res = [1.2, 1.3, 1.4, 1.5, 1.6];
% Sofia D5
% [114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134]

% Paris D8
% [165, 167, 169, 171, 173, 175, 177, 179, 181, 183, 185, 187]

% Lisbon D8
% [164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184]

for ll = 1:length(labels)
    label = labels{ll};
    idx_array = contains(folder_glob, label);
    
    if any(idx_array)
        num = numel(nonzeros(idx_array));
        ind_array = find(idx_array == 1);
        dst_files = cell(num, 1);
        for i = 1:num
            dst_files{i} = folder_glob{ind_array(i)};
        end
        
        dst_name = ExtractNames(dst_files);
        disp(dst_name);
        
        sel_array = input('Please add an array here:  ');
        
        char_array = num2str(sel_array', '%04.f');
        
        ind_array2 = zeros(size(dst_name, 1), 1);
        for i = 1:size(char_array, 1)
            cha = char_array(i, :);
            new_cha = strip(cha,'left','0');
            if any(contains(dst_name, cha))
                ind_array2 = ind_array2 + contains(dst_name, cha);
            elseif any(contains(dst_name, new_cha))
                ind_array2 = ind_array2 + contains(dst_name, new_cha);
            end
        end
        
        ind_array3 = find(ind_array2 == 1);
        list_to_read = dst_files(ind_array3);
        
        name_to_compare = ExtractNames(list_to_read);
        
        order_to_read = zeros(length(list_to_read), 1);
        for i = 1:length(list_to_read)
            cha = char_array(i, :);
            new_cha = strip(cha,'left','0');
            if any(contains(dst_name, cha))
                order_to_read(i) = find(contains(name_to_compare, cha) == 1);
            elseif any(contains(dst_name, new_cha))
                order_to_read(i) = find(contains(name_to_compare, new_cha) == 1);
            end
        end
        
        save_dir = GetFullPath(cat(2, base_dir, '/../img_LEG_resize/'));
        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end
        
        %fname_list =  ExtractNames(list_to_read);
        for i = 1:length(list_to_read)
            f = list_to_read{order_to_read(i)};
            [whatsinit, slice_data, ~] = dicom23D(f);
            img_cropped = CropAroundHeart_NoCentroid(whatsinit);
            
            scale_array = slice_data.PixelSpacing(1) ./ target_res;
            
%             for n = 1:length(scale_array)
%                 re_img = imresize(img_cropped, scale_array(n), 'bilinear');
%                 figure();
%                 subplot(1,2,1);
%                 imagesc(img_cropped); axis image;
%                 subplot(1,2,2);
%                 imagesc(re_img); axis image;
%             end
            
            if size(img_cropped, 3)  == 1
                for n = 1:length(scale_array)
                    re_img = imresize(img_cropped, scale_array(n), 'bilinear');
                    re_img = re_img / max(re_img(:));
                    re_img = adapthisteq(re_img);
                    f_to_save = cat(2, save_dir, label, '_SAX', num2str(i, '%02.f'), '_resize', num2str(target_res(n), '%.2f'),  '.png');
                    imwrite(mat2gray(re_img), f_to_save);
                end
                img_cropped = img_cropped / max(img_cropped(:));
                img_cropped = adapthisteq(img_cropped);
                %figure();
                %imagesc(img_cropped); colormap gray; axis image;
                f_to_save = cat(2, save_dir, label, '_SAX', num2str(i, '%02.f'), '_orig', num2str(slice_data.PixelSpacing(1), '%.2f'), '.png');
                imwrite(mat2gray(img_cropped), f_to_save);
                
            else
                
                disp("Not Yet Implemented");
                
            end
        end
        
    else
        disp('No Matching labels found');
        disp(label);
    end
end