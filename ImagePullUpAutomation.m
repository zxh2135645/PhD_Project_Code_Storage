clear all;
close all;
% Image Pullup automation
% Can not read all files in
% Is it worthwhile?

% I think it's necessary, let me work on it 05/27/2020
% Goal is to automate globing and export gray colormap images
% It doesn't take care of manual windowing
%% 20P10, 20P11 Pig scan

addpath('D:\src\function');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

% labels = {'SHMOLLI', 'T1MAP', 'MAG', 'PSIR', 'T2STAR', 'T2MAP'};
% labels for 20PXX
labels = {'T1MAP', 'MAG', 'PSIR', 'T2STAR', 'T2MAP', 'T2Mapping', 'T2_MAP'};
% % T2* weighted image
% labels for Zhengzhou
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P10
% SHMOLLI
% [96, 100, 103, 106, 109, 112, 115, 118, 142, 145, 127, 130, 133, 136, 139]
% T1MAP
% [42, 45, 51, 54, 57, 60, 63, 66, 69, 72, 75, 304, 81, 84, 87]
% MAG
% [320, 322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 350]
% PSIR
% [321, 323, 325, 327, 329, 331, 333, 335, 337, 339, 341, 343, 345, 347, 349, 351]
% T2STAR
% [211, 213, 215, 217, 219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P11
% SHMOLLI
% [88, 94, 97, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 131]
% T1MAP
% [40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79]
% MAG
% [236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 256, 258, 260, 262, 264]
% PSIR
% [237, 239, 241, 243, 245, 247, 249, 251, 253, 255, 257, 259, 261, 263, 265]
% T2STAR
% [185, 187, 189, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P10 wk8
% SHMOLLI
% [104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146]
% T1MAP
% [43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88]
% T2STAR
% [219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247]
% T2MAP
% [206, 209, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194, 197, 200, 203]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P09 wk8
% SHMOLLI
% [100, 103, 106, 109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139]
% T1MAP
% [45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84]
% T2STAR
% [203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 227, 229]
% T2MAP
% [154, 157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P11 wk8
% SHMOLLI
% [101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143]
% T1MAP
% [43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85]
% T2STAR
% [210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238]
% T2MAP
% [158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194, 197, 200]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhengzhou Qin Tie Gang
% MAG
% [48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 72, 74, 76, 78, 80, 82, 84, 86, 88]
% PSIR
% [49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 73, 75, 77, 79, 81, 83, 85, 87, 89]
% T2STAR
% [43]
% T2Mapping
% [39, 40, 41, 42]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = {'T1MAP', 'T2STAR', 'T2_MAP'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMRP resolution study - 8-week data
% 20P40 8wk
% T1MOLLI
% [46, 49, 52]
% MAG
% [247 249 251]
% PSIR
% [248 250 252]
% T2star
% [194 196 198 200]
% T2Mapping
% [149 152 155 158]
% T2_MAP weighted image
% [197]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P48 8wk
% T1MOLLI
% [46, 49, 52]
% MAG
% [247 249 251]
% PSIR
% [248 250 252]
% T2star
% [194 196 198 200]
% T2Mapping
% [149 152 155 158]
% T2_MAP weighted image
% [197]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        save_dir = GetFullPath(cat(2, base_dir, '\..\img\'));
        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end
        
        %fname_list =  ExtractNames(list_to_read);
        for i = 1:length(list_to_read)
            f = list_to_read{order_to_read(i)};
            whatsinit = dicom23D(f);
            img_cropped = CropAroundHeart_NoCentroid(whatsinit);
            
            if size(img_cropped, 3)  == 1
                img_cropped = img_cropped / max(img_cropped(:));
                img_cropped = adapthisteq(img_cropped);
                %fname = fname_list{order_to_read(i)};
                
                %figure();
                %imagesc(img_cropped); colormap gray; axis image;
                
                f_to_save = cat(2, save_dir, label, '_SAX', num2str(i, '%02.f'), '.png');
                imwrite(mat2gray(img_cropped), f_to_save);
            else
                for slc = 1:size(img_cropped, 3)
                    img_cropped(:,:,slc) = img_cropped(:,:,slc) / max(max(img_cropped(:,:,slc)));
                    img_cropped(:,:,slc) = adapthisteq(img_cropped(:,:,slc));
                    
                    if strcmp(label, 'T2Mapping') || strcmp(label, 'T2_MAP')
                        f_to_save = cat(2, save_dir, label, '_SAX', num2str(i, '%02.f'), '_TE', num2str(slc, '%02.f'), '.png');
                        imwrite(mat2gray(img_cropped(:,:,slc)), f_to_save);
                    else
                        f_to_save = cat(2, save_dir, label, '_SAX', num2str(slc, '%02.f'), '.png');
                        imwrite(mat2gray(img_cropped(:,:,slc)), f_to_save);
                    end
                end
            end
        end
        
    else
        disp('No Matching labels found');
        disp(label);
    end
end


