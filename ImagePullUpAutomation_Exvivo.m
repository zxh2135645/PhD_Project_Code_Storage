clear all;
close all;
% Image Pullup automation for exvivo scan which I'm expecting 3D acquistion
% for images

%% 20P10, 20P11 Pig scan

addpath('D:\src\function');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
% T1 is IR T1 Map
% T1MAP is T1 MOLLI
% T2 is T2 anatomical map
% T2MAP is T2 transversal coded
% T2STAR is longer TE/higher resolution T2STAR map
% _T2STAR is for invivo t2* map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P10 Exvivo
% T1_SIEMENS_MAPIT
% [19]
% T1MAP
% [24]
% T2_SIEMENS_MAPIT
% [35]
% T2MAP
% [27]
% T2STAR_SIEMENS_MAPIT
% [21]
% _T2STAR
% [29]
% mGRE
% [8]
% T1 TSE
% [9]
% FATSAT
% [11]
% MT
% [10]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P09 Exvivo
% T1_SIEMENS_MAPIT
% [51]
% T1MAP
% [21]
% T2_SIEMENS_MAPIT
% [16]
% T2MAP
% [24]
% T2STAR_SIEMENS_MAPIT
% [18]
% _T2STAR
% [26]
% mGRE
% [8]
% T1 TSE
% [9]
% FATSAT
% [11]
% MT
% [10]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20P11 Exvivo
% T1_SIEMENS_MAPIT
% [14]
% T1MAP
% [40]
% T2_SIEMENS_MAPIT
% [16]
% T2MAP
% [43]
% T2STAR_SIEMENS_MAPIT
% [18]
% _T2STAR
% [45]
% mGRE
% [8]
% T1 TSE
% [9]
% FATSAT
% [11]
% MT
% [10]
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
            ind_array2 = ind_array2 + contains(dst_name, cha);
        end
        
        ind_array3 = find(ind_array2 == 1);
        list_to_read = dst_files(ind_array3);
        
        name_to_compare = ExtractNames(list_to_read);
        
        order_to_read = zeros(length(list_to_read), 1);
        for i = 1:length(list_to_read)
            order_to_read(i) = find(contains(name_to_compare, char_array(i, :)) == 1);
        end
        
        save_dir = GetFullPath(cat(2, base_dir, '\..\img\'));
        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end
        
        %fname_list =  ExtractNames(list_to_read);
        for i = 1:length(list_to_read)
            f = list_to_read{order_to_read(i)};
            whatsinit = dicom23D(f);
            if strcmp(label, 'MGRE')
                whatsinit = reshape(whatsinit, size(whatsinit, 1), size(whatsinit, 2), 8, []); % 8 TE in mGRE 
            end
            clear img_cropped img_cropped_2d
            
            if size(whatsinit,3) == 1
                img_cropped = CropAroundHeart_NoCentroid(whatsinit);
                img_cropped = img_cropped / max(img_cropped(:));
                img_cropped = adapthisteq(img_cropped);
                %fname = fname_list{order_to_read(i)};
                
                %figure();
                %imagesc(img_cropped); colormap gray; axis image;
                
                f_to_save = cat(2, save_dir, label, '_SAX', num2str(i, '%02.f'), '.png');
                imwrite(mat2gray(img_cropped), f_to_save);
                
            elseif size(whatsinit, 4) == 1
                img_cropped = CropAroundHeart_NoCentroid(whatsinit);
                for slc = 1:size(img_cropped, 3)
                    img_cropped_2d = img_cropped(:,:,slc);
                    img_cropped_2d = img_cropped_2d / max(img_cropped_2d(:));
                    img_cropped_2d = adapthisteq(img_cropped_2d);
                    f_to_save = cat(2, save_dir, label, '_SAX', num2str(slc, '%02.f'), '.png');
                    imwrite(mat2gray(img_cropped_2d), f_to_save);
                end
            else
                % slc_idx = round(3/8 * size(whatsinit, 4));
                slc_idx = round(6/8 * size(whatsinit, 4));
                for tp = 1:size(whatsinit, 3)
                    img_2D = whatsinit(:,:,tp, slc_idx);
                    img_cropped = CropAroundHeart_NoCentroid(img_2D);
                    img_cropped = img_cropped / max(img_cropped(:));
                    img_cropped = adapthisteq(img_cropped);
                    f_to_save = cat(2, save_dir, label, '_SLC', num2str(slc_idx, '%02.f'),  '_TP', num2str(tp, '%02.f'), '.png');
                    imwrite(mat2gray(img_cropped), f_to_save);
                end
                

            end
        end
        
    else
        disp('No Matching labels found');
        disp(label);
    end
end


