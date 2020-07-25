clear all;
close all;
clc
current_dir = pwd;
%% 
addpath('..\function\');
base_dir = GetFullPath(cat(2, current_dir, '\..\..\T1_Fat_Project\Data\'));

name_glob = glob(cat(2, base_dir, '*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'T1', 'LGE', 'T2star', 'FF'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'freeROI'};
output_label = {'Mask'};

time_label = {'0D', '0D_occl', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};

for n = 1:length(Names)
    name = Names{n};
    for tp = 1:length(time_label)
        time_point = time_label{tp};
        clear orig_img
        for sql = 1:(length(sequence_label)-1)
            label = sequence_label{sql};
            
            f_glob = glob(cat(2, base_dir, name, '\', name, '_', time_point, '\', label, '\*'));
            
            
            if ~isempty(f_glob)
                [list_to_read, order_to_read] = NamePicker(f_glob);
                
                save_dir = GetFullPath(cat(2, base_dir, name, '\', name, '_', time_point,  '\', output_label{1}, '\'));
                if ~exist(save_dir, 'dir')
                    mkdir(save_dir)
                end
                
                clear whatsinit myo_mask
                
                if ~strcmp(label, 'T2star')
                    for i = 1:length(order_to_read)
                        f = list_to_read{order_to_read(i)};
                        whatsinit(:,:,i) = dicom23D(f);
                    end
                    orig_img.(label) = whatsinit;
                    
                    %% draw myocardium
                    f_myo = cat(2, save_dir, label, '_', anatomy_label{2}, '.mat');
                    if ~exist(f_myo, 'file')
                        myo_mask = zeros(size(whatsinit));
                        for i = 1:length(order_to_read)
                            disp('Epicardium: ');
                            figure('Position', [100 0 1600 1600]); imagesc(whatsinit(:,:,i)); axis image;
                            epi = drawpolygon(gca);
                            
                            disp('Endocardium: ');
                            figure('Position', [100 0 1600 1600]); imagesc(whatsinit(:,:,i)); axis image;
                            endo = drawpolygon(gca);
                            
                            epi_mask = createMask(epi);
                            endo_mask = createMask(endo);
                            
                            myo_mask(:,:,i) = epi_mask - endo_mask;
                        end
                        
                        close all;
                        save(f_myo, 'myo_mask');
                    else
                        disp(['It was done already', name, '_', time_point]);
                    end
                else
                    for i = 1:length(order_to_read)
                        f = list_to_read{order_to_read(i)};
                        whatsinit(:,:,:,i) = dicom23D(f);
                    end
                    orig_img.(label) = whatsinit;
                    
                     %% draw myocardium
                    f_myo = cat(2, save_dir, label, '_', anatomy_label{2}, '.mat');
                    if ~exist(f_myo, 'file')
                        myo_mask = zeros(size(whatsinit,1), size(whatsinit,2), size(whatsinit,4));
                        %myo_mask = cell(1, length(whatsinit));
                        for i = 1:length(order_to_read)
                            disp('Epicardium: ');
                            figure('Position', [100 0 1600 1600]); imagesc(whatsinit(:,:,2,i)); axis image;
                            epi = drawpolygon(gca);
                            
                            disp('Endocardium: ');
                            figure('Position', [100 0 1600 1600]); imagesc(whatsinit(:,:,2,i)); axis image;
                            endo = drawpolygon(gca);
                            
                            epi_mask = createMask(epi);
                            endo_mask = createMask(endo);
                            
                            myo_mask(:,:,i) = epi_mask - endo_mask;
                            %myo_mask{i} = epi_mask - endo_mask;
                        end
                        
                        close all;
                        save(f_myo, 'myo_mask');
                    else
                        disp(['It was done already', name, '_', time_point]);
                    end
                    
                end
            end
        end
        
        save(GetFullPath(cat(2, base_dir, name, '\', name, '_', time_point,  '\', 'orig_img.mat')), 'orig_img');
    end
end
