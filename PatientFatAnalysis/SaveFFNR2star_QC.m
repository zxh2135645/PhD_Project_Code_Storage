clear all;
close all;
clc;
current_dir = pwd;
% Patient data configuration for Khalid
%% 
addpath('../function/');
addpath('../AHA16Segment/');
base_dir = uigetdir;

contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = cell(length(contour_glob), 1); 
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'/');
    Names{i} = strings{end-1};
end

sequence_label = {'LGE', 'T2star'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

name_check = {'484060000001'};
name_idx_list = linspace(1, length(Names), length(Names)); % initialize with incremental add

if length(name_check) == 1
    starting_point = find(strcmp(name_check, Names),1);
else
    name_idx_list = zeros(1, length(name_check));
    for n = 1:length(name_check)
        % Check an array of names
        name_idxo = find(strcmp(name_check(n), Names),1);
        name_idx_list(n) = name_idxo;
    end
end

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
    };


output_label = {'LGE', 'T2star'};
save_dir = GetFullPath(cat(2, base_dir, '/Analysis/'));
data_save_dir = cat(2, base_dir, '/data/');

time_points = {'BL', 'BL2', 'FU'};

label_lge = sequence_label{1};
label_t2star = sequence_label{2};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Read excel file
T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'));
id_array = T.AnonymizationID;
hemo_array = T.withOrWithout;


%% Array initialization
name_label = {};
slice_count = 1;
vec = @(x) x(:);

%for n = 1:(length(Names)-1)
for n = 28:28
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    name_data_save_dir = cat(2, data_save_dir, name);
    if ~exist(name_data_save_dir, 'dir')
        mkdir(name_data_save_dir);
    end
    % tp_count = 0;
    
    name_for_table_searching = insertAfter(name, 6, '-');
    row = find(contains(id_array,name_for_table_searching));
    
    IMH_cell = table2cell(T(row, 13)); % IMH
    IMH = IMH_cell{1};
    for tp = 1:length(time_points)
    %for tp = 1:3
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            % tp_count = tp_count+1;
            
            % FF
            myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
            
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat')); % glob_names
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat')); % slc_array
            slc_array_t2star = slc_array;
            
            ff_map = cell(1, length(glob_names));
            num_array = ExtractNum(glob_names);
            for f = 1:length(ff_map)
                % series_name = glob_names{f};
                ff_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(num_array(f)), '.mat'));
                ff_map{f} = load(ff_glob{1}, 'fwmc_ff');
            end
            
            % convert ff_map to matrix
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            for f = 1:length(ff_map)
                ff(:,:,f) = ff_map{f}.fwmc_ff;
            end

            % No need to reorder for analysis
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(num_array(f)), '.mat'));
                r2star_map{f} = load(r2star_glob{1}, 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
            for f = 1:length(r2star_map)
                r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
            end
            
            
            tp_dir2 = cat(2, name_save_dir, '/', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            
            figure();
            for slc = 1:size(ff,3)
                subplot(1,2,1);
                imagesc(ff(:,:,slc)); axis image; caxis([0 50]); title('FF map'); axis off; 
                subplot(1,2,2);
                imagesc(r2star(:,:,slc)); axis image; caxis([0 100]); title('R2star map'); axis off;
                saveas(gcf, cat(2, tp_dir2, 'FFNR2star_slc', num2str(slc) ,'.png'))
            end

        end
    end
end