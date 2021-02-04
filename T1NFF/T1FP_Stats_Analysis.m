close all;
clear all;
% Statistical analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% data_storage_rim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% p_storage
% T_t1
% T_ff
% T_r2star
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../function/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn'};
%time_points = {'0D_baseline','1D', '7D', '28D', '8WK', '6MO', '9MO', '1YR', '15YR'};
time_points = {'8WK', '12WK', '14WK' '6MO', '9MO', '1YR', '15YR'};

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

sequence_label = {'T1', 'T2star', 'LGE'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
load(cat(2, metrics_save_dir, 'data_storage_rim.mat'));

%% Statistically t-test
p_storage = struct;
for n = 1:length(Names)
    name = Names{n};
    p_storage(n).Name = name;
    p_storage(n).nolong = struct;
    
    tp_count = 1;
    for tp = 1:length(time_points)
        time_point = time_points{tp};
        
        
        if ~isempty(data_storage_rim(n).data(end-tp+1).time_point)
            p_storage(n).nolong(tp_count).time_point = time_point;
            remote_r2star_array = data_storage_rim(n).data(end-tp+1).remote_r2star_array;
            roi_r2star_array = data_storage_rim(n).data(end-tp+1).roi_r2star_array;
            remote_r2star_array(remote_r2star_array < 0) = 0;
            roi_r2star_array(roi_r2star_array < 0) = 0;
            [p_storage(n).nolong(tp_count).ff,h,stats] = ranksum(data_storage_rim(n).data(end-tp+1).roi_ff_array, ...
                data_storage_rim(n).data(end-tp+1).remote_ff_array);
            [h, p_storage(n).nolong(tp_count).t1] = ttest2(data_storage_rim(n).data(end-tp+1).roi_t1_array, ...
                data_storage_rim(n).data(end-tp+1).remote_t1_array);
            [h, p_storage(n).nolong(tp_count).r2star] = ttest2(roi_r2star_array, ...
                remote_r2star_array);
            tp_count = tp_count + 1;
        end
        
    end
end

%% ANOVA
for n = 1:length(Names)
    name = Names{n};
    p_storage(n).long = struct;
    
    tp_count = 1;
    strength_roi_ff = [];
    alloy_roi_ff = {};
    strength_roi_t1 = [];
    alloy_roi_t1 = {};
    strength_roi_r2star = [];
    alloy_roi_r2star = {};
    
    strength_remote_ff = [];
    alloy_remote_ff = {};
    strength_remote_t1 = [];
    alloy_remote_t1 = {};
    strength_remote_r2star = [];
    alloy_remote_r2star = {};
    
    for tp = 1:length(time_points)
        time_point = time_points{tp};
        
        
        if ~isempty(data_storage_rim(n).data(end-tp+1).time_point)
            p_storage(n).nolong(tp_count).time_point = time_point;
            remote_r2star_array = data_storage_rim(n).data(end-tp+1).remote_r2star_array;
            roi_r2star_array = data_storage_rim(n).data(end-tp+1).roi_r2star_array;
            remote_r2star_array(remote_r2star_array < 0) = 0;
            roi_r2star_array(roi_r2star_array < 0) = 0;
            
            strength_roi_ff = [strength_roi_ff, data_storage_rim(n).data(end-tp+1).roi_ff_array];
            alloy_roi_ff = [alloy_roi_ff, repmat({time_point},1,length(data_storage_rim(n).data(end-tp+1).roi_ff_array))];
            strength_remote_ff = [strength_remote_ff, data_storage_rim(n).data(end-tp+1).remote_ff_array];
            alloy_remote_ff = [alloy_remote_ff, repmat({time_point},1,length(data_storage_rim(n).data(end-tp+1).remote_ff_array))];
            
            strength_roi_t1 = [strength_roi_t1, data_storage_rim(n).data(end-tp+1).roi_t1_array];
            alloy_roi_t1 = [alloy_roi_t1, repmat({time_point},1,length(data_storage_rim(n).data(end-tp+1).roi_t1_array))];
            strength_remote_t1 = [strength_remote_t1, data_storage_rim(n).data(end-tp+1).remote_t1_array];
            alloy_remote_t1 = [alloy_remote_t1, repmat({time_point},1,length(data_storage_rim(n).data(end-tp+1).remote_t1_array))];
            
            strength_roi_r2star = [strength_roi_r2star, roi_r2star_array];
            alloy_roi_r2star = [alloy_roi_r2star, repmat({time_point},1,length(roi_r2star_array))];
            strength_remote_r2star = [strength_remote_r2star,  remote_r2star_array];
            alloy_remote_r2star = [alloy_remote_r2star, repmat({time_point},1,length(remote_r2star_array))];
            
            tp_count = tp_count + 1;
        end
        
    end
    
    [p_storage(n).long.roi_ff,tbl] = anova1(strength_roi_ff, alloy_roi_ff,'off');
    [p_storage(n).long.roi_t1,tbl] = anova1(strength_roi_t1, alloy_roi_t1,'off');
    [p_storage(n).long.roi_r2star,tbl] = anova1(strength_roi_r2star, alloy_roi_r2star,'off');
    
    [p_storage(n).long.remote_ff,tbl] = anova1(strength_remote_ff, alloy_remote_ff,'off');
    [p_storage(n).long.remote_t1,tbl] = anova1(strength_remote_t1, alloy_remote_t1,'off');
    [p_storage(n).long.remote_r2star,tbl] = anova1(strength_remote_r2star, alloy_remote_r2star,'off');
end

save(cat(2, metrics_save_dir, 'p_storage.mat'), 'p_storage');

%% Make tables showing mean and sd
load(cat(2, metrics_save_dir, 'data_storage_rim.mat'));
col1 = {'Merry_MI', 'Merry_Remote', 'Ryn_MI', 'Ryn_Remote', 'Mojave_MI', 'Mojave_Remote', 'Sahara_MI', 'Sahara_Remote', ...
    'ZZ_MI', 'ZZ_Remote',  'Tina_MI', 'Tina_Remote', 'Sunny_MI', 'Sunny_Remote' 'Queenie_MI', 'Queenie_Remote', ...
    'Hope_MI', 'Hope_Remote', 'Gobi_MI', 'Gobi_Remote', 'Felicity_MI', 'Felicity_Remote', 'Evelyn_MI', 'Evelyn_Remote'}';

% time_point = time_points';
wk8_t1 = cell(length(col1), 1);
wk12_t1 = cell(length(col1), 1);
wk14_t1 = cell(length(col1), 1);
mo6_t1 = cell(length(col1), 1);
mo9_t1 = cell(length(col1), 1);
yr1_t1 = cell(length(col1), 1);
yr15_t1 = cell(length(col1), 1);

wk8_ff = cell(length(col1), 1);
wk12_ff = cell(length(col1), 1);
wk14_ff = cell(length(col1), 1);
mo6_ff = cell(length(col1), 1);
mo9_ff = cell(length(col1), 1);
yr1_ff = cell(length(col1), 1);
yr15_ff = cell(length(col1), 1);

wk8_r2star = cell(length(col1), 1);
wk12_r2star = cell(length(col1), 1);
wk14_r2star = cell(length(col1), 1);
mo6_r2star = cell(length(col1), 1);
mo9_r2star = cell(length(col1), 1);
yr1_r2star = cell(length(col1), 1);
yr15_r2star = cell(length(col1), 1);

T_t1 = table(col1, wk8_t1, wk12_t1, wk14_t1, mo6_t1, mo9_t1, yr1_t1, yr15_t1);
T_ff = table(col1, wk8_ff, wk12_ff, wk14_ff, mo6_ff, mo9_ff, yr1_ff, yr15_ff);
T_r2star = table(col1, wk8_r2star, wk12_r2star, wk14_r2star, mo6_r2star, mo9_r2star, yr1_r2star, yr15_r2star);

for i = 1:length(data_storage_rim)
    for tp = 1:length(time_points)
        time_point = time_points{tp};
        if ~isempty(data_storage_rim(i).data(end-tp+1).time_point)
            mean_roi_t1 = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_roi_t1,1));
            sd_roi_t1 = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_roi_t1,1));
            mean_remote_t1 = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_remote_t1,1));
            sd_remote_t1 = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_remote_t1,1));
            
            mean_roi_ff = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_roi_ff,1));
            sd_roi_ff = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_roi_ff,1));
            mean_remote_ff = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_remote_ff,1));
            sd_remote_ff = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_remote_ff,1));
            
            mean_roi_r2star = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_roi_r2star,1));
            sd_roi_r2star = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_roi_r2star,1));
            mean_remote_r2star = num2str(round(for_analysis_rim(i).metrics(end-tp+1).mean_remote_r2star,1));
            sd_remote_r2star = num2str(round(for_analysis_rim(i).metrics(end-tp+1).sd_remote_r2star,1));
            
            T_t1(2*(i-1)+1, tp+1) = {cat(2, mean_roi_t1, ' +/- ', sd_roi_t1)};
            T_t1(2*(i-1)+2, tp+1) = {cat(2, mean_remote_t1, ' +/- ', sd_remote_t1)};
            T_ff(2*(i-1)+1, tp+1) = {cat(2, mean_roi_ff, ' +/- ', sd_roi_ff)};
            T_ff(2*(i-1)+2, tp+1) = {cat(2, mean_remote_ff, ' +/- ', sd_remote_ff)};
            T_r2star(2*(i-1)+1, tp+1) = {cat(2, mean_roi_r2star, ' +/- ', sd_roi_r2star)};
            T_r2star(2*(i-1)+2, tp+1) = {cat(2, mean_remote_r2star, ' +/- ', sd_remote_r2star)};
        else
            T_t1(2*(i-1)+1, tp+1) = {' '};
            T_t1(2*(i-1)+2, tp+1) = {' '};
            T_ff(2*(i-1)+1, tp+1) = {' '};
            T_ff(2*(i-1)+2, tp+1) = {' '};
            T_r2star(2*(i-1)+1, tp+1) = {' '};
            T_r2star(2*(i-1)+2, tp+1) = {' '};
        end
    end
end

T_t1.Properties.VariableNames = {'Subjects', '8WK', '12WK', '14WK' '6MO', '9MO', '1YR', '15YR'};
T_ff.Properties.VariableNames = {'Subjects', '8WK', '12WK', '14WK' '6MO', '9MO', '1YR', '15YR'};
T_r2star.Properties.VariableNames = {'Subjects', '8WK', '12WK', '14WK' '6MO', '9MO', '1YR', '15YR'};

save(cat(2, metrics_save_dir, 'Table_t1.mat'), 'T_t1');
save(cat(2, metrics_save_dir, 'Table_ff.mat'), 'T_ff');
save(cat(2, metrics_save_dir, 'Table_r2star.mat'), 'T_r2star');