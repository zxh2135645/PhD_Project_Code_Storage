clear all;
close all;
addpath('../function/');
addpath('../AHA16Segment/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));

Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};

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

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');


label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};

%% Main body
data_for_R = struct;
%for n = 1:length(Names)

for n = 1:1
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    
    name = Names{n};
    data_for_R(n).Name = name;
    data_for_R(n).metrics = struct;
    
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    name_data_save_dir = cat(2, data_save_dir, name);
    if ~exist(name_data_save_dir, 'dir')
        mkdir(name_data_save_dir);
    end
    
    %for tp = 1:length(time_points)
    for tp = 1:1
        time_point = time_points{end-tp+1};
        data_for_R(n).metrics(tp).time_point = time_point;
        
        tp_f = cat(2, base_dir, '/data/',  name, '/', 'MIChordAnalysis_', name, '_', time_point, '.mat');
        if ~exist(tp_f, 'file')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            load(tp_f);
            hemo_p_t1 = MI_Chord_Analysis2(end).Mipix_hemo_p;
            hemo_n_t1 = MI_Chord_Analysis2(end).Mipix_hemo_n;
            hemo_p_ff = MI_Chord_Analysis2(end).Mipix2_hemo_p;
            hemo_n_ff = MI_Chord_Analysis2(end).Mipix2_hemo_n;
            hemo_p_r2star = MI_Chord_Analysis2(end).Mipix3_hemo_p;
            hemo_n_r2star = MI_Chord_Analysis2(end).Mipix3_hemo_n;
            mi_t1 = {};
            mi_ff = {};
            mi_r2star = {};
            slice_skip = [];
            
            ele_t1 = zeros(size(hemo_p_t1, 1), 1);
            ele_ff = zeros(size(hemo_p_t1, 1), 1);
            slc_count = 1;
            for slc = 1:size(hemo_p_t1, 2)
                for i = 1:size(hemo_p_t1, 1)
                    ele_t1(i) = length([hemo_p_t1{i, slc}; hemo_n_t1{i, slc}]);
                    ele_ff(i) = length([hemo_p_ff{i, slc}; hemo_n_ff{i, slc}]);
                end
                if sum(ele_t1~=0) ~= sum(ele_ff~=0)
                    disp(cat(2, 'slice skipped: ', num2str(slc)));
                    slice_skip = [slice_skip, slc];
                else
                    for i = 1:size(hemo_p_t1, 1)
                        mi_t1{i, slc_count} = [hemo_p_t1{i, slc}; hemo_n_t1{i, slc}];
                        mi_ff{i, slc_count} = [hemo_p_ff{i, slc}; hemo_n_ff{i, slc}];
                        mi_r2star{i, slc_count} = [hemo_p_r2star{i, slc}; hemo_n_r2star{i, slc}];
                    end
                    slc_count = slc_count + 1;
                end
            end
            
            mi_t1(:,2) = circshift(mi_t1(:,2), -1);
            % mi_t1(:,3) = circshift(mi_t1(:,3), -1);
            
            data_for_R(n).metrics(tp).mi_t1 = mi_t1;
            data_for_R(n).metrics(tp).mi_ff = mi_ff;
            data_for_R(n).metrics(tp).mi_r2star = mi_r2star;
            
            mean_mi_t1 = cell(size(mi_t1));
            sd_mi_t1 = cell(size(mi_t1));
            mean_mi_ff = cell(size(mi_t1));
            mean_mi_r2star = cell(size(mi_t1));
            
            for i = 1:size(mi_t1, 1)
                for slc = 1:size(mi_t1, 2)
                    mean_mi_t1{i,slc} = mean(mi_t1{i,slc});
                    sd_mi_t1{i,slc} = std(mi_t1{i,slc});
                    
                    mi_ff_temp = mi_ff{i,slc};
                    mi_ff_temp(mi_ff_temp < 0) = 0;
                    mi_ff_temp(mi_ff_temp > 100) = 100;
                    
                    mi_r2star_temp = mi_r2star{i,slc};
                    mi_r2star_temp(mi_r2star_temp < 0) = 0;
                    
                    mean_mi_ff{i,slc} = mean(mi_ff_temp);
                    mean_mi_r2star{i,slc} = mean(mi_r2star_temp);
                end
            end
            
            data_for_R(n).metrics(tp).mean_mi_t1 = mean_mi_t1;
            data_for_R(n).metrics(tp).sd_mi_t1 = sd_mi_t1;
            data_for_R(n).metrics(tp).mean_mi_ff = mean_mi_ff;
            data_for_R(n).metrics(tp).mean_mi_r2star = mean_mi_r2star;
        end
        
        % Remote
        tp_f_rem = cat(2, base_dir, '/data/',  name, '/', 'RemoteChordAnalysis_', name, '_', time_point, '.mat');
        if ~exist(tp_f_rem, 'file')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            load(tp_f_rem);
            remote_t1 = Remote_Chord_Analysis2(end).Remotepix;
            remote_ff = Remote_Chord_Analysis2(end).Remotepix2;
            remote_r2star = Remote_Chord_Analysis2(end).Remotepix3;
                  
            slc_count = 1;
            remote_t1(:,slice_skip) = [];
            remote_ff(:,slice_skip) = [];
            remote_r2star(:,slice_skip) = [];
            
            data_for_R(n).metrics(tp).remote_t1 = remote_t1;
            data_for_R(n).metrics(tp).remote_ff = remote_ff;
            data_for_R(n).metrics(tp).remote_r2star = remote_r2star;
            
            mean_remote_t1 = cell(size(remote_t1));
            sd_remote_t1 = cell(size(remote_t1));
            mean_remote_ff = cell(size(remote_t1));
            mean_remote_r2star = cell(size(remote_t1));
            
            for i = 1:size(remote_t1, 1)
                for slc = 1:size(remote_t1, 2)
                    mean_remote_t1{i,slc} = mean(remote_t1{i,slc});
                    sd_remote_t1{i,slc} = std(remote_t1{i,slc});
                    
                    remote_ff_temp = remote_ff{i,slc};
                    remote_ff_temp(remote_ff_temp < 0) = 0;
                    remote_ff_temp(remote_ff_temp > 100) = 100;
                    
                    remote_r2star_temp = remote_r2star{i,slc};
                    remote_r2star_temp(remote_r2star_temp < 0) = 0;
                    
                    mean_remote_ff{i,slc} = mean(remote_ff_temp);
                    mean_remote_r2star{i,slc} = mean(remote_r2star_temp);
                end
            end
            
            data_for_R(n).metrics(tp).mean_remote_t1 = mean_remote_t1;
            data_for_R(n).metrics(tp).sd_remote_t1 = sd_remote_t1;
            data_for_R(n).metrics(tp).mean_remote_ff = mean_remote_ff;
            data_for_R(n).metrics(tp).mean_remote_r2star = mean_remote_r2star;
        end
    end
end


metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
save(cat(2, metrics_save_dir, 'data_for_R.mat'), 'data_for_R');

%% 
chords =size(data_for_R(1).metrics(1).mi_t1, 1);
slices = size(data_for_R(1).metrics(1).mi_t1, 2);
ele_t1 = zeros(chords, 1);
ele_ff = zeros(chords, 1);
figure();
for j = 1:slices
    for i = 1:chords
        ele_t1(i) = length(data_for_R(1).metrics(1).mi_t1{i,j});
        ele_ff(i) = length(data_for_R(1).metrics(1).mi_ff{i,j});
    end
    plot(ele_t1); hold on;
    plot(ele_ff);
    legend({'T1', 'FF'});
    pause(1);
    hold off;
end


