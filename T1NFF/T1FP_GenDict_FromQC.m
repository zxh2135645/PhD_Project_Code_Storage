clear all;
close all;
% Generate a dictionary for excluding low quality images

addpath('../function/');
base_dir = uigetdir;


Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn'};

time_points = {'8WK', '12WK', '14WK' '6MO', '9MO', '1YR', '15YR'};

sequence_label = {'T1', 'T2star', 'LGE'};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

%% Iterate names
pre_QualControl = struct;

for n = 1:length(Names)
    name = Names{n};
    pre_QualControl(n).Name = name;
    pre_QualControl(n).status = struct;
    for tp = 1:length(time_points)
        time_point = time_points{tp};
        name_glob = glob(cat(2, base_dir, '/img/', name, '/overview/', name, '_', time_point, '_Slice*'));

        if ~isempty(name_glob)
            num_array = zeros(1, length(name_glob));
            for i = 1:length(name_glob)
                strs = name_glob{i};
                strings = strsplit(strs, '/');
                strings2 = strings{end};
                strings3 = split(strings2, '_');
                A = strings3{end};
                B = regexp(A,'\d*','Match');
                for ii= 1:length(B)
                    if ~isempty(B{ii})
                        Num(ii,1)=str2double(B{ii}(end));
                    else
                        Num(ii,1)=NaN;
                    end
                end
                num_array(i) = Num;
            end
            
            pre_QualControl(n).status(tp).time_point = time_point;
            for i = 1:length(name_glob)
               slc_loc = cat(2, 'Slice', num2str(num_array(i)));
               pre_QualControl(n).status(tp).(slc_loc) = 1;
            end
        else
            pre_QualControl(n).status(tp).time_point = time_point;
        end
    end
end

%% To exclude unqualified images
% Temporaly remove
% ZZ_1YR: slice1 (T1 map issue)
% pre_QualControl(5).status(6).Slice1 = 0;
% Tina good
% Sahara good
% Ryn 8WK slice1 (T1 map issue)
% Ryn 8WK 4 slices, 9MO,1YR 5 slices
% pre_QualControl(2).status(1).Slice1 = 0;
% Mojave 8WK Slice3-T2*map, Slice2-T1map
% Merry good
% Hope good
% Gobi good


% Permernantly remove
% Sunny_8WK bad T2* quality
pre_QualControl(7).status(1).Slice1 = 0;
% Queenie 6MO
pre_QualControl(8).status(4).Slice1 = 0;
% Evelyn 6MO slice1 and slice2
pre_QualControl(12).status(4).Slice1 = 0;
pre_QualControl(12).status(4).Slice2 = 0;
% Felicity 6MO slice3
pre_QualControl(11).status(4).Slice3 = 0;


% Reorder Gobi
% pre_QualControl(10).status(1).Slice1 = 0;
% TODO
% special case for Gobi and Ryn
pre_QualControl(2).status(1).Slice5 = [];
pre_QualControl(10).status(1).Slice1 = [];

save(cat(2, metrics_save_dir, 'pre_QualControl.mat'), 'pre_QualControl');