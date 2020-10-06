% Plotting SNR vs Fitting residual

addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FitResults_<Avg0016, Avg0001, Invivo>.mat from T2star_fitting_main.m
% SNR_<Avg0016, Avg0001, Invivo>.mat from T2star_SNR_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
label = labels{5};

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name_cell = {'18P95'};
avg_num_cell = {'Avg0016', 'Avg0001', 'Invivo'};

snr_all_remote = cell(length(avg_num_cell), length(subject_name_cell));
snr_all_air = cell(length(avg_num_cell), length(subject_name_cell));

res_all = cell(length(avg_num_cell), length(subject_name_cell));

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    for j = 1:length(avg_num_cell)
        avg_name = avg_num_cell{j};
        % subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
        
        snr = load(cat(2, subject_data_dir, 'SNR_', avg_name, '.mat'));
        snr_all_remote{j,i} = snr.SNR.snr_remote;
        snr_all_air{j,i} = snr.SNR.snr_air;
        
        FitResults = load(cat(2, subject_data_dir, 'FitResults_', avg_name, '.mat'));
        
        if j == 1
            
            res_all{j,i} = zeros(28, 1);
            for k = 1:length(res_all{j,i})
                res = FitResults.FitResults_struct(k).FitResults.res;
                res(isnan(res)) = 0;
                res_all{j,i}(k) = mean(nonzeros(res));
            end
        elseif j == 2
            res_all{j,i} = zeros(28, 1);
            for k = 1:length(res_all{j,i})
                res = FitResults.FitResults_struct(k).FitResults.res;
                res(isnan(res)) = 0;
                res_all{j,i}(k) = mean(nonzeros(res));
            end
        elseif j == 3
            res_all{j,i} = zeros(20, 1);
            for k = 1:length(res_all{j,i})
                res = FitResults.FitResults_struct(k).FitResults.res;
                res(isnan(res)) = 0;
                res_all{j,i}(k) = mean(nonzeros(res));
            end
        end
        
    end
end


%% For one case
snr_avg16_remote = mean(snr_all_remote{1,1}, 2);
snr_avg01_remote = mean(snr_all_remote{2,1}, 2);
snr_invivo_remote = mean(snr_all_remote{3,1}, 2);
snr_avg16_air = mean(snr_all_air{1,1}, 2);
snr_avg01_air = mean(snr_all_air{2,1}, 2);
snr_invivo_air = mean(snr_all_air{3,1}, 2);

%% Plot
figure();
scatter(snr_avg16_remote, res_all{1}, 72, 'filled');
hold on;
scatter(snr_avg01_remote, res_all{2}, 72, 'filled');

scatter(snr_invivo_remote, res_all{3}, 72, 'filled');
grid on;

legend({'Avg0016', 'Avg0001', 'Invivo'});
set(gca, 'FontSize', 16);

xlabel('SNR (unitless)'); ylabel('Fitting Residual');