clear all;
close all;
addpath('../function/');
base_dir = uigetdir;
data_glob = glob(cat(2, base_dir, '/data/Demographic_Metrics.mat'));
load(data_glob{1});
%% Demographic_Metrics.mat (unrimmed)
hetero = struct;
count = 1;
for i = 1:length(metrics)
    for tp = 1:length(metrics(i).TimePoints)
        if ~isempty(fieldnames(metrics(i).TimePoints(tp)))
            if ~isempty(metrics(i).TimePoints(tp).time_point)
                hetero(count).name = cat(2, metrics(i).name, '_', metrics(i).TimePoints(tp).time_point);
                SliceAnalysis = metrics(i).TimePoints(tp).SliceAnalysis;
                hetero(count).Slice = struct;
                for slc = 1:length(SliceAnalysis)
                    cov_roi = SliceAnalysis(slc).std_t1_roi ./ SliceAnalysis(slc).mean_t1_roi;
                    cov_remote = SliceAnalysis(slc).std_t1_remote ./ SliceAnalysis(slc).mean_t1_remote;
                    hetero(count).Slice(slc).cov_roi = cov_roi;
                    hetero(count).Slice(slc).cov_remote = cov_remote;
                    hetero(count).Slice(slc).mean_ff_roi = SliceAnalysis(slc).mean_ff_roi;
                    hetero(count).Slice(slc).mean_ff_remote = SliceAnalysis(slc).mean_ff_remote;
                    hetero(count).Slice(slc).mean_r2star_roi = SliceAnalysis(slc).mean_r2star_roi;
                    hetero(count).Slice(slc).mean_r2star_remote = SliceAnalysis(slc).mean_r2star_remote;
                end

                count = count + 1;
            end
        end
    end
end

%% for analysis_rim_slices.mat (rimmed)
time_points = {'6MO', '9MO', '1YR', '15YR'};
data_glob = glob(cat(2, base_dir, '/data/Demographic_Metrics_rim_unaware_normalizedinROIMean_N_Entropy_bin40.mat'));
load(data_glob{1});

hetero_rim = struct;
count = 1;
for i = 1:length(metrics)
    for tp = 1:length(metrics(i).TimePoints)
        if ~isempty(fieldnames(metrics(i).TimePoints(tp)))
            if ~isempty(metrics(i).TimePoints(tp).time_point)
                hetero_rim (count).name = cat(2, metrics(i).name, '_', metrics(i).TimePoints(tp).time_point);
                SliceAnalysis = metrics(i).TimePoints(tp).SliceAnalysis;
                HeteroAnalysis = metrics(i).TimePoints(tp).HeteroAnalysis;
                hetero_rim (count).Slice = struct;
                for slc = 1:length(SliceAnalysis)
                    cov_roi = SliceAnalysis(slc).std_t1_roi ./ SliceAnalysis(slc).mean_t1_roi;
                    cov_remote = SliceAnalysis(slc).std_t1_remote ./ SliceAnalysis(slc).mean_t1_remote;
                    hetero_rim (count).Slice(slc).cov_roi = cov_roi;
                    hetero_rim (count).Slice(slc).cov_remote = cov_remote;
                    hetero_rim (count).Slice(slc).mean_ff_roi = SliceAnalysis(slc).mean_ff_roi;
                    hetero_rim (count).Slice(slc).mean_ff_remote = SliceAnalysis(slc).mean_ff_remote;
                    hetero_rim (count).Slice(slc).mean_r2star_roi = SliceAnalysis(slc).mean_r2star_roi;
                    hetero_rim (count).Slice(slc).mean_r2star_remote = SliceAnalysis(slc).mean_r2star_remote;

                    hetero_rim (count).Slice(slc).hetero_roi = HeteroAnalysis(slc).hetero_roi;
                    hetero_rim (count).Slice(slc).hetero_remote = HeteroAnalysis(slc).hetero_remote;

                    hetero_rim (count).Slice(slc).std_ff_roi = SliceAnalysis(slc).std_ff_roi;
                    hetero_rim (count).Slice(slc).std_r2star_roi = SliceAnalysis(slc).std_r2star_roi;

                    hetero_rim (count).Slice(slc).entropy_roi = HeteroAnalysis(slc).entropy_roi;
                    hetero_rim (count).Slice(slc).entropy_remote = HeteroAnalysis(slc).entropy_remote;

                end

                count = count + 1;
            end
        end
    end
end


%% for analysis_rim_slices.mat (rimmed, LGE)
time_points = {'6MO', '9MO', '1YR', '15YR'};
data_glob = glob(cat(2, base_dir, '/data/Demographic_Metrics_rim_unaware_LGE.mat'));
load(data_glob{1});

hetero_rim = struct;
count = 1;
for i = 1:length(metrics)
    for tp = 1:length(metrics(i).TimePoints)
        if ~isempty(fieldnames(metrics(i).TimePoints(tp)))
            if ~isempty(metrics(i).TimePoints(tp).time_point)
                hetero_rim (count).name = cat(2, metrics(i).name, '_', metrics(i).TimePoints(tp).time_point);
                SliceAnalysis = metrics(i).TimePoints(tp).SliceAnalysis;
                HeteroAnalysis = metrics(i).TimePoints(tp).HeteroAnalysis;
                hetero_rim (count).Slice = struct;
                for slc = 1:length(SliceAnalysis)
                    cov_roi = SliceAnalysis(slc).std_lge_roi ./ SliceAnalysis(slc).mean_lge_roi;
                    cov_remote = SliceAnalysis(slc).std_lge_remote ./ SliceAnalysis(slc).mean_lge_remote;
                    hetero_rim (count).Slice(slc).cov_roi = cov_roi;
                    hetero_rim (count).Slice(slc).cov_remote = cov_remote;
                    hetero_rim (count).Slice(slc).mean_lge_roi = SliceAnalysis(slc).mean_lge_roi;
                    hetero_rim (count).Slice(slc).mean_lge_remote = SliceAnalysis(slc).mean_lge_remote;
                    
                    hetero_rim (count).Slice(slc).hetero_roi = HeteroAnalysis(slc).hetero_roi;
                    hetero_rim (count).Slice(slc).hetero_remote = HeteroAnalysis(slc).hetero_remote;

                    hetero_rim (count).Slice(slc).std_lge_roi = SliceAnalysis(slc).std_lge_roi;
                end

                count = count + 1;
            end
        end
    end
end


%% for analysis_rim_slices.mat (rimmed, LGE, normalized)
time_points = {'6MO', '9MO', '1YR', '15YR'};
data_glob = glob(cat(2, base_dir, '/data/Demographic_Metrics_rim_unaware_LGE_normalizedinROIMean_N_Entropy.mat'));
load(data_glob{1});

hetero_rim = struct;
count = 1;
for i = 1:length(metrics)
    for tp = 1:length(metrics(i).TimePoints)
        if ~isempty(fieldnames(metrics(i).TimePoints(tp)))
            if ~isempty(metrics(i).TimePoints(tp).time_point)
                hetero_rim (count).name = cat(2, metrics(i).name, '_', metrics(i).TimePoints(tp).time_point);
                SliceAnalysis = metrics(i).TimePoints(tp).SliceAnalysis;
                HeteroAnalysis = metrics(i).TimePoints(tp).HeteroAnalysis;
                hetero_rim (count).Slice = struct;
                for slc = 1:length(SliceAnalysis)
                    cov_roi = SliceAnalysis(slc).std_lge_roi ./ SliceAnalysis(slc).mean_lge_roi;
                    cov_remote = SliceAnalysis(slc).std_lge_remote ./ SliceAnalysis(slc).mean_lge_remote;
                    hetero_rim (count).Slice(slc).cov_roi = cov_roi;
                    hetero_rim (count).Slice(slc).cov_remote = cov_remote;
                    hetero_rim (count).Slice(slc).mean_lge_roi = SliceAnalysis(slc).mean_lge_roi;
                    hetero_rim (count).Slice(slc).mean_lge_remote = SliceAnalysis(slc).mean_lge_remote;
                    
                    hetero_rim (count).Slice(slc).hetero_roi = HeteroAnalysis(slc).hetero_roi;
                    hetero_rim (count).Slice(slc).hetero_remote = HeteroAnalysis(slc).hetero_remote;

                    hetero_rim (count).Slice(slc).std_lge_roi = SliceAnalysis(slc).std_lge_roi;

                    hetero_rim (count).Slice(slc).entropy_roi = HeteroAnalysis(slc).entropy_roi;
                    hetero_rim (count).Slice(slc).entropy_remote = HeteroAnalysis(slc).entropy_remote;
                end

                count = count + 1;
            end
        end
    end
end