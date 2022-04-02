close all
clear all

addpath('../function/');
base_dir = uigetdir;

%Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16'};

time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
time_points = {'6MO', '9MO', '1YR', '15YR'};

time_points = {'7D'};
time_points = {'8WK', '12WK', '14WK'};
time_points = {'4MO', '6MO'};
time_points = {'9MO', '1YR', '15YR'};
time_points = {'1YR', '15YR'};

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');
%%
kmean_input = [];
count = 0;
slc_num_old = 0;
name_labels = {};
for n = 1:length(Names)
    %for n = 4:4
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    
    name_data_save_dir = cat(2, data_save_dir, name);
    
    for tp = 1:length(time_points)
        %for tp = 3:3
        time_point = time_points{end-tp+1};
        start_fname = cat(2, name_data_save_dir, '/FF_Stratify_', name, '_', time_point, '.mat');
        tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
        if ~exist(tp_dir2, 'dir')
            mkdir(tp_dir2);
        end
        if exist(start_fname, 'file')
            load(start_fname);
            display(cat(2, name, ' ', time_point));
            sizes = size(strat_cell);
            
            figure();
            nn = ceil(sqrt(sizes(2)));
            for col = 1:sizes(2)
                subplot(nn, nn, col);
                for row = 1:sizes(1)
                    strat_sort = sort(strat_cell{row, col});
                    len_strat = length(strat_sort);
                    plot(strat_sort);
                    hold on;
                    if 0.1*len_strat > 1
                        ten_perc = round(0.1*len_strat);
                        max_k = maxk(strat_sort, ten_perc);
                        min_k = mink(strat_sort, ten_perc);
                        mean_k = mean(strat_sort);
                        sd_k = std(strat_sort);
                    end
                    
                    if isnan(mean(strat_sort))
                        mean_k = 0;
                        sd_k = 0;
                        max_k = 0;
                        min_k = 0;
                    end
                    
                    kmean_input(col+slc_num_old, ((row-1)*4+1):((row-1)*4+4)) = [mean_k,sd_k,min_k(end),max_k(end)];
                    
                end
                legend({'0-5 %','5-10 %', '10-15 %', '15-20 %', '20-25 %', '25-30 %', '>30 %'}, 'Location', 'SouthEast');
                title(['Slice = ', num2str(col)]);
                
                name_labels{col+slc_num_old} = cat(2, name, '_', time_point, '_Slice', num2str(col));
                
            end
            
            saveas(gcf, cat(2, tp_dir2, 'Strat_Curve_', name, '_', time_point , '.png'));
            slc_num_old = slc_num_old + col;
            count = count + 1;
        end
        
        
    end
    close all;
    
end

%%
s_array = zeros(9, 1);
for k = 2:10
    clust = kmeans(kmean_input, k);
    %clust(isnan(clust)) = [];
    s = silhouette(kmean_input,clust);
    s(isnan(s)) = [];
    s_array(k-1) = mean(s);
end

figure();
bar([2:10], s_array);
xlabel('Cluster #');
ylabel('Sihouette Coefficient');

%% 
name_labels2 = name_labels';
[idx,C,sumd,D] = kmeans(kmean_input, 2);

% temp_save = 
% dbscan
idx_dbscan = dbscan(kmean_input,2000,5);

