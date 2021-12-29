clear all;
close all;

addpath('../function/');
base_dir = uigetdir;
%% No EndoEpi
%Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16'};

time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
time_points = {'6MO', '9MO', '1YR', '15YR'};

time_points = {'7D'};
%time_points = {'8WK', '12WK', '14WK'};
%time_points = {'4MO', '6MO'};
%time_points = {'9MO', '1YR', '15YR'};

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

compound = [];
compound_remote = [];
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
        chord_values_fname = cat(2, name_data_save_dir, '/Chord_values_', name, '_', time_point, '.mat');
        
        if exist(chord_values_fname, 'file')
            load(chord_values_fname);
            
            %mean_ff_1d = nonzeros(mean_ff_array);
            mean_ff_1d = nonzeros(mean_ff_array);
            mean_r2star_1d = nonzeros(mean_r2star_array);
            mean_t1_1d = nonzeros(mean_t1_array);
            
            mean_ff_hemo_1d = nonzeros(mean_ff_hemo_array);
            mean_r2star_hemo_1d = nonzeros(mean_r2star_hemo_array);
            mean_t1_hemo_1d = nonzeros(mean_t1_hemo_array);
            
            mean_t1_remote_1d = nonzeros(mean_t1_array_remote);
            mean_ff_remote_1d = nonzeros(mean_ff_array_remote);
            mean_r2star_remote_1d = nonzeros(mean_r2star_array_remote);
            
            mean_t1_1d_n_hemo = [mean_t1_1d;mean_t1_hemo_1d];
            mean_ff_1d_n_hemo = [mean_ff_1d;mean_ff_hemo_1d];
            mean_r2star_1d_n_hemo = [mean_r2star_1d;mean_r2star_hemo_1d];
            
            %             temp = [mean_t1_1d, mean_ff_1d, mean_r2star_1d];
            %             compound = [compound; temp];
            mean_ff_1d_n_hemo(mean_r2star_1d_n_hemo > 200) = [];
            mean_ff_1d_n_hemo(mean_r2star_1d_n_hemo < 0) = [];
            mean_t1_1d_n_hemo(mean_r2star_1d_n_hemo > 200) = [];
            mean_t1_1d_n_hemo(mean_r2star_1d_n_hemo < 0) = [];
            
            mean_t1_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_t1_remote_1d(mean_r2star_remote_1d<0) = [];
            mean_ff_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_ff_remote_1d(mean_r2star_remote_1d<0) = [];
            
            mean_r2star_1d_n_hemo(mean_r2star_1d_n_hemo > 200) = [];
            mean_r2star_1d_n_hemo(mean_r2star_1d_n_hemo < 0) = [];
            mean_r2star_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_r2star_remote_1d(mean_r2star_remote_1d<0) = [];
            
            temp = [mean_t1_1d_n_hemo,mean_ff_1d_n_hemo,mean_r2star_1d_n_hemo];
            temp_remote = [mean_t1_remote_1d, mean_ff_remote_1d, mean_r2star_remote_1d];
            compound = [compound; temp];
            compound_remote = [compound_remote; temp_remote];
        end
    end
    
end
% 2D plot with labels of 
% %%
% figure();
% plot3(compound(:,1),compound(:,2),compound(:,3),'o');
%%
compound_norm = zeros(size(compound));
compound_norm(:,1) = compound(:,1)/max(compound(:,1));
compound_norm(:,2) = compound(:,2)/max(compound(:,2));
compound_norm(:,3) = compound(:,3)/max(compound(:,3));
s_array = zeros(9, 1);
for k = 2:10
    clust = kmeans(compound_norm, k);
    %clust(isnan(clust)) = [];
    s = silhouette(compound_norm,clust);
    s(isnan(s)) = [];
    s_array(k-1) = mean(s);
end

figure();
bar([2:10], s_array);
xlabel('Cluster #');
ylabel('Sihouette Coefficient');
%% 



[idx,C,sumd,D] = kmeans(compound_norm, 3);

figure();
plot3(compound(idx==1,1),compound(idx==1,2),compound(idx==1,3),'ro');
hold on;
plot3(compound(idx==2,1),compound(idx==2,2),compound(idx==2,3),'bo');
plot3(compound(idx==3,1),compound(idx==3,2),compound(idx==3,3),'ko');

%%
[idx,C,sumd,D] = kmeans(compound_norm, 9);

figure();
plot3(compound(idx==1,1),compound(idx==1,2),compound(idx==1,3),'ro');
hold on;
plot3(compound(idx==2,1),compound(idx==2,2),compound(idx==2,3),'bo');
plot3(compound(idx==3,1),compound(idx==3,2),compound(idx==3,3),'ko');
plot3(compound(idx==4,1),compound(idx==4,2),compound(idx==4,3),'go');
plot3(compound(idx==5,1),compound(idx==5,2),compound(idx==5,3),'co');
plot3(compound(idx==6,1),compound(idx==6,2),compound(idx==6,3),'mo');
plot3(compound(idx==7,1),compound(idx==7,2),compound(idx==7,3),'o','Color', [0 0.4470 0.7410]);
plot3(compound(idx==8,1),compound(idx==8,2),compound(idx==8,3),'o','Color', [0.8500 0.3250 0.0980]);
plot3(compound(idx==9,1),compound(idx==9,2),compound(idx==9,3),'o','Color', [0.9290 0.6940 0.1250]);

%% Parse endo-epi chords 
clear all; close all;

addpath('../function/');
base_dir = uigetdir;
%Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16'};

time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
%time_points = {'6MO', '9MO', '1YR', '15YR'};

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

%%
%time_points = {'6MO', '9MO'};
time_points = {'7D'};
time_points = {'8WK', '12WK', '14WK'};
time_points = {'4MO', '6MO'};
time_points = {'9MO', '1YR', '15YR'};
vec = @(x) x(:);
compound = [];
compound_remote = [];
compound_sd = [];
compound_sd_remote = [];
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
        chord_values_fname = cat(2, name_data_save_dir, '/Chord_values2_', name, '_', time_point, '.mat');
        
        if exist(chord_values_fname, 'file')
            load(chord_values_fname);
            
            mean_t1_1d_endo = nonzeros(mean_t1_array_endo);
            mean_ff_1d_endo = nonzeros(mean_ff_array_endo);
            mean_r2star_1d_endo = nonzeros(mean_r2star_array_endo);
            mean_t1_1d_epi = nonzeros(mean_t1_array_epi);
            mean_ff_1d_epi = nonzeros(mean_ff_array_epi);
            mean_r2star_1d_epi = nonzeros(mean_r2star_array_epi);
            
            mean_t1_hemo_1d_endo = nonzeros(mean_t1_hemo_array_endo);
            mean_ff_hemo_1d_endo = nonzeros(mean_ff_hemo_array_endo);
            mean_r2star_hemo_1d_endo = nonzeros(mean_r2star_hemo_array_endo);
            mean_t1_hemo_1d_epi = nonzeros(mean_t1_hemo_array_epi);
            mean_ff_hemo_1d_epi = nonzeros(mean_ff_hemo_array_epi);
            mean_r2star_hemo_1d_epi = nonzeros(mean_r2star_hemo_array_epi);
            
            mean_t1_remote_1d_endo = nonzeros(mean_t1_array_remote_endo);
            mean_ff_remote_1d_endo = nonzeros(mean_ff_array_remote_endo);
            mean_r2star_remote_1d_endo = nonzeros(mean_r2star_array_remote_endo);
            mean_t1_remote_1d_epi = nonzeros(mean_t1_array_remote_epi);
            mean_ff_remote_1d_epi = nonzeros(mean_ff_array_remote_epi);
            mean_r2star_remote_1d_epi = nonzeros(mean_r2star_array_remote_epi);
            
            mean_t1_1d = [mean_t1_1d_endo;mean_t1_1d_epi;mean_t1_hemo_1d_endo;mean_t1_hemo_1d_epi];
            mean_ff_1d = [mean_ff_1d_endo;mean_ff_1d_epi;mean_ff_hemo_1d_endo;mean_ff_hemo_1d_epi];
            mean_r2star_1d = [mean_r2star_1d_endo;mean_r2star_1d_epi;mean_r2star_hemo_1d_endo;mean_r2star_hemo_1d_epi];
            
            mean_t1_remote_1d = [mean_t1_remote_1d_endo;mean_t1_remote_1d_epi];
            mean_ff_remote_1d = [mean_ff_remote_1d_endo;mean_ff_remote_1d_epi];
            mean_r2star_remote_1d = [mean_r2star_remote_1d_endo;mean_r2star_remote_1d_epi];
            
            mean_ff_1d(mean_r2star_1d > 200) = [];
            mean_ff_1d(mean_r2star_1d < 0) = [];
            mean_t1_1d(mean_r2star_1d > 200) = [];
            mean_t1_1d(mean_r2star_1d < 0) = [];
            
            mean_t1_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_t1_remote_1d(mean_r2star_remote_1d<0) = [];
            mean_ff_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_ff_remote_1d(mean_r2star_remote_1d<0) = [];
            
            
            
            % SD
            sd_t1_1d_endo = vec(sd_t1_array_endo(sd_t1_array_endo ~= 1000));
            sd_ff_1d_endo = vec(sd_ff_array_endo(sd_ff_array_endo ~= 1000));
            sd_r2star_1d_endo = vec(sd_r2star_array_endo(sd_r2star_array_endo ~= 1000));
            sd_t1_1d_epi = vec(sd_t1_array_epi(sd_t1_array_epi ~= 1000));
            sd_ff_1d_epi = vec(sd_ff_array_epi(sd_ff_array_epi ~= 1000));
            sd_r2star_1d_epi = vec(sd_r2star_array_epi(sd_r2star_array_epi ~= 1000));
            
            sd_t1_hemo_1d_endo = vec(sd_t1_hemo_array_endo(sd_t1_hemo_array_endo ~= 1000));
            sd_ff_hemo_1d_endo = vec(sd_ff_hemo_array_endo(sd_ff_hemo_array_endo ~= 1000));
            sd_r2star_hemo_1d_endo = vec(sd_r2star_hemo_array_endo(sd_r2star_hemo_array_endo ~= 1000));
            sd_t1_hemo_1d_epi = vec(sd_t1_hemo_array_epi(sd_t1_hemo_array_epi ~= 1000));
            sd_ff_hemo_1d_epi = vec(sd_ff_hemo_array_epi(sd_ff_hemo_array_epi ~= 1000));
            sd_r2star_hemo_1d_epi = vec(sd_r2star_hemo_array_epi(sd_t1_hemo_array_epi ~= 1000));
            
            sd_t1_remote_1d_endo = vec(sd_t1_array_remote_endo(sd_t1_array_remote_endo ~= 1000));
            sd_ff_remote_1d_endo = vec(sd_ff_array_remote_endo(sd_ff_array_remote_endo ~= 1000));
            sd_r2star_remote_1d_endo = vec(sd_r2star_array_remote_endo(sd_r2star_array_remote_endo ~= 1000));
            sd_t1_remote_1d_epi = vec(sd_t1_array_remote_epi(sd_t1_array_remote_epi ~= 1000));
            sd_ff_remote_1d_epi = vec(sd_ff_array_remote_epi(sd_ff_array_remote_epi ~= 1000));
            sd_r2star_remote_1d_epi = vec(sd_r2star_array_remote_epi(sd_r2star_array_remote_epi ~= 1000));
            
            sd_t1_1d = [sd_t1_1d_endo;sd_t1_1d_epi;sd_t1_hemo_1d_endo;sd_t1_hemo_1d_epi];
            sd_ff_1d = [sd_ff_1d_endo;sd_ff_1d_epi;sd_ff_hemo_1d_endo;sd_ff_hemo_1d_epi];
            sd_r2star_1d = [sd_r2star_1d_endo;sd_r2star_1d_epi;sd_r2star_hemo_1d_endo;sd_r2star_hemo_1d_epi];
            
            sd_t1_remote_1d = [sd_t1_remote_1d_endo;sd_t1_remote_1d_epi];
            sd_ff_remote_1d = [sd_ff_remote_1d_endo;sd_ff_remote_1d_epi];
            sd_r2star_remote_1d = [sd_r2star_remote_1d_endo;sd_r2star_remote_1d_epi];
            
            sd_ff_1d(mean_r2star_1d > 200) = [];
            sd_ff_1d(mean_r2star_1d < 0) = [];
            sd_t1_1d(mean_r2star_1d > 200) = [];
            sd_t1_1d(mean_r2star_1d < 0) = [];
            sd_r2star_1d(mean_r2star_1d > 200) = [];
            sd_r2star_1d(mean_r2star_1d < 0) = [];
            
            sd_t1_remote_1d(mean_r2star_remote_1d>200) = [];
            sd_t1_remote_1d(mean_r2star_remote_1d<0) = [];
            sd_ff_remote_1d(mean_r2star_remote_1d>200) = [];
            sd_ff_remote_1d(mean_r2star_remote_1d<0) = [];
            sd_r2star_remote_1d(mean_r2star_remote_1d>200) = [];
            sd_r2star_remote_1d(mean_r2star_remote_1d<0) = [];
            
            mean_r2star_1d(mean_r2star_1d > 200) = [];
            mean_r2star_1d(mean_r2star_1d < 0) = [];
            mean_r2star_remote_1d(mean_r2star_remote_1d>200) = [];
            mean_r2star_remote_1d(mean_r2star_remote_1d<0) = [];
            
            temp = [mean_t1_1d,  mean_ff_1d, mean_r2star_1d];
            temp_remote = [mean_t1_remote_1d, mean_ff_remote_1d, mean_r2star_remote_1d];
            compound = [compound; temp];
            compound_remote = [compound_remote; temp_remote];
            
            temp = [sd_t1_1d,  sd_ff_1d, sd_r2star_1d];
            temp_remote = [sd_t1_remote_1d, sd_ff_remote_1d, sd_r2star_remote_1d];
            compound_sd = [compound_sd; temp];
            compound_sd_remote = [compound_sd_remote; temp_remote];
        end
    end
    
end


compound_t1 = compound(:,1);
compound_ff = compound(:,2);
compound_r2star = compound(:,3);
compound_t1_remote = compound_remote(:,1);
compound_ff_remote = compound_remote(:,2);
compound_r2star_remote = compound_remote(:,3);
%%
figure();
plot3(compound(:,1),compound(:,2),compound(:,3),'o');
ylim([0 70]);
hold on;
plot3(compound_remote(:,1),compound_remote(:,2),compound_remote(:,3),'rx');
xlim([800 1700]);

figure();
plot3(compound_sd(:,1),compound_sd(:,2),compound_sd(:,3),'o');


compound_sd_t1 = compound_sd(:,1);
compound_sd_t1_label = compound_sd_t1<100;

figure();
plot3(compound_t1(compound_sd_t1_label),compound_ff(compound_sd_t1_label),compound_r2star(compound_sd_t1_label), 'o');
ylim([0 70]);

mdl = fitlm(compound_ff(compound_sd_t1_label), compound_r2star(compound_sd_t1_label));
%% 2D plot
figure();
subplot(2,2,1);
plot(compound(:,1),compound(:,2), 'o'); xlabel('T1 (ms)'); ylabel('FF (%)');xlim([800 1700]); ylim([0 70]);
hold on; 
plot(compound_remote(:,1),compound_remote(:,2), 'o'); xlabel('T1 (ms)'); ylabel('FF (%)');xlim([800 1700]); ylim([0 70]);

subplot(2,2,2);
plot(compound(:,1),compound(:,3), 'o'); xlabel('T1 (ms)'); ylabel('R2star (s^{-1})');xlim([800 1700]); ylim([0 200]);
hold on;
plot(compound_remote(:,1),compound_remote(:,3), 'o'); xlabel('T1 (ms)'); ylabel('R2star (s^{-1})');xlim([800 1700]); ylim([0 200]);
subplot(2,2,3);
plot(compound(:,2),compound(:,3), 'o'); xlabel('FF (%)'); ylabel('R2star (s^{-1})');xlim([0 70]); ylim([0 200]);
hold on;
plot(compound_remote(:,2),compound_remote(:,3), 'o'); xlabel('FF (%)'); ylabel('R2star (s^{-1})');xlim([0 70]); ylim([0 200]);
subplot(2,2,4);
histogram(compound(:,1)); hold on;
histogram(compound_remote(:,1));
xlabel('T1 (ms)');
legend({'MI', 'Remote'})
%%
compound_norm = zeros(size(compound));
compound_norm(:,1) = compound(:,1)/max(compound(:,1));
compound_norm(:,2) = compound(:,2)/max(compound(:,2));
compound_norm(:,3) = compound(:,3)/max(compound(:,3));
s_array = zeros(9, 1);
for k = 2:10
    clust = kmeans(compound_norm, k);
    %clust(isnan(clust)) = [];
    s = silhouette(compound_norm,clust);
    s(isnan(s)) = [];
    s_array(k-1) = mean(s);
end

figure();
bar([2:10], s_array);
xlabel('Cluster #');
ylabel('Sihouette Coefficient');


[idx,C,sumd,D] = kmeans(compound_norm, 3);

figure();
plot3(compound(idx==1,1),compound(idx==1,2),compound(idx==1,3),'ro');
hold on;
plot3(compound(idx==2,1),compound(idx==2,2),compound(idx==2,3),'bo');
plot3(compound(idx==3,1),compound(idx==3,2),compound(idx==3,3),'go');

[idx,C,sumd,D] = kmeans(compound_norm, 6);

figure();
plot3(compound(idx==1,1),compound(idx==1,2),compound(idx==1,3),'ro');
hold on;
plot3(compound(idx==2,1),compound(idx==2,2),compound(idx==2,3),'bo');
plot3(compound(idx==3,1),compound(idx==3,2),compound(idx==3,3),'ko');
plot3(compound(idx==4,1),compound(idx==4,2),compound(idx==4,3),'go');
plot3(compound(idx==5,1),compound(idx==5,2),compound(idx==5,3),'co');
plot3(compound(idx==6,1),compound(idx==6,2),compound(idx==6,3),'mo');
%%
r2star_label = zeros(size(compound, 1),1);
r2star_label = r2star_label + double(compound(:,3) <= 20);
r2star_label = r2star_label + double(compound(:,3) <= 40);
r2star_label = r2star_label + double(compound(:,3) <= 60);
r2star_label = r2star_label + double(compound(:,3) <= 80);
r2star_label = r2star_label + double(compound(:,3) <= 100);
r2star_label = r2star_label + double(compound(:,3) <= 200);


figure();
t1 = compound(:,1);
ff = compound(:,2);
plot(t1(r2star_label == 6),ff(r2star_label == 6),'o');
hold on;
plot(t1(r2star_label == 5),ff(r2star_label == 5),'+');
plot(t1(r2star_label == 4),ff(r2star_label == 4),'*');
plot(t1(r2star_label == 3),ff(r2star_label == 3),'^');
plot(t1(r2star_label == 2),ff(r2star_label == 2),'s');
plot(t1(r2star_label == 1),ff(r2star_label == 1),'p');
ylim([0 70]);

%% Parse pixel-wise 
clear all; close all;

addpath('../function/');
base_dir = uigetdir;
%Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16'};

time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
%time_points = {'6MO', '9MO', '1YR', '15YR'};

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

%%
%time_points = {'6MO', '9MO'};
time_points = {'7D'};
%time_points = {'8WK', '12WK', '14WK'};
%time_points = {'4MO', '6MO'};
%time_points = {'9MO', '1YR', '15YR'};
vec = @(x) x(:);
compound = [];
compound_remote = [];
compound_sd = [];
compound_sd_remote = [];
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
        chord_values_fname = cat(2, name_data_save_dir, '/Chord_values_pixelwise_', name, '_', time_point, '.mat');
        
        if exist(chord_values_fname, 'file')
            load(chord_values_fname);
            
            mean_ff_1d = nonzeros(mean_ff_array);
            mean_r2star_1d = nonzeros(mean_r2star_array);
            mean_t1_1d = nonzeros(mean_t1_array);
            
            mean_ff_hemo_1d = nonzeros(mean_ff_hemo_array);
            mean_r2star_hemo_1d = nonzeros(mean_r2star_hemo_array);
            mean_t1_hemo_1d = nonzeros(mean_t1_hemo_array);
            
            mean_t1_remote_1d = nonzeros(mean_t1_array_remote);
            mean_ff_remote_1d = nonzeros(mean_ff_array_remote);
            mean_r2star_remote_1d = nonzeros(mean_r2star_array_remote);
            
            mean_t1_1d_n_hemo = [mean_t1_1d;mean_t1_hemo_1d];
            mean_ff_1d_n_hemo = [mean_ff_1d;mean_ff_hemo_1d];
            mean_r2star_1d_n_hemo = [mean_r2star_1d;mean_r2star_hemo_1d];
            
            %             temp = [mean_t1_1d, mean_ff_1d, mean_r2star_1d];
            %             compound = [compound; temp];
            mean_ff_1d_n_hemo(mean_r2star_1d_n_hemo > 200 | mean_r2star_1d_n_hemo < 0) = [];
            mean_t1_1d_n_hemo(mean_r2star_1d_n_hemo > 200 | mean_r2star_1d_n_hemo < 0) = [];
            
            mean_t1_remote_1d(mean_r2star_remote_1d>200 | mean_r2star_remote_1d<0) = [];
            mean_ff_remote_1d(mean_r2star_remote_1d>200 | mean_r2star_remote_1d<0) = [];
            
            mean_r2star_1d_n_hemo(mean_r2star_1d_n_hemo > 200 | mean_r2star_1d_n_hemo < 0) = [];
            mean_r2star_remote_1d(mean_r2star_remote_1d>200 | mean_r2star_remote_1d<0) = [];
            
            temp = [mean_t1_1d_n_hemo,mean_ff_1d_n_hemo,mean_r2star_1d_n_hemo];
            temp_remote = [mean_t1_remote_1d, mean_ff_remote_1d, mean_r2star_remote_1d];
            compound = [compound; temp];
            compound_remote = [compound_remote; temp_remote];
        end
    end
    
end


compound_t1 = compound(:,1);
compound_ff = compound(:,2);
compound_r2star = compound(:,3);
compound_t1_remote = compound_remote(:,1);
compound_ff_remote = compound_remote(:,2);
compound_r2star_remote = compound_remote(:,3);

%% 2D plot
figure();
subplot(2,2,1);
plot(compound(:,1),compound(:,2), 'o'); xlabel('T1 (ms)'); ylabel('FF (%)');xlim([800 1700]); ylim([0 70]);
hold on; 
plot(compound_remote(:,1),compound_remote(:,2), 'o'); xlabel('T1 (ms)'); ylabel('FF (%)');xlim([800 1700]); ylim([0 70]);

subplot(2,2,2);
plot(compound(:,1),compound(:,3), 'o'); xlabel('T1 (ms)'); ylabel('R2star (s^{-1})');xlim([800 1700]); ylim([0 200]);
hold on;
plot(compound_remote(:,1),compound_remote(:,3), 'o'); xlabel('T1 (ms)'); ylabel('R2star (s^{-1})');xlim([800 1700]); ylim([0 200]);
subplot(2,2,3);
plot(compound(:,2),compound(:,3), 'o'); xlabel('FF (%)'); ylabel('R2star (s^{-1})');xlim([0 70]); ylim([0 200]);
hold on;
plot(compound_remote(:,2),compound_remote(:,3), 'o'); xlabel('FF (%)'); ylabel('R2star (s^{-1})');xlim([0 70]); ylim([0 200]);
subplot(2,2,4);
histogram(compound(:,1)); hold on;
histogram(compound_remote(:,1));
xlabel('T1 (ms)');
legend({'MI', 'Remote'})