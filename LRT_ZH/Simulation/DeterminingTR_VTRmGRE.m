close all;
clear all;

addpath('../../function/');
addpath('../../lib_EPGX/')
addpath('../../EPGX-src/')
addpath('../../BlochSimDemo/');
addpath('../../M219/');


%% After Discuss with Randy  (This is the true main 07/26/2025)
% No need to iterate TI, but change the N (#ofSegment)
df = 0;
N_ti = 10;
N_array = [100:20:300];
% T1_array = 300:2:800;
T1_array = 200:5:1200;

flip_prep_array = [pi];
Nshot = 1;
TI = 10.5;
TR1 = 3.6;
TR = 13.6 + 3.6;
TR2 = TR - TR1;
alpha = 5 / 180 * pi;
dt = 0.1;
T2 = 45;
clear Mz_recov_TE_dict;
Mz_recov_TE_dict = cell(length(T1_array), length(N_array), length(flip_prep_array));
for n = 1:length(N_array)
    N = N_array(n);
for tff = 1:length(flip_prep_array)
    flip_prep = flip_prep_array(tff);
        for t11 = 1:length(T1_array)
            T1 = T1_array(t11);
            disp(['T1: ', num2str(T1), ' ms']);
            disp(['N: ', num2str(N), ' seg']);
            %tic;
            total_time = (TI*Nshot + TR * N) * N_ti;
            t = [dt:dt:total_time];
            
            Nstep_TI = TI / dt;
            Nstep_mgre1 = round(TR1 * N / dt / Nshot);
            Nstep_mgre2 = round(TR2 * N / dt / Nshot);
            E1 = exp(-dt / T1);
            E2 = exp(-dt / T2);
            Mz_recov_array = [];
            Mxy_recov_array = [];
            M0 = [0 0 cos(flip_prep)];
            Mz0 = M0(3);
            Mxy0 = M0(1);
            Nprep = 5;
            for p = 1:Nprep
                for s = 1:Nshot
                    Mz_TI_recov_array = zeros(1, Nstep_TI);
                    Mxy_TI_recov_array = zeros(1, Nstep_TI);
                    for i = 1:Nstep_TI
                        if i == 1
                            Mz_TI_recov_array(i) = Mz0;
                            Mxy_TI_recov_array(i) = Mxy0;
                        else
                            Mz_TI_recov_array(i) = 1 - E1 + Mz_TI_recov_array(i-1) * E1;
                            Mxy_TI_recov_array(i) = Mxy_TI_recov_array(i-1)*E2;
                        end
                    end
                    
                    % TR1 (Navigator line)
                    Mz1 = Mz_TI_recov_array(end);
                    Mxy1 = Mxy_TI_recov_array(end);
                    Mz_mgre_recov_array = zeros(1, Nstep_mgre1);
                    Mxy_mgre_recov_array = zeros(1, Nstep_mgre1);
                    Nstep_TR1 = TR1/dt;
                    for i = 1:Nstep_mgre1
                        
                        if mod(i-1, Nstep_TR1) == 0
                            if i == 1
                                Mz_mgre_recov_array(i) = 1 - E1 + Mz1 * cos(alpha) * E1;
                                Mxy_mgre_recov_array(i) = Mz1 * sin(alpha) * E2;
                            else
                                Mz_mgre_recov_array(i) = 1 - E1 + Mz_mgre_recov_array(i-1) * cos(alpha) * E1;
                                Mxy_mgre_recov_array(i) = Mz_mgre_recov_array(i-1) * sin(alpha) * E2;
                            end
                        else
                            Mz_mgre_recov_array(i) = 1 - E1 + Mz_mgre_recov_array(i-1) * E1;
                            Mxy_mgre_recov_array(i) = Mxy_mgre_recov_array(i-1) * E2;
                        end
                    end


                    % TR2 (mGRE lines)
                    Mz2 = Mz_mgre_recov_array(end);
                    Mxy2 = Mxy_mgre_recov_array(end);
                    Mz_mgre_recov_array2 = zeros(1, Nstep_mgre2);
                    Mxy_mgre_recov_array2 = zeros(1, Nstep_mgre2);
                    Nstep_TR2 = TR2/dt;
                    for i = 1:Nstep_mgre2

                        if mod(i-1, Nstep_TR2) == 0
                            if i == 1
                                Mz_mgre_recov_array2(i) = 1 - E1 + Mz2 * cos(alpha) * E1;
                                Mxy_mgre_recov_array2(i) = Mz2 * sin(alpha) * E2;
                            else
                                Mz_mgre_recov_array2(i) = 1 - E1 + Mz_mgre_recov_array2(i-1) * cos(alpha) * E1;
                                Mxy_mgre_recov_array2(i) = Mz_mgre_recov_array2(i-1) * sin(alpha) * E2;
                            end
                        else
                            Mz_mgre_recov_array2(i) = 1 - E1 + Mz_mgre_recov_array2(i-1) * E1;
                            Mxy_mgre_recov_array2(i) = Mxy_mgre_recov_array2(i-1) * E2;
                        end
                    end


                    Mz0 = cos(flip_prep)*Mz_mgre_recov_array2(end);
                    Mxy0 = cos(flip_prep)*Mxy_mgre_recov_array2(end);
                    Mz_recov_array = [Mz_recov_array, Mz_TI_recov_array, Mz_mgre_recov_array, Mz_mgre_recov_array2];
                    Mxy_recov_array = [Mxy_recov_array, Mxy_TI_recov_array, Mxy_mgre_recov_array, Mxy_mgre_recov_array2];
                end
            end
            % TE_array = [1.4, 3.4, 5.4, 7.4, 9.4, 11.4, 13.4, 15.4];
            TE_array = [1.6]; % Only looking at first-echo
            Nstep_TE_array = round(TE_array / dt);
            Nstep_TR = round(TR/dt);
            Nstep_TE_indices = [];
            for i = 1:Nprep
                for k = 1:Nshot
                    for j = 1:(N/Nshot)
                        Nstep_TE_indices = [Nstep_TE_indices, (Nshot*(i-1)+k-1)*(N/Nshot)*Nstep_TR + (Nshot*(i-1)+k)*Nstep_TI + (j-1)*Nstep_TR + Nstep_TE_array];
                    end
                end
            end
            %toc;
            Mz_recov_TE = Mz_recov_array(Nstep_TE_indices);
            Mz_recov_TE_dict{t11, n, tff} = Mz_recov_TE;
        end
    end
end

%% By Analytical solution (Tianle)

Modd_ss =  M0 * (1 - exp(-TR2/T1) + exp(-TR2/T1) * cos(alpha) - exp(-(TR1+TR2)/T1)*cos(alpha)) / (1 - exp(-(TR1+TR2)/T1)*(cos(alpha).^2));
Meven_ss = M0 * (1 - exp(-TR1/T1) + exp(-TR1/T1) * cos(alpha) - exp(-(TR1+TR2)/T1)*cos(alpha)) / (1 - exp(-(TR1+TR2)/T1)*(cos(alpha).^2));

T1 = 600;
N_array = [20:2:50];
N = 300;
TE = 1.6;
S_2n = zeros(N,3);
for n = 1:N
    %S_2n(n,:) = (Meven_ss + Modd_ss * exp(-TR1/T1) * cos(alpha) * (-2) * (exp(-(TR1+TR2)/T1)*(cos(alpha)^2))^(n-1));
    S_2n(n,:) = (Meven_ss + Modd_ss * exp(-TR1/T1) * cos(alpha) * (-2) * (exp(-(TR1+TR2)/T1)*(cos(alpha)^2))^(n-1)) * exp(-TE/T2) * sin(alpha);
end

figure();
plot(S_2n);
%% Adding noise 100 times
% figure(); plot(Mz_recov_TE_dict{1, 1, 1});
iters = 100;
Mz_recov_TE_dict_noise = cell(length(T1_array), length(N_array), length(flip_prep_array), iters);
T1_array2 = [200, 300, 400, 500, 600, 700, 800, 900, 1000];
for its = 1:iters
    for n = 1:length(N_array)
        N = N_array(n);
        for tff = 1:length(flip_prep_array)
            flip_prep = flip_prep_array(tff);
            for t11 = 1:length(T1_array2)
                T1 = T1_array2(t11);
                idx = find(T1 == T1_array);
                % Mz_recov_TE_dict_noise{t11, n, tff, its} = Mz_recov_TE_dict{idx, n, tff} + 0.1*randn(1,length(Mz_recov_TE_dict{idx, n, tff}));
                Mz_recov_TE_dict_noise{t11, n, tff, its} = Mz_recov_TE_dict{idx, n, tff} + 0.1*randn(1,length(Mz_recov_TE_dict{idx, n, tff}));
            end
        end
    end
end

%% Plot
Nseg = 200;
n = find(N_array == Nseg);
N = N_array(n);
total_time = (TI*Nshot + TR * N) * N_ti;
t = [dt:dt:total_time];
figure();
T1 = 800;
n1 = find(T1_array == T1)
n2 = find(T1_array2 == T1);
plot(Mz_recov_TE_dict_noise{n2, n, 1, 1});
hold on;
plot(Mz_recov_TE_dict{n1, n, 1}, 'LineWidth', 1.5);
ylim([-1.2 1.2]);


figure();
for i = 1:length(T1_array2)
    subplot(3,3,i)
    T1 = T1_array2(i);
    n1 = find(T1_array == T1)
    n2 = find(T1_array2 == T1);
    plot(Mz_recov_TE_dict_noise{n2, n, 1, 1});
    hold on;
    plot(Mz_recov_TE_dict{n1, n, 1}, 'LineWidth', 1.5);
    ylim([-1.2 1.2]);
    title(cat(2, 'T1 = ', num2str(T1), ' ms'));
end
%%
T1 = 800;
n1 = find(T1_array == T1)
n2 = find(T1_array2 == T1);

N = N_array(n);
total_time = (TI*Nshot + TR * N) * N_ti;
t = [dt:dt:total_time];
figure();
plot(abs(Mz_recov_TE_dict_noise{n2, n, 1, 1}(end-N+1:end)));
hold on;
plot(abs(Mz_recov_TE_dict{n1, n, 1}(end-N+1:end)), 'LineWidth', 1.5);
ylim([0 1.2]);


N_cutoff = 50;
SNR = mean(abs(Mz_recov_TE_dict_noise{n2, n, 1, 1}(end-N_cutoff+1:end))) / std(abs(Mz_recov_TE_dict_noise{n2, n, 1, 1}(end-N_cutoff+1:end)))



%% Dictionary Matching
SumSQRT = @(x,y) sum(sqrt((x-y).^2))/length(x);

error_mat_cat = [];
best_fit_mat = zeros(length(T1_array2), length(N_array), length(flip_prep_array), iters);
for t111 = 1:length(T1_array2)
    Mz_recov_TE_noise = cell(length(N_array), length(flip_prep_array), iters);
    for n = 1:length(N_array)
        N = N_array(n);
        for tff = 1:length(flip_prep_array)
            flip_prep = flip_prep_array(tff);
            for its = 1:iters
                temp = abs(Mz_recov_TE_dict_noise{t111,n,tff,its}(end-N+1:end)); % Changed to absolute value
                % Added for a more realistic case
                %if tff == 2
                %    idx = find(temp == min(temp));
                %    temp(1:idx) = -temp(1:idx);
                %end
                Mz_recov_TE_noise{n,tff,its} = temp;
            end
        end
    end
    
    RMS_mat = zeros(length(T1_array), length(N_array), length(flip_prep_array), iters);
    for n = 1:length(N_array)
        N = N_array(n);
        for tff = 1:length(flip_prep_array)
            flip_prep = flip_prep_array(tff);
            for t11 = 1:length(T1_array)
                T1 = T1_array(t11);
                for its = 1:iters
                    RMS_mat(t11, n, tff, its) = SumSQRT(Mz_recov_TE_noise{n,tff,its}, abs(Mz_recov_TE_dict{t11,n,tff}(end-N+1:end)));
                end
            end
        end
    end
    
    for n = 1:length(N_array)
        N = N_array(n);
        for tff = 1:length(flip_prep_array)
            flip_prep = flip_prep_array(tff);
            for its = 1:iters
                best_fit_mat(t111, n,tff,its) = T1_array(find(RMS_mat(:,n,tff,its) == min(RMS_mat(:,n,tff,its))));
            end
        end
    end
%     error_mat = abs(best_fit_mat - T1_array2(t111));
%     subplot(2,2,t111);
%     imagesc(error_mat); colorbar; caxis([0 50]);
%     title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
end

%% Plot precision and accuracy
best_fit_mat_stats = zeros(length(T1_array2), length(N_array), length(flip_prep_array), 2);
error_mat = zeros(length(T1_array2), length(N_array), length(flip_prep_array), 2);

for t111 = 1:length(T1_array2)
    for n = 1:length(N_array)
        for tff = 1:length(flip_prep_array)
            best_fit_mat_stats(t111,n,tff,1) = mean(squeeze(best_fit_mat(t111,n,tff,:)));
            best_fit_mat_stats(t111,n,tff,2) = std(squeeze(best_fit_mat(t111,n,tff,:)));
        end
    end
    error_mat(t111,:,:,1) = abs(best_fit_mat_stats(t111,:,:,1) - T1_array2(t111));
    error_mat(t111,:,:,2) = best_fit_mat_stats(t111,:,:,2);
end


figure();
ytlbl = [160:20:700];
yt = ytlbl;
xt = [1, 6, 11];
xtlbl = {'100', '200', '300'};
for t111 = 1:length(T1_array2)
    subplot(3,3,t111);
    imagesc(squeeze(error_mat(t111,:,:,1))); colorbar; caxis([0 10]);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl);
end

figure();
for t111 = 1:length(T1_array2)
    subplot(3,3,t111);
    imagesc(squeeze(error_mat(t111,:,:,2))); colorbar; caxis([0 30]);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl);
end

%% Plot curves to find optimal
figure();
xtlbl = [160:20:700];
xt = 1:length(xtlbl);
for t111 = 1:length(T1_array2)
    subplot(3,3,t111);
    plot(error_mat(t111,:,1,1)); ylim([0 30]); hold on;
    y = medfilt1(error_mat(t111,:,1,1), 5);
    plot(y);
    ylabel('Absolute Error (ms)'); xlabel('# of Segments');
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
end

%%
figure();
xtlbl = [100:10:500];
xt = 1:length(xtlbl);
seg_array = [100:10:500];
idx_array = 1:length(seg_array);
for t111 = 1:length(T1_array2)
    subplot(3,3,t111);
    plot(error_mat(t111,:,1,2)); ylim([0 35]); hold on; %ylim([0 0.05]);% 
    y_temp = [repmat(error_mat(t111,1,1,2), [1, 10]), error_mat(t111,:,1,2)];
    y_med = medfilt1(y_temp, 10);
    y_med = y_med(11:end);
    plot(y_med, 'LineWidth', 1.5);
    y = mean(error_mat(t111,end-21:end,1,2));
    y_sd = std(error_mat(t111,end-21:end,1,2));
    yline(y);
    yline(y+2*y_sd, 'k--', 'LineWidth', 1);
    yline(y-2*y_sd, 'k--', 'LineWidth', 1);
    xx = idx_array(find(y_med < y+2*y_sd == 1, 1));
    xline(xx);
    xxlbl = seg_array(find(y_med < y+2*y_sd == 1, 1));
    ylabel('Absolute Error (ms)'); xlabel('# of Segments');
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
    text(xx+5, 25, ['X = ', num2str(xxlbl)]);
end