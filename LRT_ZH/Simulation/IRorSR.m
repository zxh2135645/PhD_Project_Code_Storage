close all;
clear all;

addpath('../../function/');
addpath('../../lib_EPGX/')
addpath('../../EPGX-src/')
addpath('../../BlochSimDemo/');
addpath('../../M219/');

%% T1 mapping MOLLI
lSegment = 192;
TR = 17.1;
TI_array = 10.5 + TR * lSegment * ((1:10) - 1);
figure();
b1 = 750;
PhaseEnc = 192;
num_rampup = 0;

HR = 80;
window = round(60*1000 / HR);
RAMP_DOWN = 0;
acq_win = TR*(PhaseEnc+num_rampup+RAMP_DOWN);


line([0,10],[0,0],'Color', 'black')
hold on;

% 3 bSSFP readout
for i = 1:length(TI_array)
    trigger = TI_array(i);
    plot([trigger trigger]/1000, [-b1 b1], 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8])
    rectangle('Position', [trigger/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
end

xlabel('Time (s)')
ylabel('B_1 (Hz)')

%% Fig. 2
TI_array = 10.5 + TR * lSegment * ((1:10) - 1);
num_rampup = 0;
% npulse = 58 + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1000;
T2 = 45;
df = 0; %Hz


prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
M0 = [0 0 1]';
TD = trigger; % ms

RAMP_DOWN = 0;
npulse = 192 + num_rampup + RAMP_DOWN;

%%% User inputs for adiabatic pulse:
adiabatic.mu = 5;   % Phase modulation parameter [dimensionless]
adiabatic.beta1 = 750;   % Frequency modulation parameter [rad/s]
% adiabatic.pulseWidth = 10*2;   % For temporal resolution 0.1ms
adiabatic.pulseWidth = 10.24*2;  % RF pulse duration [ms] % According to siemens 3T
adiabatic.A0 = 0.12;
%% Main Body (8 echo)
figure('Position', [100 100 1600 1000]);
T1_array = [400, 500, 600, 700];

for ttt = 1:length(T1_array)
    T1 = T1_array(ttt);
    TI_array = 10.5 + TR * lSegment * ((1:5) - 1);
    flip1 = pi/2;
    flip2 = pi;
    M01 = [0 0 cos(flip1)];
    Mz01 = M01(3);
    M02 = [0 0 cos(flip2)];
    Mz02 = M02(3);
    df = 0;
    N_ti = 5;
    
    % T2 = 45;
    N = 192;
    %M = zeros(3,N);
    %phi = ([1:N]/N-0.5 ) * 4*pi;
    %phi = zeros(1,N);
    TR = 17.1;
    %TE = TR/2;
    alpha = 5 / 180 * pi;
    dt = TR;
    E1 = exp(-dt / T1);
    
    E10 = exp(-TI_array(1)/T1);
    
    
    Mz_recov_array1 = zeros(1, N);
    Mz_mgre_array1 = zeros(1, N);
    
    Mz1_recov1 = 1 - E10 + Mz01*E10;
    Mz1_mgre1 = 1 - E10 + Mz01*E10;
    
    Mz_recov_total1 = zeros(1, N*N_ti);
    Mz_mgre_total1 = zeros(1, N*N_ti);
    
    
    Mz_recov_array2 = zeros(1, N);
    Mz_mgre_array2 = zeros(1, N);
    
    Mz1_recov2 = 1 - E10 + Mz02*E10;
    Mz1_mgre2 = 1 - E10 + Mz02*E10;
    
    Mz_recov_total2 = zeros(1, N*N_ti);
    Mz_mgre_total2 = zeros(1, N*N_ti);
    for tt = 1:N_ti
        for i = 1:N
            if i == 1
                Mz_recov_array1(i) = Mz1_recov1;
            else
                Mz_recov_array1(i) = 1 - E1 + Mz_recov_array1(i-1) * E1;
            end
        end
        Mz1_recov1 = cos(flip1)*Mz_recov_array1(end);
        Mz_recov_total1((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array1;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array1(i) = Mz1_mgre1;
            else
                Mz_mgre_array1(i) = 1 - E1 + Mz_mgre_array1(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre1 = cos(flip1)*Mz_mgre_array1(end);
        Mz_mgre_total1((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array1;
        
        
        for i = 1:N
            if i == 1
                Mz_recov_array2(i) = Mz1_recov2;
            else
                Mz_recov_array2(i) = 1 - E1 + Mz_recov_array2(i-1) * E1;
            end
        end
        Mz1_recov2 = cos(flip2)*Mz_recov_array2(end);
        Mz_recov_total2((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array2;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array2(i) = Mz1_mgre2;
            else
                Mz_mgre_array2(i) = 1 - E1 + Mz_mgre_array2(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre2 = cos(flip2)*Mz_mgre_array2(end);
        Mz_mgre_total2((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array2;
    end
%     
%     t = (1:N* N_ti) * TR;
%     figure();
%     plot(t/1000, Mz_mgre_total1,'LineWidth', 2);
%     hold on;
%     plot(t/1000, Mz_mgre_total2, 'LineWidth', 2);
%     legend({'SR', 'IR'}, 'Location', 'SouthEast');
%     grid on;
%     ylabel('M_z'); xlabel('Time (s)');
    
    %
    N = 96;
    N_ti = 10;
    M03 = [0 0 cos(flip1)];
    Mz03 = M03(3);
    Mz_recov_array3 = zeros(1, N);
    Mz_mgre_array3 = zeros(1, N);
    
    Mz1_recov3 = 1 - E10 + Mz03*E10;
    Mz1_mgre3 = 1 - E10 + Mz03*E10;
    
    Mz_recov_total3 = zeros(1, N*N_ti);
    Mz_mgre_total3 = zeros(1, N*N_ti);
    for tt = 1:N_ti
        
        for i = 1:N
            if i == 1
                Mz_recov_array3(i) = Mz1_recov3;
            else
                Mz_recov_array3(i) = 1 - E1 + Mz_recov_array3(i-1) * E1;
            end
        end
        Mz1_recov3 = cos(flip1)*Mz_recov_array3(end);
        Mz_recov_total3((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array3;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array3(i) = Mz1_mgre3;
            else
                Mz_mgre_array3(i) = 1 - E1 + Mz_mgre_array3(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre3 = cos(flip1)*Mz_mgre_array3(end);
        Mz_mgre_total3((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array3;
    end
    
    t = (1:N* N_ti) * TR;
    subplot(2,2,ttt);
    plot(t/1000, Mz_mgre_total1,'LineWidth', 3);
    hold on;
    plot(t/1000, Mz_mgre_total2, 'LineWidth', 2);
    plot(t/1000, Mz_mgre_total3, 'LineWidth', 1.5);
    legend({'SR 1-shot', 'IR', 'SR 2-shot'}, 'Location', 'SouthEast');
    grid on;
    ylabel('M_z'); xlabel('Time (s)');
    title(['T1 = ', num2str(T1), ' ms']);
end

%% Main Body 2 (6 echo)
figure('Position', [100 100 1600 1000]);
T1_array = [400, 500, 600, 700];
for ttt = 1:length(T1_array)
    T1 = T1_array(ttt);
    TI_array = 10.5 + TR * lSegment * ((1:5) - 1);
    flip1 = pi/2;
    flip2 = pi;
    M01 = [0 0 cos(flip1)];
    Mz01 = M01(3);
    M02 = [0 0 cos(flip2)];
    Mz02 = M02(3);
    df = 0;
    
    % T2 = 45;
    N = 192;
    %M = zeros(3,N);
    %phi = ([1:N]/N-0.5 ) * 4*pi;
    %phi = zeros(1,N);
    TR = 11.3;
    %TE = TR/2;
    alpha = 5 / 180 * pi;
    dt = TR;
    E1 = exp(-dt / T1);
    E10 = exp(-TI_array(1)/T1);
    N_ti = length(TI_array);
    
    Mz_recov_array1 = zeros(1, N);
    Mz_mgre_array1 = zeros(1, N);
    
    Mz1_recov1 = 1 - E10 + Mz01*E10;
    Mz1_mgre1 = 1 - E10 + Mz01*E10;
    
    Mz_recov_total1 = zeros(1, N*N_ti);
    Mz_mgre_total1 = zeros(1, N*N_ti);
    
    
    Mz_recov_array2 = zeros(1, N);
    Mz_mgre_array2 = zeros(1, N);
    
    Mz1_recov2 = 1 - E10 + Mz02*E10;
    Mz1_mgre2 = 1 - E10 + Mz02*E10;
    
    Mz_recov_total2 = zeros(1, N*N_ti);
    Mz_mgre_total2 = zeros(1, N*N_ti);
    for tt = 1:N_ti
        
        for i = 1:N
            if i == 1
                Mz_recov_array1(i) = Mz1_recov1;
            else
                Mz_recov_array1(i) = 1 - E1 + Mz_recov_array1(i-1) * E1;
            end
        end
        Mz1_recov1 = cos(flip1)*Mz_recov_array1(end);
        Mz_recov_total1((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array1;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array1(i) = Mz1_mgre1;
            else
                Mz_mgre_array1(i) = 1 - E1 + Mz_mgre_array1(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre1 = cos(flip1)*Mz_mgre_array1(end);
        Mz_mgre_total1((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array1;
        
        
        for i = 1:N
            if i == 1
                Mz_recov_array2(i) = Mz1_recov2;
            else
                Mz_recov_array2(i) = 1 - E1 + Mz_recov_array2(i-1) * E1;
            end
        end
        Mz1_recov2 = cos(flip2)*Mz_recov_array2(end);
        Mz_recov_total2((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array2;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array2(i) = Mz1_mgre2;
            else
                Mz_mgre_array2(i) = 1 - E1 + Mz_mgre_array2(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre2 = cos(flip2)*Mz_mgre_array2(end);
        Mz_mgre_total2((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array2;
    end
%     
%     t = (1:N* N_ti) * TR;
%     figure();
%     plot(t/1000, Mz_mgre_total1,'LineWidth', 2);
%     hold on;
%     plot(t/1000, Mz_mgre_total2, 'LineWidth', 2);
%     legend({'SR', 'IR'}, 'Location', 'SouthEast');
%     grid on;
%     ylabel('M_z'); xlabel('Time (s)');
    
    %
    N = 96;
    N_ti = 10;
    M03 = [0 0 cos(flip1)];
    Mz03 = M03(3);
    Mz_recov_array3 = zeros(1, N);
    Mz_mgre_array3 = zeros(1, N);
    
    Mz1_recov3 = 1 - E10 + Mz03*E10;
    Mz1_mgre3 = 1 - E10 + Mz03*E10;
    
    Mz_recov_total3 = zeros(1, N*N_ti);
    Mz_mgre_total3 = zeros(1, N*N_ti);
    for tt = 1:N_ti
        
        for i = 1:N
            if i == 1
                Mz_recov_array3(i) = Mz1_recov3;
            else
                Mz_recov_array3(i) = 1 - E1 + Mz_recov_array3(i-1) * E1;
            end
        end
        Mz1_recov3 = cos(flip1)*Mz_recov_array3(end);
        Mz_recov_total3((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_recov_array3;
        
        for i = 1:N
            if i == 1
                Mz_mgre_array3(i) = Mz1_mgre3;
            else
                Mz_mgre_array3(i) = 1 - E1 + Mz_mgre_array3(i-1) * cos(alpha) * E1;
            end
        end
        Mz1_mgre3 = cos(flip1)*Mz_mgre_array3(end);
        Mz_mgre_total3((1+(tt-1)*N):(N+(tt-1)*N)) = Mz_mgre_array3;
    end
    
    t = (1:N* N_ti) * TR;
    subplot(2,2,ttt);
    plot(t/1000, Mz_mgre_total1,'LineWidth', 3);
    hold on;
    plot(t/1000, Mz_mgre_total2, 'LineWidth', 2);
    plot(t/1000, Mz_mgre_total3, 'LineWidth', 1.5);
    legend({'SR 1-shot', 'IR', 'SR 2-shot'}, 'Location', 'SouthEast');
    grid on;
    ylabel('M_z'); xlabel('Time (s)');
    title(['T1 = ', num2str(T1), ' ms']);
end

%% A more precise simulation
df = 0;
N_ti = 5;
N = 192;
T1_array = 300:10:800;

flip_prep_array = [pi/2, pi];
Nshot_array = [1, 2];
TI_array = [10:10:1000];
TR = 17.1;
alpha = 5 / 180 * pi;
dt = 0.1;
T2 = 45;
clear Mz_recov_TE_dict;

for tff = 1:length(flip_prep_array)
    flip_prep = flip_prep_array(tff);
    for tnn = 1:length(Nshot_array)
        Nshot = Nshot_array(tnn);
        for t11 = 1:length(T1_array)
            T1 = T1_array(t11);
            for tii = 1:length(TI_array)
                TI = TI_array(tii);
                disp(['T1: ', num2str(T1), ' ms']);
                disp(['TI: ', num2str(TI), ' ms']);
                %tic;
                total_time = (TI*Nshot + TR * N) * N_ti;
                t = [dt:dt:total_time];
                
                Nstep_TI = TI / dt;
                Nstep_mgre = TR * N / dt / Nshot;
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
                        
                        Mz1 = Mz_TI_recov_array(end);
                        Mxy1 = Mxy_TI_recov_array(end);
                        Mz_mgre_recov_array = zeros(1, Nstep_mgre);
                        Mxy_mgre_recov_array = zeros(1, Nstep_mgre);
                        Nstep_TR = TR/dt;
                        for i = 1:Nstep_mgre
                            
                            if mod(i-1, Nstep_TR) == 0
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
                        Mz0 = cos(flip_prep)*Mz_mgre_recov_array(end);
                        Mxy0 = cos(flip_prep)*Mxy_mgre_recov_array(end);
                        Mz_recov_array = [Mz_recov_array, Mz_TI_recov_array, Mz_mgre_recov_array];
                        Mxy_recov_array = [Mxy_recov_array, Mxy_TI_recov_array, Mxy_mgre_recov_array];
                    end
                end
                TE_array = [1.4, 3.4, 5.4, 7.4, 9.4, 11.4, 13.4, 15.4];
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
                Mz_recov_TE_dict(:,:,tii,t11,tnn,tff) = Mz_recov_TE;
            end
        end
    end
end
%%
figure(); plot(Mz_recov_array, 'LineWidth', 1.5); hold on;
for i = 1:Nshot*Nprep
    x = i*Nstep_TI + (i-1) * Nstep_mgre;
    plot([x x], [-1 1], 'k', 'LineWidth', 0.5);
    x_init = (i-1)*Nstep_TI + (i-1) * Nstep_mgre;
    plot([x_init x_init], [-1 1], 'k', 'LineWidth', 0.5);
end

figure(); plot(Mxy_recov_array, 'LineWidth', 1.5);
%% TE dict generation (Not necessary)
TE_array = [1.4, 3.4, 5.4, 7.4, 9.4, 11.4, 13.4, 15.4];
Nstep_TE_array = round(TE_array / dt);
Nstep_TR = round(TR/dt);
clear Nstep_TE_indices_dict;
for tff = 1:length(flip_prep_array)
    flip_prep = flip_prep_array(tff);
    for tnn = 1:length(Nshot_array)
        Nshot = Nshot_array(tnn);
        for tii = 1:length(TI_array)
            TI = TI_array(tii);
            Nstep_TI = round(TI/dt);
            Nstep_TE_indices = [];
            for i = 1:Nprep
                for k = 1:Nshot
                    for j = 1:(N/Nshot)
                        Nstep_TE_indices = [Nstep_TE_indices, (Nshot*(i-1)+k-1)*(N/Nshot)*Nstep_TR + (Nshot*(i-1)+k)*Nstep_TI + (j-1)*Nstep_TR + Nstep_TE_array];
                    end
                end
            end
            Nstep_TE_indices_dict(:,:,tii,tnn) = Nstep_TE_indices;
        end
    end
end
%% 
figure();, plot(Mz_recov_TE_dict(:,:,100,26,1,1));ylim([-1 1]);
SumSQRT = @(x,y) sum(sqrt((x-y).^2));
%Mz_recov_TE = Mz_recov_array(Nstep_TE_indices);
%figure(); plot(t(Nstep_TE_indices)/1000, Mz_recov_TE + 0.1*randn(1,length(Mz_recov_TE)), 'LineWidth', 1.5);

clear Mz_recov_TE_dict_noise;
T1_array2 = [400, 500, 600, 700];
for tff = 1:length(flip_prep_array)
    flip_prep = flip_prep_array(tff);
    for tnn = 1:length(Nshot_array)
        Nshot = Nshot_array(tnn);
        for t11 = 1:length(T1_array2)
            T1 = T1_array2(t11);
            idx = find(T1 == T1_array);
            for tii = 1:length(TI_array)
                TI = TI_array(tii);
                Mz_recov_TE_dict_noise(:,:,tii,t11,tnn,tff) = Mz_recov_TE_dict(:,:,tii,idx,tnn,tff) + 0.1*randn(1,length(Nstep_TE_indices));
            end
        end
    end
end
%%
t11 = 4;
Mz_recov_TE_noise = Mz_recov_TE_dict_noise(:,:,:,t11,:,:);
RMS_mat = zeros(length(TI_array), length(T1_array), length(Nshot_array), length(flip_prep_array));
for tff = 1:length(flip_prep_array)
    flip_prep = flip_prep_array(tff);
    for tnn = 1:length(Nshot_array)
        Nshot = Nshot_array(tnn);
        for tii = 1:length(TI_array)
            TI = TI_array(tii);
            for t11 = 1:length(T1_array)
                T1 = T1_array(t11);
                RMS_mat(tii, t11, tnn, tff) = SumSQRT(Mz_recov_TE_noise(:,:,tii,1,tnn,tff), Mz_recov_TE_dict(:,:,tii,t11,tnn,tff));         
            end
        end
    end
end

RMS_min = min(RMS_mat(:));
%%
idx = find(RMS_mat == RMS_min);
[r, c, v, l] = ind2sub(size(RMS_mat),idx);
figure();, plot(Mz_recov_TE_dict_noise(:,:,r,1,v,l));ylim([-1 1]);
figure();, plot(Mz_recov_TE_dict(:,:,r,c,v,l));ylim([-1 1]);
% r is TI
% c is T1
% v is Nshot 1 - 1shot, 2 - 2shots
% l is flip_prep: 1 - SR, 2 - IR
%y = awgn(Mz_recov_TE,10,'measured');

%figure(); plot(t(Nstep_TE_indices),Mxy_recov_array(Nstep_TE_indices), 'LineWidth', 1.5);
%one_rep_Mz(1:N) = Mz_recov_array;

%% After Discuss with Randy
% No need to iterate TI, but change the N (#ofSegment)
df = 0;
N_ti = 5;
N_array = [64:2:256];
T1_array = 300:2:800;

flip_prep_array = [pi/2, pi];
Nshot = 1;
TI = 10.5;
TR = 17.1;
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
            Nstep_mgre = round(TR * N / dt / Nshot);
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
                    
                    Mz1 = Mz_TI_recov_array(end);
                    Mxy1 = Mxy_TI_recov_array(end);
                    Mz_mgre_recov_array = zeros(1, Nstep_mgre);
                    Mxy_mgre_recov_array = zeros(1, Nstep_mgre);
                    Nstep_TR = TR/dt;
                    for i = 1:Nstep_mgre
                        
                        if mod(i-1, Nstep_TR) == 0
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
                    Mz0 = cos(flip_prep)*Mz_mgre_recov_array(end);
                    Mxy0 = cos(flip_prep)*Mxy_mgre_recov_array(end);
                    Mz_recov_array = [Mz_recov_array, Mz_TI_recov_array, Mz_mgre_recov_array];
                    Mxy_recov_array = [Mxy_recov_array, Mxy_TI_recov_array, Mxy_mgre_recov_array];
                end
            end
            % TE_array = [1.4, 3.4, 5.4, 7.4, 9.4, 11.4, 13.4, 15.4];
            TE_array = [1.4]; % Only looking at first-echo
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

%% Adding noise 100 times
% figure(); plot(Mz_recov_TE_dict{1, 1, 1});
iters = 100;
Mz_recov_TE_dict_noise = cell(length(T1_array), length(N_array), length(flip_prep_array), iters);
T1_array2 = [400, 500, 600, 700];
for its = 1:iters
    for n = 1:length(N_array)
        N = N_array(n);
        for tff = 1:length(flip_prep_array)
            flip_prep = flip_prep_array(tff);
            for t11 = 1:length(T1_array2)
                T1 = T1_array2(t11);
                idx = find(T1 == T1_array);
                Mz_recov_TE_dict_noise{t11, n, tff, its} = Mz_recov_TE_dict{idx, n, tff} + 0.1*randn(1,length(Mz_recov_TE_dict{idx, n, tff}));
            end
        end
    end
end
%% Plot
n = 65;
N = N_array(n);
total_time = (TI*Nshot + TR * N) * N_ti;
t = [dt:dt:total_time];
figure();
plot(Mz_recov_TE_dict_noise{1, n, 1, 1});
hold on;
plot(Mz_recov_TE_dict{51, n, 1}, 'LineWidth', 1.5);
ylim([-1.2 1.2]);
%% Plot abs value
n = 65;
N = N_array(n);
total_time = (TI*Nshot + TR * N) * N_ti;
t = [dt:dt:total_time];
figure();
plot(abs(Mz_recov_TE_dict_noise{4, n, 2, 1}(end-N+1:end)));
hold on;
plot(abs(Mz_recov_TE_dict{201, n, 2}(end-N+1:end)), 'LineWidth', 1.5);
ylim([0 1.2]);
%% Plot inverted abs value
n = 65;
N = N_array(n);
total_time = (TI*Nshot + TR * N) * N_ti;
t = [dt:dt:total_time];
figure();
Mz_recov_TE_noise_abs = abs(Mz_recov_TE_dict_noise{4, n, 2, 1}(end-N+1:end));
idx = find(Mz_recov_TE_noise_abs == min(Mz_recov_TE_noise_abs));
Mz_recov_TE_noise_abs(1:idx) = -Mz_recov_TE_noise_abs(1:idx);
plot(Mz_recov_TE_noise_abs);
hold on;
plot(Mz_recov_TE_dict{201, n, 2}(end-N+1:end), 'LineWidth', 1.5);
ylim([-1.2 1.2]);
%% Dictionary Matching
SumSQRT = @(x,y) sum(sqrt((x-y).^2))/length(x);
figure();
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
ytlbl = 64:16:256;
yt = 1:8:length(64:2:256);
xt = [1,2];
xtlbl = {'SR', 'IR'};
for t111 = 1:length(T1_array2)
    subplot(2,2,t111);
    imagesc(squeeze(error_mat(t111,:,:,1))); colorbar; caxis([0 10]);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl);
end

figure();
for t111 = 1:length(T1_array2)
    subplot(2,2,t111);
    imagesc(squeeze(error_mat(t111,:,:,2))); colorbar; caxis([0 30]);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl);
end

%% Plot curves to find optimal
figure();
xtlbl = 64:32:256;
xt = 1:16:length(64:2:256);
for t111 = 1:length(T1_array2)
    subplot(2,2,t111);
    plot(error_mat(t111,:,2,1)); ylim([0 30]); hold on;
    y = medfilt1(error_mat(t111,:,2,1), 20);
    plot(y);
    ylabel('Absolute Error (ms)'); xlabel('# of Segments');
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
    title(['T1 = ', num2str(T1_array2(t111)), ' ms']);
end

figure();
xtlbl = 64:32:256;
xt = 1:16:length(64:2:256);
seg_array = 64:2:256;
idx_array = 1:length(seg_array);
for t111 = 1:length(T1_array2)
    subplot(2,2,t111);
    plot(error_mat(t111,:,2,2)); ylim([0 35]); hold on; %ylim([0 0.05]);% 
    y_temp = [repmat(error_mat(t111,1,2,2), [1, 10]), error_mat(t111,:,2,2)];
    y_med = medfilt1(y_temp, 10);
    y_med = y_med(11:end);
    plot(y_med, 'LineWidth', 1.5);
    y = mean(error_mat(t111,end-33:end,2,2));
    y_sd = std(error_mat(t111,end-33:end,2,2));
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
%% Save mat
save_results = ['IRorSR_SimPhantom.mat']
%clear twix.obj;
%cd(fid_path);
save_dir = GetFullPath(cat(2, mainpath, '/../../../Data/Results/Simulation/'));
if ~exist(save_dir, 'dir')
   mkdir(save_dir); 
end
save(cat(2, save_dir, save_results), 'Mz_recov_TE_dict', 'Mz_recov_TE_dict_noise', 'T1_array', 'N_array', 'flip_prep_array', 'its', 'RMS_mat', 'best_fit_mat',...
    'best_fit_mat_stats', 'error_mat');
%[Atr,Btr] = freeprecess(TR,T1,T2,df);
%[Ate,Bte] = freeprecess(TE,T1,T2,df);
%Rflip = yrot(alpha);
% for k = 1:N
%     [M1sig,M1] = gssignal(alpha,T1,T2,TE,TR,df,phi(k));
%     M(:,k)=M1;
% end
%Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);

%[t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, prep, Mz0, df, RAMP_DOWN)