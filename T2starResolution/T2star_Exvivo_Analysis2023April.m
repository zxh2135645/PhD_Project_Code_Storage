clear all;
close all;

% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92, 18P94_Exvivo3, 18P95, 17P73, 20P48
width_array = [2.40,      0.91,  1.21,  0.53,  0.56,  0.59,          0.92,  0.59,  1.10,  0.69];
res_array_in = [0.8, 1.0, 1.3, 1.6, 2.1];
res_array_thru = [2 4 6 8];
res_mat = res_array_thru' * (res_array_in.^2);

auc_mat = [0.50	  0.48	0.63	0.46	0.59	0.72	0.47	0.50	0.58	0.64;
           0.74	  0.67	0.77	0.50	0.59	0.71	1.00	0.46	0.64	0.65;
           0.80	  0.67	0.79	0.60	0.46	0.67	1.00	0.46	0.84	0.55;
           0.88	  0.91	0.84	0.67	0.73	0.58	1.00	0.77	0.82	0.65;
           0.73	  0.75	0.75	0.50	0.50	0.58	0.89	0.50	0.65	0.50;
           0.89	  0.74	0.83	0.50	0.69	0.70	0.47	0.60	0.64	0.65;
           0.92	  0.76	0.78	0.50	0.73	0.58	0.47	0.55	0.63	0.80;
           0.96	  0.89	0.79	0.48	0.69	0.75	0.47	0.74	0.64	0.60;
           0.92	  0.88	0.76	0.46	0.68	0.58	1.00	0.62	0.79	0.60;
           0.71	  0.69	0.75	0.50	0.50	0.67	0.89	0.50	0.55	0.50;
           0.89	  0.68	0.86	0.50	0.50	0.63	0.45	0.60	0.68	0.80;
           0.93	  0.66	0.73	0.50	0.65	0.67	0.45	0.70	0.72	0.80;
           0.93	  0.67	0.66	0.50	0.64	0.58	0.45	0.58	0.69	0.67;
           0.92	  0.64	0.65	0.50	0.64	0.50	0.47	0.60	0.64	0.65;
           0.77	  0.50	0.68	0.50	0.45	0.50	0.44	0.50	0.58	0.50;
           1.00	  0.63	0.74	0.50	0.64	0.44	0.45	0.55	0.65	0.80;
           0.92	  0.61	0.55	0.50	0.49	0.58	0.47	0.43	0.55	0.71;
           0.90	  0.57	0.67	0.50	0.59	0.50	0.47	0.50	0.64	0.67;
           0.92	  0.73	0.70	0.50	0.58	0.50	0.50	0.50	0.62	0.69;
           0.67	  0.48	0.60	0.50	0.41	0.50	0.44	0.50	0.60	0.47];

auc_mat_avg16 = [1	  0.99	0.97	0.95	0.96	0.86	0.77	0.97	0.99	0.95    0.86    0.89    0.79    0.67    0.92    0.88    0.92    0.87    0.97    0.80    0.71    0.95    0.93    0.89    0.97    0.96    0.82    0.71;
           1	  0.98 	0.97	0.97	0.74	0.70	0.92	0.93	0.80	0.77    0.75    0.70    0.67    0.73    0.82    0.72    0.75    0.63    0.51    0.59    0.62    0.67    0.61    0.60    0.53    0.53    0.42    0.55;
           1	  1 	0.92	0.97	0.88	0.87	0.85	0.94	0.97	0.94    0.88    0.88    0.82    0.82    0.82    0.80    0.79    0.85    0.75    0.69    0.70    0.59    0.67    0.65    0.61    0.64    0.57    0.61;
           1	  1 	1   	0.80	0.78	1   	0.48	0.92	0.77	0.75    0.63    0.44    0.50    0.50    0.65    0.40    0.46    0.46    0.48    0.50    0.50    0.60    0.42    0.46    0.46    0.48    0.50    0.50;
           1	  0.94	0.96	0.94	0.51	0.51	0.41	0.92	0.96	0.79    0.84    0.67    0.73    0.73    0.83    0.79    0.65    0.74    0.60    0.61    0.51    0.75    0.64    0.69    0.70    0.59    0.57    0.61;
           1	  0.82	0.89	0.86	0.85	0.66	0.58	0.90	0.96	0.94    0.98    0.69    0.81    0.85    0.89    0.85    0.90    0.85    0.81    0.77    0.82    0.83    0.88    0.80    0.86    0.82    0.72    0.83;
           1	  1 	0.95	1   	0.74	0.95	0.78	0.95	0.95	0.95    0.95    0.39    0.42    0.42    0.89    0.89    0.84    0.39    0.45    0.47    0.44    0.74    0.68    0.84    0.45    0.45    0.42    0.50;
           1	  1 	0.98   	0.90	0.89	0.82	0.65	0.95	0.97	0.76    0.87    0.88    0.65    0.60    0.78    0.81    0.76    0.64    0.60    0.60    0.60    0.68    0.62    0.51    0.50    0.55    0.58    0.50;
           1	  0.96	0.97	0.95	0.92	0.86	0.64	0.84	0.84	0.84    0.83    0.82    0.81    0.65    0.81    0.74    0.76    0.76    0.74    0.78    0.63    0.72    0.71    0.73    0.68    0.68    0.70    0.58;
           1   	  0.98	0.83	0.83	0.76	0.68	0.58	0.93	0.97	0.81    0.77    0.85    0.75    0.54    0.83    0.89    0.84    0.65    0.73    0.71    0.68    0.81    0.80    0.78    0.71    0.66    0.73    0.76;];

SNR_mat = [0.31	  3.37	4.11	4.14	5.28	4.54	3.65	4.61	3.85	5.03;
           4.71	  6.25	7.50	6.32	7.93	7.57	6.20	7.56	6.87	7.07;
           6.56	  6.41	8.52	10.44	8.25	9.56	9.31	9.32	11.47	9.22;
           9.01	  10.19	11.76	9.78	14.62	13.12	12.31	12.03	10.39	13.74;
           11.88  13.62	18.26	15.15	13.63	18.55	14.51	16.03	13.90	15.98;
           5.79	  6.50	8.49	6.50	8.56	10.43	8.57	8.22	8.87	7.98;
           10.24  10.78	11.22	12.28	14.18	14.32	10.58	13.23	11.78	15.13;
           13.54  10.91	16.29	13.45	19.67	17.90	12.60	16.97	14.36	14.74;
           14.85  14.71	19.02	14.30	19.47	22.11	21.61	22.66	15.52	17.13;
           21.76  18.21	28.13	20.13	25.55	25.58	27.07	28.56	27.15	19.72;
           7.86	  9.42	11.14	10.86	11.89	15.13	11.65	11.57	12.77	13.46;
           12.28  13.42	16.94	17.85	16.46	19.01	16.56	16.20	17.10	17.88;
           15.31  14.98	21.75	17.17	22.00	25.92	20.78	20.29	20.93	19.90;
           19.98  20.40	26.01	15.43	25.84	24.81	22.13	23.48	19.57	20.60;
           19.22  15.99	34.14	23.29	41.57	22.51	28.23	23.83	43.96	25.30;
           11.04  11.80	16.48	12.43	17.84	15.95	14.82	13.02	16.42	14.39;
           13.34  15.81	19.82	21.83	20.62	20.67	17.19	13.61	17.83	20.13;
           17.71  15.19	23.27	18.53	29.56	26.02	25.83	13.82	22.57	24.63;
           20.20  20.68	28.96	18.21	36.49	24.73	30.66	14.05	19.79	27.10;
           26.67  23.17	29.97	25.62	50.76	27.09	30.36	16.22	33.85	21.99]; 


sens_mat = [0.00	0.00	0.08	0.00	0.00	0.64	0.00	0.00	0.00	0.00;
            0.38	0.00	0.25	0.00	0.25	0.71	0.00	0.00	0.30	0.20;
            0.69	0.17	0.33	0.33	0.00	0.71	0.00	0.00	1.00	0.10;
            0.77	0.83	0.75	0.33	0.50	0.64	1.00	0.60	0.95	0.30;
            0.46	0.50	0.67	0.00	0.00	0.64	1.00	0.00	0.40	0.00;
            0.54	0.17	0.17	0.00	0.25	0.64	0.00	0.00	0.15	0.20;
            0.77	0.50	0.25	0.00	0.50	0.64	0.00	0.00	0.40	0.60;
            0.92	0.67	0.67	0.00	0.25	0.79	0.00	0.60	0.90	0.20;
            0.85	0.83	0.67	0.00	0.50	0.64	1.00	0.40	0.85	0.20;
            0.69	0.50	0.67	0.00	0.00	0.71	1.00	0.00	0.35	0.00;
            0.62	0.17	0.42	0.00	0.00	0.57	0.00	0.00	0.20	0.50;
            0.85	0.33	0.50	0.00	0.00	0.64	0.00	0.00	0.70	0.70;
            0.92	0.33	0.50	0.00	0.25	0.64	0.00	0.00	0.65	0.40;
            0.92	0.67	0.67	0.00	0.50	0.57	0.00	0.20	0.40	0.30;
            0.54	0.00	0.50	0.00	0.00	0.57	0.00	0.00	0.25	0.00;
            1.00	0.17	0.58	0.00	0.25	0.50	0.00	0.00	0.15	0.50;
            0.69	0.00	0.42	0.00	0.00	0.64	0.00	0.40	0.50	0.60;
            0.85	0.17	0.50	0.00	0.25	0.57	0.00	0.00	0.55	0.50;
            0.92	0.67	0.58	0.00	0.50	0.57	0.00	0.00	0.45	0.40;
            0.62	0.17	0.42	0.00	0.00	0.57	0.00	0.00	0.20	0.00];
%%
[rho, pval] = corr(SNR_mat', auc_mat');
%% Scatter plot
width_mat = repmat(width_array, [length(res_array_in)*length(res_array_thru),1]);

sz = 72;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat(1,:), sz, auc_mat(1,:),'filled');
hold on;
p2 = scatter(width_mat(2,:), SNR_mat(2,:), sz, auc_mat(2,:),'filled');
p3 = scatter(width_mat(3,:), SNR_mat(3,:), sz, auc_mat(3,:),'filled');
p4 = scatter(width_mat(4,:), SNR_mat(4,:), sz, auc_mat(4,:),'filled');
p5 = scatter(width_mat(5,:), SNR_mat(5,:), sz, auc_mat(5,:),'filled');

p6 = scatter(width_mat(6,:), SNR_mat(6,:), sz, auc_mat(6,:),'filled');
p7 = scatter(width_mat(7,:), SNR_mat(7,:), sz, auc_mat(7,:),'filled');
p8 = scatter(width_mat(8,:), SNR_mat(8,:), sz, auc_mat(8,:),'filled');
p9 = scatter(width_mat(9,:), SNR_mat(9,:), sz, auc_mat(9,:),'filled');
p10 = scatter(width_mat(10,:), SNR_mat(10,:), sz, auc_mat(10,:),'filled');

p11 = scatter(width_mat(11,:), SNR_mat(11,:), sz, auc_mat(11,:),'filled');
p12 = scatter(width_mat(12,:), SNR_mat(12,:), sz, auc_mat(12,:),'filled');
p13 = scatter(width_mat(13,:), SNR_mat(13,:), sz, auc_mat(13,:),'filled');
p14 = scatter(width_mat(14,:), SNR_mat(14,:), sz, auc_mat(14,:),'filled');
p15 = scatter(width_mat(15,:), SNR_mat(15,:), sz, auc_mat(15,:),'filled');

p16 = scatter(width_mat(16,:), SNR_mat(16,:), sz, auc_mat(16,:),'filled');
p17 = scatter(width_mat(17,:), SNR_mat(17,:), sz, auc_mat(17,:),'filled');
p18 = scatter(width_mat(18,:), SNR_mat(18,:), sz, auc_mat(18,:),'filled');
p19 = scatter(width_mat(19,:), SNR_mat(19,:), sz, auc_mat(19,:),'filled');
p20 = scatter(width_mat(20,:), SNR_mat(20,:), sz, auc_mat(20,:),'filled');
caxis([0.4 1]);

%%
figure();
p1 = plot3(width_mat, auc_mat, SNR_mat, 'o');
%%
figure();
p1 = plot3(width_mat, sens_mat, SNR_mat, 'o');
%% Average over thru-plane
SNR_mat_thru2 = (SNR_mat(1,:) + SNR_mat(2,:) + SNR_mat(3,:) + SNR_mat(4,:) + SNR_mat(5,:))/5;
SNR_mat_thru4 = (SNR_mat(6,:) + SNR_mat(7,:) + SNR_mat(8,:) + SNR_mat(9,:) + SNR_mat(10,:))/5;
SNR_mat_thru6 = (SNR_mat(11,:) + SNR_mat(12,:) + SNR_mat(13,:) + SNR_mat(14,:) + SNR_mat(15,:))/5;
SNR_mat_thru8 = (SNR_mat(16,:) + SNR_mat(17,:) + SNR_mat(18,:) + SNR_mat(19,:) + SNR_mat(20,:))/5;

auc_mat_thru2 = (auc_mat(1,:) + auc_mat(2,:) + auc_mat(3,:) + auc_mat(4,:) + auc_mat(5,:))/5;
auc_mat_thru4 = (auc_mat(6,:) + auc_mat(7,:) + auc_mat(8,:) + auc_mat(9,:) + auc_mat(10,:))/5;
auc_mat_thru6 = (auc_mat(11,:) + auc_mat(12,:) + auc_mat(13,:) + auc_mat(14,:) + auc_mat(15,:))/5;
auc_mat_thru8 = (auc_mat(16,:) + auc_mat(17,:) + auc_mat(18,:) + auc_mat(19,:) + auc_mat(20,:))/5;

figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_thru2, sz, auc_mat_thru2,'filled');
hold on;
p2 = scatter(width_mat(1,:), SNR_mat_thru4, sz, auc_mat_thru4,'filled');
p3 = scatter(width_mat(1,:), SNR_mat_thru6, sz, auc_mat_thru6,'filled');
p4 = scatter(width_mat(1,:), SNR_mat_thru8, sz, auc_mat_thru8,'filled');

%% Average over in-plane
SNR_mat_in08 = (SNR_mat(1,:) + SNR_mat(6,:) + SNR_mat(11,:) + SNR_mat(16,:))/4;
SNR_mat_in10 = (SNR_mat(2,:) + SNR_mat(7,:) + SNR_mat(12,:) + SNR_mat(17,:))/4;
SNR_mat_in13 = (SNR_mat(3,:) + SNR_mat(8,:) + SNR_mat(13,:) + SNR_mat(18,:))/4;
SNR_mat_in16 = (SNR_mat(4,:) + SNR_mat(9,:) + SNR_mat(14,:) + SNR_mat(19,:))/4;
SNR_mat_in21 = (SNR_mat(5,:) + SNR_mat(10,:) + SNR_mat(15,:) + SNR_mat(20,:))/4;

auc_mat_in08 = (auc_mat(1,:) + auc_mat(6,:) + auc_mat(11,:) + auc_mat(16,:))/4;
auc_mat_in10 = (auc_mat(2,:) + auc_mat(7,:) + auc_mat(12,:) + auc_mat(17,:))/4;
auc_mat_in13 = (auc_mat(3,:) + auc_mat(8,:) + auc_mat(13,:) + auc_mat(18,:))/4;
auc_mat_in16 = (auc_mat(4,:) + auc_mat(9,:) + auc_mat(14,:) + auc_mat(19,:))/4;
auc_mat_in21 = (auc_mat(5,:) + auc_mat(10,:) + auc_mat(15,:) + auc_mat(20,:))/4;

figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_in08, sz, auc_mat_in08,'filled');
hold on;
p2 = scatter(width_mat(1,:), SNR_mat_in10, sz, auc_mat_in10,'filled');
p3 = scatter(width_mat(1,:), SNR_mat_in13, sz, auc_mat_in13,'filled');
p4 = scatter(width_mat(1,:), SNR_mat_in16, sz, auc_mat_in16,'filled');
p5 = scatter(width_mat(1,:), SNR_mat_in21, sz, auc_mat_in21,'filled');
caxis([0.5 1]);

figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), auc_mat_in08, sz, SNR_mat_in08,'filled');
hold on;
p2 = scatter(width_mat(1,:), auc_mat_in10, sz, SNR_mat_in10,'filled');
p3 = scatter(width_mat(1,:), auc_mat_in13, sz, SNR_mat_in13,'filled');
p4 = scatter(width_mat(1,:), auc_mat_in16, sz, SNR_mat_in16,'filled');
p5 = scatter(width_mat(1,:), auc_mat_in21, sz, SNR_mat_in21,'filled');
caxis([5 30]);
%% Width vs AUC (Averaged over spatial resolution: width vs averaged AUC)
width_array = [2.40,      0.91,  1.21,  0.53,  0.56,  0.59,    0.92,  0.59,  1.10,  0.69];
sz = 192;
figure('Position', [0 100 300 600]);
% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
auc_avg = [0.85	 0.68	0.73	0.51	0.59	0.60	0.61	0.56	0.66	0.65];
auc_best = [1.00	 0.9120	 0.8583 	0.6667	0.7292	0.7500	1.00	    0.7667	0.8375	0.8048];
% 0.8*0.8*8, 1.6*1.6*2, 0.8*0.8*6, 1.6*1.6*2, 1*1*4, 1.3*1.3*4  1,1.3,1.6*2
% and 1.6*1.6*4, 1.6*1.6*2, 1.3*1.3*2, 0.8*0.8*8

scatter(width_array, auc_avg, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylim([.5 1]);
grid on;

X = [ones(length(width_array'), 1) width_array'];
b = X \ auc_avg';
yCalc1 = X*b;
hold on;
scatter(width_array([1 8]), auc_avg([1 8]), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
plot(width_array, yCalc1, 'k', 'LineWidth', 3);

x = width_array.';
y = auc_avg';
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.8423

g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.8766
%set(gca, 'YDir','reverse')

%% Width vs AUC (Averaged over spatial resolution - width vs highest AUC)
width_array = [2.40,      0.91,  1.21,  0.53,  0.56,  0.59,    0.92,  0.59,  1.10,  0.69];
sz = 192;
figure('Position', [0 100 300 600]);
% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
%auc_avg = [0.85	 0.68	0.73	0.51	0.59	0.60	0.61	0.56	0.66	0.65];
auc_best = [1.00	 0.9120	 0.8583 	0.6667	0.7292	0.7500	1.00	    0.7667	0.8375	0.8048];
% 0.8*0.8*8, 1.6*1.6*2, 0.8*0.8*6, 1.6*1.6*2, 1*1*4, 1.3*1.3*4  1,1.3,1.6*2
% and 1.6*1.6*4, 1.6*1.6*2, 1.3*1.3*2, 0.8*0.8*8

scatter(width_array, auc_best, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylim([.5 1]);
grid on;

X = [ones(length(width_array'), 1) width_array'];
b = X \ auc_best';
yCalc1 = X*b;
hold on;
scatter(width_array([1 8]), auc_best([1 8]), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
plot(width_array, yCalc1, 'k', 'LineWidth', 3);

x = width_array.';
y = auc_best';
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.8423

g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.8766
%set(gca, 'YDir','reverse')

%% CNR vs AUC (Averaged over subjects) 2024/07/17
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
CNR = [0.271152731, 0.403216524, 0.662198867, 1.059054132, 0.275836912; 0.362920365, 0.593043431, 1.202343766, 0.7565595, 0.439118987;0.352409393, 0.903281609, 0.680035797, 0.644388966, 0.163023767; 0.581640052, 0.367325563, 0.335793349, 0.303284038, 0.225562753];
sz = 192;
figure('Position', [0 100 600 600]);
% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
auc_avg = AUC_invivo_allavg(:);
CNR = CNR(:);
% 0.8*0.8*8, 1.6*1.6*2, 0.8*0.8*6, 1.6*1.6*2, 1*1*4, 1.3*1.3*4  1,1.3,1.6*2
% and 1.6*1.6*4, 1.6*1.6*2, 1.3*1.3*2, 0.8*0.8*8

scatter(CNR, auc_avg, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylim([.5 1]);
grid on;

X = [ones(length(CNR), 1) CNR];
b = X \ auc_avg;
yCalc1 = X*b;
hold on;
scatter(CNR([10 13]), auc_avg([10 13]), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
plot(CNR, yCalc1, 'k', 'LineWidth', 3);

x = CNR;
y = auc_avg;
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.6722

g = fittype('a+b*x');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), x]\y]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.7115
%set(gca, 'YDir','reverse')


%% CNR vs AUC (Averaged over subjects) 2024/07/17
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
CNR = [0.271152731, 0.403216524, 0.662198867, 1.059054132, 0.275836912; 0.362920365, 0.593043431, 1.202343766, 0.7565595, 0.439118987;0.352409393, 0.903281609, 0.680035797, 0.644388966, 0.163023767; 0.581640052, 0.367325563, 0.335793349, 0.303284038, 0.225562753];
CNR_exvivo_avg16 = [4.917567048, 4.142151951, 3.65411296, 3.390405178, 2.930137319, 2.347969506, 1.761432212;...
    3.327271233, 3.247062677, 2.695875634, 2.757874984, 2.497432825, 2.029112228, 1.734881351;...
    2.478948297, 2.075045777, 1.775212634, 1.702296113, 1.83721189, 1.231803074, 1.249931583;...
    1.770604591, 1.45535285, 1.229800818, 1.178872845, 1.210659654, 1.074993639, 0.789078152];
AUC_exvivo_avg16 = [1, 0.96, 0.94, 0.92, 0.86, 0.8, 0.71; 0.92, 0.88, 0.88, 0.85, 0.8, 0.76, 0.69; ...
    0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.7; 0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];


sz = 192;
figure('Position', [0 100 600 600]);
% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
auc_avg = AUC_invivo_allavg(:);
CNR = CNR(:);
% 0.8*0.8*8, 1.6*1.6*2, 0.8*0.8*6, 1.6*1.6*2, 1*1*4, 1.3*1.3*4  1,1.3,1.6*2
% and 1.6*1.6*4, 1.6*1.6*2, 1.3*1.3*2, 0.8*0.8*8

scatter(CNR, auc_avg, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
hold on;
scatter(CNR_exvivo_avg16(:), AUC_exvivo_avg16(:), sz, 'filled', 'MarkerFaceColor', [87 160 211]/255);


ylim([.5 1]);
grid on;

X = [ones(length(CNR), 1) CNR];
b = X \ auc_avg;
yCalc1 = X*b;

scatter(CNR([10 13]), auc_avg([10 13]), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
scatter(CNR_exvivo_avg16(1), AUC_exvivo_avg16(1), sz, 'filled', 'MarkerFaceColor', [253,190,133]/255);
plot(CNR, yCalc1, '--k', 'LineWidth', 3);

x = CNR;
y = auc_avg;
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.6722

g = fittype('a+b*x');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), x]\y]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.7115
%set(gca, 'YDir','reverse')

xlim([0 5]);

X = [ones(length(CNR_exvivo_avg16(:)), 1) CNR_exvivo_avg16(:)];
b = X \ AUC_exvivo_avg16(:);
yCalc1 = X*b;
x = CNR_exvivo_avg16(:);
y = AUC_exvivo_avg16(:);
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.9080
plot(CNR_exvivo_avg16(:), yCalc1, '--k', 'LineWidth', 3);


x = [CNR;CNR_exvivo_avg16(:)];
y = [auc_avg;AUC_exvivo_avg16(:)];
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(0, 5, 100);

%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2)

%% Ver. 3
%% CNR vs AUC (Averaged over subjects) 2024/07/17
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
CNR = [0.271152731, 0.403216524, 0.662198867, 1.059054132, 0.275836912; 0.362920365, 0.593043431, 1.202343766, 0.7565595, 0.439118987;0.352409393, 0.903281609, 0.680035797, 0.644388966, 0.163023767; 0.581640052, 0.367325563, 0.335793349, 0.303284038, 0.225562753];
CNR_exvivo_avg16 = [4.917567048, 4.142151951, 3.65411296, 3.390405178, 2.930137319, 2.347969506, 1.761432212;...
    3.327271233, 3.247062677, 2.695875634, 2.757874984, 2.497432825, 2.029112228, 1.734881351;...
    2.478948297, 2.075045777, 1.775212634, 1.702296113, 1.83721189, 1.231803074, 1.249931583;...
    1.770604591, 1.45535285, 1.229800818, 1.178872845, 1.210659654, 1.074993639, 0.789078152];
AUC_exvivo_avg16 = [1, 0.96, 0.94, 0.92, 0.86, 0.8, 0.71; 0.92, 0.88, 0.88, 0.85, 0.8, 0.76, 0.69; ...
    0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.7; 0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];


sz = 192;
figure('Position', [0 100 400 400]);
% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
auc_avg = AUC_invivo_allavg(:);
CNR = CNR(:);
% 0.8*0.8*8, 1.6*1.6*2, 0.8*0.8*6, 1.6*1.6*2, 1*1*4, 1.3*1.3*4  1,1.3,1.6*2
% and 1.6*1.6*4, 1.6*1.6*2, 1.3*1.3*2, 0.8*0.8*8

scatter(CNR, auc_avg, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
hold on;

ylim([.5 1]);
grid on;

X = [ones(length(CNR), 1) CNR];
b = X \ auc_avg;
yCalc1 = X*b;

scatter(CNR([10 13]), auc_avg([10 13]), sz, 'filled', 'MarkerFaceColor', [253,190,133]/255);
scatter(CNR_exvivo_avg16(1), AUC_exvivo_avg16(1), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
plot(CNR, yCalc1, 'k', 'LineWidth', 3);

x = CNR;
y = auc_avg;
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.6722

g = fittype('a+b*x');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), x]\y]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.7115
%set(gca, 'YDir','reverse')

xlim([0 5]);

%% Voxel size vs AUC
auc_mat_D2 = mean(auc_mat,2);
figure();
scatter(res_mat(:), auc_mat_D2);

auc_mat_in08 = (auc_mat(1,:) + auc_mat(6,:) + auc_mat(11,:) + auc_mat(16,:))/4;
auc_mat_in10 = (auc_mat(2,:) + auc_mat(7,:) + auc_mat(12,:) + auc_mat(17,:))/4;
auc_mat_in13 = (auc_mat(3,:) + auc_mat(8,:) + auc_mat(13,:) + auc_mat(18,:))/4;
auc_mat_in16 = (auc_mat(4,:) + auc_mat(9,:) + auc_mat(14,:) + auc_mat(19,:))/4;
auc_mat_in21 = (auc_mat(5,:) + auc_mat(10,:) + auc_mat(15,:) + auc_mat(20,:))/4;

%%
figure(); scatter(res_array_in, res_mat(1,:), sz, 'filled');
hold on;
scatter(res_array_in, res_mat(2,:), sz, 'filled');
scatter(res_array_in, res_mat(3,:), sz, 'filled');
scatter(res_array_in, res_mat(4,:), sz, 'filled');

figure(); scatter(res_mat(1,:), res_array_in, sz, 'filled');
hold on;
scatter(res_mat(2,:), res_array_in, sz, 'filled');
scatter(res_mat(3,:), res_array_in, sz, 'filled');
scatter(res_mat(4,:), res_array_in, sz, 'filled');

%%
SNR_mat_D2 = mean(SNR_mat, 2);
figure(); scatter(res_array_in, auc_mat_D2(1:5), sz, SNR_mat_D2(1:5), 'filled');
hold on;
scatter(res_array_in, auc_mat_D2(6:10), sz, SNR_mat_D2(6:10), 'filled');
scatter(res_array_in, auc_mat_D2(11:15), sz, SNR_mat_D2(11:15), 'filled');
scatter(res_array_in, auc_mat_D2(16:20), sz, SNR_mat_D2(16:20), 'filled');

%% T2* value in different voxel size
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';
idx_mat = repmat(1:5, [4,1]) + 5*repmat([0:3]', [1,5]);

vol_array = vol_mat(:);
idx_array = idx_mat(:);

[B,I] = sort(vol_array)
BI = idx_array(I)

% 6, 11,12,17,18,23

AUC_avg16 = [1,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan];
AUC_invivo = [0.57,0.68,0.68,0.77,0.64;0.69,0.76,0.76,0.75,0.66;0.76,0.76,0.71,0.69,0.6;0.74,0.69,0.66,0.66,0.58];
%05/29/2023
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
AUC_invivo_subjectavg = [0.56,0.67,0.68,0.78,0.64;0.67,0.67,0.70,0.73,0.62;0.65,0.68,0.63,0.62,0.54;0.64,0.58,0.60,0.62,0.52];

AUC_avg16_array = AUC_avg16(:);
AUC_invivo_array = AUC_invivo(:);

figure(); plot(B, AUC_avg16_array(I),'o'); grid on;
figure(); plot(B, AUC_invivo_array(I),'o'); grid on;

%% Flatten the voxel size, plot vs AUC (version 2 for SCMR)
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 300 600]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);


ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
B_avg16 = B;
B_avg16(1) = 0.3 * 0.3 * 2;
plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
plot(X, Y_pred); %ylim([0.5 1])
Y_lb = modelfun(ci(:,1), X);
Y_ub = modelfun(ci(:,2), X);
plot(X, Y_lb);
plot(X, Y_ub);
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_invivo(1,:),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_invivo(2,:),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_invivo(3,:),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_invivo(4,:),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%% Flatten the voxel size, plot vs AUC (version 2 for SCMR) different size for PPT
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 300 300]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);


ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
B_avg16 = B;
B_avg16(1) = 0.3 * 0.3 * 2;
plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
plot(X, Y_pred); %ylim([0.5 1])
Y_lb = modelfun(ci(:,1), X);
Y_ub = modelfun(ci(:,2), X);
plot(X, Y_lb);
plot(X, Y_ub);
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_invivo(1,:),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_invivo(2,:),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_invivo(3,:),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_invivo(4,:),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');
%% Flatten the voxel size, plot vs AUC (Avg 16)

inplane_res =   [0.3, 0.6, 0.8, 1.0, 1.3, 1.6, 2.1];
vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_exvivo_allavg = [1, 0.96, 0.94, 0.92, 0.86, 0.8, 0.71; 0.92, 0.88, 0.88, 0.85, 0.8, 0.76, 0.69; 0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.7; 0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];
AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 400 400]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);

ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
B_avg16 = B;
B_avg16(1) = 0.3 * 0.3 * 2;
plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
plot(X, Y_pred); %ylim([0.5 1])
Y_lb = modelfun(ci(:,1), X);
Y_ub = modelfun(ci(:,2), X);
plot(X, Y_lb);
plot(X, Y_ub);
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_exvivo_allavg(1,:),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_exvivo_allavg(2,:),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_exvivo_allavg(3,:),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_exvivo_allavg(4,:),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%% Width vs AUC (Averaged over spatial resolution)
width_array = [2.40, 0.91,  1.21,  0.53,  0.56,  0.59, 0.92,  0.59,  1.10,  0.69];
sz = 192;
figure('Position', [0 100 300 600]);


% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92,    18P94_Exvivo3, 18P95,    17P73,   20P48
auc_avg = [0.8905	 0.7204	    0.8019 0.6231 0.7208  0.8354	0.7108	      0.7369	0.7641	 0.7804];


scatter(width_array, auc_avg, sz, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylim([.5 1]);
grid on;

X = [ones(length(width_array'), 1) width_array'];
b = X \ auc_avg';
yCalc1 = X*b;
hold on;
scatter(width_array([1 8]), auc_avg([1 8]), sz, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
plot(width_array, yCalc1, 'k', 'LineWidth', 3);

x = width_array.';
y = auc_avg';
Rsq2 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) % R^2 0.8423

g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(0.5, 2.5, 50);
%plot(xx,f0(xx),'--k', 'LineWidth', 3);

Rsq2 = 1 - sum((y - f0(x)).^2)/sum((y - mean(y)).^2) % R^2 0.8766

%% Flatten the voxel size, plot vs AUC from Simulation - need to load auc_array first

inplane_res =   [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0];
vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

slc = 1;
figure('Position', [100 0 800 800]);
for i = 1:size(auc_array, 3)
    subplot(2,2,i);
    AUC_exvivo_allavg = auc_array(1:length(inplane_res),:,i,slc).';
    AUC_invivo = AUC_invivo_allavg;
    Y = AUC_avg16_array(I);
    X = B;
    tbl = table(X, Y);
    modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
    %modelfun = @(b,x) b(1) + b(2)*x.^b(3);
    beta0 = [0 0 0];
    mdl = fitnlm(tbl,modelfun,beta0);

    ci = coefCI(mdl);
    b = mdl.Coefficients.Estimate;

    Y_pred = modelfun(b, X)
    %B_avg16 = B;
    %B_avg16(1) = 0.3 * 0.3 * 2;
    %plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
    hold on;
    plot(X, Y_pred); %ylim([0.5 1])
    Y_lb = modelfun(ci(:,1), X);
    Y_ub = modelfun(ci(:,2), X);
    plot(X, Y_lb);
    plot(X, Y_ub);
    %set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    %    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
    %set(plotHandles_auc(:,1), 'Visible','off');
    plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_exvivo_allavg(1,:),'o');
    set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
        'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
    plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_exvivo_allavg(2,:),'square');
    set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
        'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
    plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_exvivo_allavg(3,:),'diamond');
    set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
        'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
    plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_exvivo_allavg(4,:),'^');
    set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
        'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
    xlim([-1, 36]);
    ylim([0.5 1.0]);
    %set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
    set(gca, 'XTickLabels', []);
    set(gca, 'XTick',[0 10 20 30 36]);
    set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
    set(gca, 'YTickLabels', []);
    set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
    set(gca,'TickDir','out', 'YGrid', 'on');
    set(gca,'box','off');
end