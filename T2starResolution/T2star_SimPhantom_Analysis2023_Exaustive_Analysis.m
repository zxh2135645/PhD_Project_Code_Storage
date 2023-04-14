clear all;
close all;
% Load SimPhantom_2023_analysis.mat
addpath('../function/');

base_dir = uigetdir;
fname = 'SimPhantom_2023_analysis.mat';
f_to_read = cat(2, base_dir, '/', fname);
load(f_to_read);

SNR_mat_04 = SimPhantom_2023_analysis.SNR_mat_04;
SNR_mat_06 = SimPhantom_2023_analysis.SNR_mat_06;
SNR_mat_08 = SimPhantom_2023_analysis.SNR_mat_08;
SNR_mat_10 = SimPhantom_2023_analysis.SNR_mat_10;
SNR_mat_20 = SimPhantom_2023_analysis.SNR_mat_20;
SNR_mat_30 = SimPhantom_2023_analysis.SNR_mat_30;
res_array = SimPhantom_2023_analysis.res_array;
sigma_array = SimPhantom_2023_analysis.sigma_array;
width_max = SimPhantom_2023_analysis.width_max;
accu_matt2_04 = SimPhantom_2023_analysis.accu_matt2_04;
accu_matt2_06 = SimPhantom_2023_analysis.accu_matt2_06;
accu_matt2_08 = SimPhantom_2023_analysis.accu_matt2_08;
accu_matt2_10 = SimPhantom_2023_analysis.accu_matt2_10;
accu_matt2_20 = SimPhantom_2023_analysis.accu_matt2_20;
accu_matt2_30 = SimPhantom_2023_analysis.accu_matt2_30;

auc_matt2_04 = SimPhantom_2023_analysis.auc_matt2_04;
auc_matt2_06 = SimPhantom_2023_analysis.auc_matt2_06;
auc_matt2_08 = SimPhantom_2023_analysis.auc_matt2_08;
auc_matt2_10 = SimPhantom_2023_analysis.auc_matt2_10;
auc_matt2_20 = SimPhantom_2023_analysis.auc_matt2_20;
auc_matt2_30 = SimPhantom_2023_analysis.auc_matt2_30;

dice_matt2_04 = SimPhantom_2023_analysis.dice_matt2_04;
dice_matt2_06 = SimPhantom_2023_analysis.dice_matt2_06;
dice_matt2_08 = SimPhantom_2023_analysis.dice_matt2_08;
dice_matt2_10 = SimPhantom_2023_analysis.dice_matt2_10;
dice_matt2_20 = SimPhantom_2023_analysis.dice_matt2_20;
dice_matt2_30 = SimPhantom_2023_analysis.dice_matt2_30;

sens_matt2_04 = SimPhantom_2023_analysis.sens_matt2_04;
sens_matt2_06 = SimPhantom_2023_analysis.sens_matt2_06;
sens_matt2_08 = SimPhantom_2023_analysis.sens_matt2_08;
sens_matt2_10 = SimPhantom_2023_analysis.sens_matt2_10;
sens_matt2_20 = SimPhantom_2023_analysis.sens_matt2_20;
sens_matt2_30 = SimPhantom_2023_analysis.sens_matt2_30;

accu_matt3_04 = SimPhantom_2023_analysis.accu_matt3_04;
accu_matt3_06 = SimPhantom_2023_analysis.accu_matt3_06;
accu_matt3_08 = SimPhantom_2023_analysis.accu_matt3_08;
accu_matt3_10 = SimPhantom_2023_analysis.accu_matt3_10;
accu_matt3_20 = SimPhantom_2023_analysis.accu_matt3_20;
accu_matt3_30 = SimPhantom_2023_analysis.accu_matt3_30;

auc_matt3_04 = SimPhantom_2023_analysis.auc_matt3_04;
auc_matt3_06 = SimPhantom_2023_analysis.auc_matt3_06;
auc_matt3_08 = SimPhantom_2023_analysis.auc_matt3_08;
auc_matt3_10 = SimPhantom_2023_analysis.auc_matt3_10;
auc_matt3_20 = SimPhantom_2023_analysis.auc_matt3_20;
auc_matt3_30 = SimPhantom_2023_analysis.auc_matt3_30;
%%
auc_matt3_avg = (auc_matt3_30+auc_matt3_20+auc_matt3_10+auc_matt3_08+auc_matt3_06+auc_matt3_04)/6;
accu_matt3_avg = (accu_matt3_30+accu_matt3_20+accu_matt3_10+accu_matt3_08+accu_matt3_06+accu_matt3_04)/6;
sens_matt2_avg = (sens_matt2_30+sens_matt2_20+sens_matt2_10+sens_matt2_08+sens_matt2_06+sens_matt2_04)/6;

figure();
imagesc(sens_matt2_avg.');axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

%%
noise_level = 10;
figure();
plot(SNR_mat_06(:,noise_level), accu_matt3_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt3_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt3_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt3_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt3_30(:,noise_level), 'LineWidth', 2);
%%
noise_level = 10;
figure();
plot(SNR_mat_04(:,noise_level), accu_matt3_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_08(:,noise_level), accu_matt3_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_10(:,noise_level), accu_matt3_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_20(:,noise_level), accu_matt3_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_30(:,noise_level), accu_matt3_30(:,noise_level), 'LineWidth', 2);

%%
noise_level = 10;
figure('Position', [0 100 900 300]);
p1 = plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level), 'LineWidth', 2);
p3 = plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level), 'LineWidth', 2);
p4 = plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level), 'LineWidth', 2);
p5 = plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level), 'LineWidth', 2);
p6 = plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level), 'LineWidth', 2);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;
ylim([0 1]); xlim([0 25]);
axis off;

noise_level = 5;
figure('Position', [0 100 900 300]);
plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level),  'LineWidth', 2, 'Color', c1);
hold on;
plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level),  'LineWidth', 2, 'Color', c2);
plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level),  'LineWidth', 2, 'Color', c3);
plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level),  'LineWidth', 2, 'Color', c4);
plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level),  'LineWidth', 2, 'Color', c5);
plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level),  'LineWidth', 2, 'Color', c6);
ylim([0 1]); xlim([0 25]);
grid on; axis off;

%%
noise_level = 10;
figure('Position', [0 100 900 300]);
p1 = plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level), 'LineWidth', 3);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level), 'LineWidth', 3);
p3 = plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level), 'LineWidth', 3);
p4 = plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level), 'LineWidth', 3);
p5 = plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level), 'LineWidth', 3);
p6 = plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level), 'LineWidth', 3);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;
ylim([0 1]); xlim([0 25]);
grid on;

noise_level = 5;
figure('Position', [0 100 900 300]);
plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level),  'LineWidth', 3, 'Color', c1);
hold on;
plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level),  'LineWidth', 3, 'Color', c2);
plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level),  'LineWidth', 3, 'Color', c3);
plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level),  'LineWidth', 3, 'Color', c4);
plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level),  'LineWidth', 3, 'Color', c5);
plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level),  'LineWidth', 3, 'Color', c6);
ylim([0 1]); xlim([0 25]);
grid on; 

%%
noise_level = 10;
figure('Position', [0 100 900 300]);
p1 = plot(SNR_mat_04(:,noise_level), accu_matt3_04(:,noise_level), 'LineWidth', 3);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level), 'LineWidth', 3);
p3 = plot(SNR_mat_08(:,noise_level), accu_matt3_08(:,noise_level), 'LineWidth', 3);
p4 = plot(SNR_mat_10(:,noise_level), accu_matt3_10(:,noise_level), 'LineWidth', 3);
p5 = plot(SNR_mat_20(:,noise_level), accu_matt3_20(:,noise_level), 'LineWidth', 3);
p6 = plot(SNR_mat_30(:,noise_level), accu_matt3_30(:,noise_level), 'LineWidth', 3);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;
ylim([0 1]); xlim([0 25]); grid on;

noise_level = 5;
figure('Position', [0 100 900 300]);
plot(SNR_mat_04(:,noise_level), accu_matt3_04(:,noise_level),  'LineWidth', 3, 'Color', c1);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level),  'LineWidth', 3, 'Color', c2);
plot(SNR_mat_08(:,noise_level), accu_matt3_08(:,noise_level),  'LineWidth', 3, 'Color', c3);
plot(SNR_mat_10(:,noise_level), accu_matt3_10(:,noise_level),  'LineWidth', 3, 'Color', c4);
plot(SNR_mat_20(:,noise_level), accu_matt3_20(:,noise_level),  'LineWidth', 3, 'Color', c5);
plot(SNR_mat_30(:,noise_level), accu_matt3_30(:,noise_level),  'LineWidth', 3, 'Color', c6);
ylim([0 1]); xlim([0 25]);
grid on; 

%%
noise_level = 10;
figure('Position', [0 100 900 300]);
p1 = plot(SNR_mat_04(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 3);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 3);
p3 = plot(SNR_mat_08(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 3);
p4 = plot(SNR_mat_10(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 3);
p5 = plot(SNR_mat_20(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 3);
p6 = plot(SNR_mat_30(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 3);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;
ylim([0.6 1]); xlim([0 25]); grid on;

noise_level = 5;
figure('Position', [0 100 900 300]);
plot(SNR_mat_04(:,noise_level), accu_matt2_04(:,noise_level),  'LineWidth', 3, 'Color', c1);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level),  'LineWidth', 3, 'Color', c2);
plot(SNR_mat_08(:,noise_level), accu_matt2_08(:,noise_level),  'LineWidth', 3, 'Color', c3);
plot(SNR_mat_10(:,noise_level), accu_matt2_10(:,noise_level),  'LineWidth', 3, 'Color', c4);
plot(SNR_mat_20(:,noise_level), accu_matt2_20(:,noise_level),  'LineWidth', 3, 'Color', c5);
plot(SNR_mat_30(:,noise_level), accu_matt2_30(:,noise_level),  'LineWidth', 3, 'Color', c6);
ylim([0.6 1]); xlim([0 25]);
grid on; 
%%
noise_level = 10;
figure();
p1 = plot(SNR_mat_04(:,noise_level), accu_matt3_04(:,noise_level), 'LineWidth', 2);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level), 'LineWidth', 2);
p3 = plot(SNR_mat_08(:,noise_level), accu_matt3_08(:,noise_level), 'LineWidth', 2);
p4 = plot(SNR_mat_10(:,noise_level), accu_matt3_10(:,noise_level), 'LineWidth', 2);
p5 = plot(SNR_mat_20(:,noise_level), accu_matt3_20(:,noise_level), 'LineWidth', 2);
p6 = plot(SNR_mat_30(:,noise_level), accu_matt3_30(:,noise_level), 'LineWidth', 2);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;

noise_level = 5;
plot(SNR_mat_04(:,noise_level), accu_matt3_04(:,noise_level), '--', 'LineWidth', 2, 'Color', c1);
plot(SNR_mat_06(:,noise_level), accu_matt3_06(:,noise_level), '--', 'LineWidth', 2, 'Color', c2);
plot(SNR_mat_08(:,noise_level), accu_matt3_08(:,noise_level), '--', 'LineWidth', 2, 'Color', c3);
plot(SNR_mat_10(:,noise_level), accu_matt3_10(:,noise_level), '--', 'LineWidth', 2, 'Color', c4);
plot(SNR_mat_20(:,noise_level), accu_matt3_20(:,noise_level), '--', 'LineWidth', 2, 'Color', c5);
plot(SNR_mat_30(:,noise_level), accu_matt3_30(:,noise_level), '--', 'LineWidth', 2, 'Color', c6);

grid on;

%% thickness vs accuracy
width_mat = repmat(width_max, [length(res_array),1]);
accu_matt3_3d = cat(3,  accu_matt3_04, accu_matt3_06, accu_matt3_08, accu_matt3_10, accu_matt3_20, accu_matt3_30);
accu_matt3_3d = permute(accu_matt3_3d, [1,3,2]);

sens_matt2_3d = cat(3,  sens_matt2_04, sens_matt2_06, sens_matt2_08, sens_matt2_10, sens_matt2_20, sens_matt2_30);
sens_matt2_3d = permute(sens_matt2_3d, [1,3,2]);

SNR_mat_3d = cat(3,  SNR_mat_04, SNR_mat_06, SNR_mat_08, SNR_mat_10, SNR_mat_20, SNR_mat_30);
SNR_mat_3d = permute(SNR_mat_3d, [1,3,2]);

accu_matt2_3d = cat(3,  accu_matt2_04, accu_matt2_06, accu_matt2_08, accu_matt2_10, accu_matt2_20, accu_matt2_30);
accu_matt2_3d = permute(accu_matt2_3d, [1,3,2]);

auc_matt2_3d = cat(3,  auc_matt2_04, auc_matt2_06, auc_matt2_08, auc_matt2_10, auc_matt2_20, auc_matt2_30);
auc_matt2_3d = permute(auc_matt2_3d, [1,3,2]);

res_level = 6;
figure();
p1 = scatter(width_mat, accu_matt3_3d(:,:,10));
%%
figure();
p1 = plot3(width_mat, sens_matt2_3d(:,:,10), SNR_mat_3d(:,:,10), 'o');
%%
sz = 72;
figure();
p1 = scatter(width_mat(1,:), sens_matt2_3d(1,:,10), sz, SNR_mat_3d(1,:,10),'filled');
hold on;
p2 = scatter(width_mat(2,:), sens_matt2_3d(2,:,10), sz, SNR_mat_3d(2,:,10),'filled');
p3 = scatter(width_mat(3,:), sens_matt2_3d(3,:,10), sz, SNR_mat_3d(3,:,10),'filled');
p4 = scatter(width_mat(4,:), sens_matt2_3d(4,:,10), sz, SNR_mat_3d(4,:,10),'filled');
p5 = scatter(width_mat(5,:), sens_matt2_3d(5,:,10), sz, SNR_mat_3d(5,:,10),'filled');
p6 = scatter(width_mat(6,:), sens_matt2_3d(6,:,10), sz, SNR_mat_3d(6,:,10),'filled');
p7 = scatter(width_mat(7,:), sens_matt2_3d(7,:,10), sz, SNR_mat_3d(7,:,10),'filled');
caxis([5 15]);

sz = 72;
figure();
p1 = scatter(width_mat(1,:), sens_matt2_3d(1,:,5), sz, SNR_mat_3d(1,:,5),'filled');
hold on;
p2 = scatter(width_mat(2,:), sens_matt2_3d(2,:,5), sz, SNR_mat_3d(2,:,5),'filled');
p3 = scatter(width_mat(3,:), sens_matt2_3d(3,:,5), sz, SNR_mat_3d(3,:,5),'filled');
p4 = scatter(width_mat(4,:), sens_matt2_3d(4,:,5), sz, SNR_mat_3d(4,:,5),'filled');
p5 = scatter(width_mat(5,:), sens_matt2_3d(5,:,5), sz, SNR_mat_3d(5,:,5),'filled');
p6 = scatter(width_mat(6,:), sens_matt2_3d(6,:,5), sz, SNR_mat_3d(6,:,5),'filled');
p7 = scatter(width_mat(7,:), sens_matt2_3d(7,:,5), sz, SNR_mat_3d(7,:,5),'filled');
caxis([5 15]);

%%
figure();
p1 = scatter(width_mat(1,:), accu_matt2_3d(1,:,10), sz, SNR_mat_3d(1,:,10),'filled');
hold on;
p2 = scatter(width_mat(2,:), accu_matt2_3d(2,:,10), sz, SNR_mat_3d(2,:,10),'filled');
p3 = scatter(width_mat(3,:), accu_matt2_3d(3,:,10), sz, SNR_mat_3d(3,:,10),'filled');
p4 = scatter(width_mat(4,:), accu_matt2_3d(4,:,10), sz, SNR_mat_3d(4,:,10),'filled');
p5 = scatter(width_mat(5,:), accu_matt2_3d(5,:,10), sz, SNR_mat_3d(5,:,10),'filled');
p6 = scatter(width_mat(6,:), accu_matt2_3d(6,:,10), sz, SNR_mat_3d(6,:,10),'filled');
p7 = scatter(width_mat(7,:), accu_matt2_3d(7,:,10), sz, SNR_mat_3d(7,:,10),'filled');
caxis([5 15]);

figure();
p1 = scatter(width_mat(1,:), accu_matt2_3d(1,:,5), sz, SNR_mat_3d(1,:,5),'filled');
hold on;
p2 = scatter(width_mat(2,:), accu_matt2_3d(2,:,5), sz, SNR_mat_3d(2,:,5),'filled');
p3 = scatter(width_mat(3,:), accu_matt2_3d(3,:,5), sz, SNR_mat_3d(3,:,5),'filled');
p4 = scatter(width_mat(4,:), accu_matt2_3d(4,:,5), sz, SNR_mat_3d(4,:,5),'filled');
p5 = scatter(width_mat(5,:), accu_matt2_3d(5,:,5), sz, SNR_mat_3d(5,:,5),'filled');
p6 = scatter(width_mat(6,:), accu_matt2_3d(6,:,5), sz, SNR_mat_3d(6,:,5),'filled');
p7 = scatter(width_mat(7,:), accu_matt2_3d(7,:,5), sz, SNR_mat_3d(7,:,5),'filled');
caxis([5 15]);

%%
figure();
p1 = scatter(width_mat(1,:), accu_matt3_3d(1,:,10), sz, SNR_mat_3d(1,:,10),'filled');
hold on;
p2 = scatter(width_mat(2,:), accu_matt3_3d(2,:,10), sz, SNR_mat_3d(2,:,10),'filled');
p3 = scatter(width_mat(3,:), accu_matt3_3d(3,:,10), sz, SNR_mat_3d(3,:,10),'filled');
p4 = scatter(width_mat(4,:), accu_matt3_3d(4,:,10), sz, SNR_mat_3d(4,:,10),'filled');
p5 = scatter(width_mat(5,:), accu_matt3_3d(5,:,10), sz, SNR_mat_3d(5,:,10),'filled');
p6 = scatter(width_mat(6,:), accu_matt3_3d(6,:,10), sz, SNR_mat_3d(6,:,10),'filled');
p7 = scatter(width_mat(7,:), accu_matt3_3d(7,:,10), sz, SNR_mat_3d(7,:,10),'filled');
caxis([5 15]);

figure();
p1 = scatter(width_mat(1,:), accu_matt3_3d(1,:,5), sz, SNR_mat_3d(1,:,5),'filled');
hold on;
p2 = scatter(width_mat(2,:), accu_matt3_3d(2,:,5), sz, SNR_mat_3d(2,:,5),'filled');
p3 = scatter(width_mat(3,:), accu_matt3_3d(3,:,5), sz, SNR_mat_3d(3,:,5),'filled');
p4 = scatter(width_mat(4,:), accu_matt3_3d(4,:,5), sz, SNR_mat_3d(4,:,5),'filled');
p5 = scatter(width_mat(5,:), accu_matt3_3d(5,:,5), sz, SNR_mat_3d(5,:,5),'filled');
p6 = scatter(width_mat(6,:), accu_matt3_3d(6,:,5), sz, SNR_mat_3d(6,:,5),'filled');
p7 = scatter(width_mat(7,:), accu_matt3_3d(7,:,5), sz, SNR_mat_3d(7,:,5),'filled');
caxis([5 15]);

%%
noise_level = 10;
figure();
p1 = scatter(width_mat(1,:), auc_matt2_3d(1,:,noise_level), sz, SNR_mat_3d(1,:,noise_level),'filled');
hold on;
p2 = scatter(width_mat(2,:), auc_matt2_3d(2,:,noise_level), sz, SNR_mat_3d(2,:,noise_level),'filled');
p3 = scatter(width_mat(3,:), auc_matt2_3d(3,:,noise_level), sz, SNR_mat_3d(3,:,noise_level),'filled');
p4 = scatter(width_mat(4,:), auc_matt2_3d(4,:,noise_level), sz, SNR_mat_3d(4,:,noise_level),'filled');
p5 = scatter(width_mat(5,:), auc_matt2_3d(5,:,noise_level), sz, SNR_mat_3d(5,:,noise_level),'filled');
p6 = scatter(width_mat(6,:), auc_matt2_3d(6,:,noise_level), sz, SNR_mat_3d(6,:,noise_level),'filled');
p7 = scatter(width_mat(7,:), auc_matt2_3d(7,:,noise_level), sz, SNR_mat_3d(7,:,noise_level),'filled');
plot(width_mat(1,:), auc_matt2_3d(1,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(2,:), auc_matt2_3d(2,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(3,:), auc_matt2_3d(3,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(4,:), auc_matt2_3d(4,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(5,:), auc_matt2_3d(5,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(6,:), auc_matt2_3d(6,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(7,:), auc_matt2_3d(7,:,noise_level), 'k', 'LineWidth', 1.5);
caxis([5 15]);


noise_level = 5;
figure();
p1 = scatter(width_mat(1,:), auc_matt2_3d(1,:,noise_level), sz, SNR_mat_3d(1,:,noise_level),'filled');
hold on;
p2 = scatter(width_mat(2,:), auc_matt2_3d(2,:,noise_level), sz, SNR_mat_3d(2,:,noise_level),'filled');
p3 = scatter(width_mat(3,:), auc_matt2_3d(3,:,noise_level), sz, SNR_mat_3d(3,:,noise_level),'filled');
p4 = scatter(width_mat(4,:), auc_matt2_3d(4,:,noise_level), sz, SNR_mat_3d(4,:,noise_level),'filled');
p5 = scatter(width_mat(5,:), auc_matt2_3d(5,:,noise_level), sz, SNR_mat_3d(5,:,noise_level),'filled');
p6 = scatter(width_mat(6,:), auc_matt2_3d(6,:,noise_level), sz, SNR_mat_3d(6,:,noise_level),'filled');
p7 = scatter(width_mat(7,:), auc_matt2_3d(7,:,noise_level), sz, SNR_mat_3d(7,:,noise_level),'filled');
plot(width_mat(1,:), auc_matt2_3d(1,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(2,:), auc_matt2_3d(2,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(3,:), auc_matt2_3d(3,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(4,:), auc_matt2_3d(4,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(5,:), auc_matt2_3d(5,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(6,:), auc_matt2_3d(6,:,noise_level), 'k', 'LineWidth', 1.5);
plot(width_mat(7,:), auc_matt2_3d(7,:,noise_level), 'k', 'LineWidth', 1.5);
caxis([5 15]);

%% 4 Dimensional plot
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure();
s1 = surf(X, Y, squeeze(sens_matt2_3d(1,:,:)), squeeze(SNR_mat_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
%s2 = surf(X, Y, squeeze(sens_matt2_3d(2,:,:)), C, 'FaceAlpha',0.5);
%s3 = surf(X, Y, squeeze(sens_matt2_3d(3,:,:)), C, 'FaceAlpha',0.5);
%s4 = surf(X, Y, squeeze(sens_matt2_3d(4,:,:)), C, 'FaceAlpha',0.5);
%s5 = surf(X, Y, squeeze(sens_matt2_3d(5,:,:)), C, 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(sens_matt2_3d(6,:,:)), squeeze(SNR_mat_3d(6,:,:)), 'FaceAlpha',0.5);
%s7 = surf(X, Y, squeeze(sens_matt2_3d(7,:,:)), C, 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';
caxis([5 15]);

%%
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure();
s1 = surf(X, Y, squeeze(SNR_mat_3d(1,:,:)), squeeze(sens_matt2_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
s2 = surf(X, Y, squeeze(SNR_mat_3d(2,:,:)), squeeze(sens_matt2_3d(2,:,:)), 'FaceAlpha',0.5);
s3 = surf(X, Y, squeeze(SNR_mat_3d(3,:,:)), squeeze(sens_matt2_3d(3,:,:)), 'FaceAlpha',0.5);
s4 = surf(X, Y, squeeze(SNR_mat_3d(4,:,:)), squeeze(sens_matt2_3d(4,:,:)), 'FaceAlpha',0.5);
s5 = surf(X, Y, squeeze(SNR_mat_3d(5,:,:)), squeeze(sens_matt2_3d(5,:,:)), 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(SNR_mat_3d(6,:,:)), squeeze(sens_matt2_3d(6,:,:)),  'FaceAlpha',0.5);
s7 = surf(X, Y, squeeze(SNR_mat_3d(7,:,:)), squeeze(sens_matt2_3d(7,:,:)), 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';



colormap(brewermap([],'*RdYlBu'));
caxis([0.5 1]);
zlim([0 15]);


width_max_ps = [0.4, 3];
res_array_ps = [0 15];
[X,Y,Z] = meshgrid(noise_level, width_max_ps, res_array_ps);
xslice = [5 10];   
yslice = [];
zslice = [];
sens_matt2_3d_ps = zeros(size(X));
ss = slice(X,Y,Z,sens_matt2_3d_ps,xslice,yslice,zslice)
%ss(1).EdgeColor = 'none';
set(ss(1),'FaceColor','interp',...
    'FaceAlpha',0.2, ...
    'LineWidth', 1.5);
%ss(2).EdgeColor = 'none';
set(ss(2),'FaceColor','interp',...
    'FaceAlpha',0.2,...
    'LineWidth', 1.5);
%%
sz = 72;
noise_level = 10;
figure();
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,12,7));
contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 1.5,'ShowText','on');
xlim([0 3]);

sz = 72;
noise_level = 5;
figure();
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,26,7));
contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 1.5,'ShowText','on');
xlim([0 3]);

%%
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure();
s1 = surf(X, Y, squeeze(SNR_mat_3d(1,:,:)), squeeze(sens_matt2_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
s2 = surf(X, Y, squeeze(SNR_mat_3d(2,:,:)), squeeze(sens_matt2_3d(2,:,:)), 'FaceAlpha',0.5);
s3 = surf(X, Y, squeeze(SNR_mat_3d(3,:,:)), squeeze(sens_matt2_3d(3,:,:)), 'FaceAlpha',0.5);
s4 = surf(X, Y, squeeze(SNR_mat_3d(4,:,:)), squeeze(sens_matt2_3d(4,:,:)), 'FaceAlpha',0.5);
s5 = surf(X, Y, squeeze(SNR_mat_3d(5,:,:)), squeeze(sens_matt2_3d(5,:,:)), 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(SNR_mat_3d(6,:,:)), squeeze(sens_matt2_3d(6,:,:)),  'FaceAlpha',0.5);
s7 = surf(X, Y, squeeze(SNR_mat_3d(7,:,:)), squeeze(sens_matt2_3d(7,:,:)), 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';
% set(s1,'FaceColor','interp');
% set(s2,'FaceColor','interp');
% set(s3,'FaceColor','interp');
% set(s4,'FaceColor','interp');
% set(s5,'FaceColor','interp');
% set(s6,'FaceColor','interp');
% set(s7,'FaceColor','interp');

colormap(brewermap([],'*RdYlBu'));
caxis([0.5 1]);
zlim([0 15]);


width_max_ps = [0.4,3];
res_array_ps = 0:15;
noise_level_ps = [1 10];
[X,Y,Z] = meshgrid(noise_level_ps, width_max_ps, res_array_ps);
xslice = [];   
yslice = [];
zslice = [5 10];
sens_matt2_3d_ps = zeros(size(X));
ss = slice(X,Y,Z,sens_matt2_3d_ps,xslice,yslice,zslice);
%ss(1).EdgeColor = 'none';
set(ss(1),'FaceColor','interp',...
    'FaceAlpha',0.2, ...
    'LineWidth', 1.5);
%ss(2).EdgeColor = 'none';
set(ss(2),'FaceColor','interp',...
    'FaceAlpha',0.2,...
    'LineWidth', 1.5);
%% Need to interp
noise_level2 = 1:0.02:10; % [1 51 101 151 201 251 301 351 401 451]
width_max2 = 0.4:0.1:3; % [1 3 5 7 17 27]
w_idx = [1 3 5 7 17 27];
n_idx = [1 51 101 151 201 251 301 351 401 451];
[X,Y] = meshgrid(noise_level, width_max);
[Xq,Yq] = meshgrid(noise_level2, width_max2);
sz = 72;
figure();
lvl1 = 5;

for ii = 1:size(SNR_mat_3d, 1)
    Vq = interp2(X,Y,squeeze(SNR_mat_3d(ii,:,:)),Xq,Yq);
    sens_interp = interp2(X,Y,squeeze(sens_matt2_3d(ii,:,:)),Xq,Yq);
    % figure();
    % surf(Xq,Yq,Vq);
    % title('Linear Interpolation Using Finer Grid');

    [M,I] = min(abs(Vq - lvl1),[],2);
    % I(w_idx)
    % noise_level2(I(w_idx))
    % width_max2(w_idx)
    % Vq(w_idx, I(w_idx))
    % sens_interp(w_idx, I(w_idx))

    sens_array = zeros(1, length(w_idx));
    for i = 1:length(w_idx)
        tmp = I(w_idx);
        sens_array(i) = sens_interp(w_idx(i), tmp(i));
    end

    hold on;
    p1 = scatter(width_max2(w_idx), noise_level2(I(w_idx)), sz, sens_array,'filled');
end