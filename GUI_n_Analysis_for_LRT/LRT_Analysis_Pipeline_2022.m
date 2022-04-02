clear all;
close all;
%%
t = [1,4,6,8,10,12,14,16,25,27,30,31]; % min
t1 = [276.4, 443.3, 482.8, 520.0, 549.0, 571.1, 581.2, 594.7, 644.1, 648.6, 669.8, 670.1];
t1_dict_fit = [451.2, 544.6, 526.7, 590.6, 591.8, 602.1, 610.5, 618.7, 616.4, 634.8, 662.9, 643.8, 640.0, 674.6, 714.0];
t1_dict_fit2 = [426.9, 518.4, 493.2, 554.8, 551.9, 564.9, 574.8, 592.3, 704.1, 594.4, 624.6, 659.6, 608.8, 714.0, 704.4];
t1_molli_fit = [329.609311073838;354.945179759665;374.097083060912;406.537596977238;405.666708376213;408.857840151081;419.263843225612;414.140586659037;452.195068504333;489.526564397767;487.898749442568;460.679939823402;761.201652069652;521.861314436800;490.300092715841]
%%
x = t';
x2 = [1:15]';
y = t1';
% This seems to be working
figure();
modelFun = @(b,x) b(1).*exp(b(2).*(x))+b(3);  
start = [-1000; -1; 100];

nlm = fitnlm(x,y,modelFun,start);
xx = linspace(0,40)';
plot(x,y,'o'); hold on;
line(xx,predict(nlm,xx),'linestyle','--','color','k')
hold off;
grid on;
xlabel('Time (min)');
ylabel('T1 (ms)')
text(8, 400, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f, R^2 = %.4f', [table2array(nlm.Coefficients(:,1)); nlm.Rsquared.Adjusted]), 'FontSize', 18);

hold on;
plot(x2, t1_dict_fit, 'o');
t1_dict_fit_med = medfilt1(t1_dict_fit,3);
p1 = plot(x2, t1_dict_fit_med, 'LineWidth', 1.5);

plot(x2, t1_molli_fit, 'o');
t1_molli_fit_med = medfilt1(t1_molli_fit,3);
p2 = plot(x2, t1_molli_fit_med, 'LineWidth', 1.5);
legend([p1, p2], 'Dict Fit', 'MOLLI Fit');
set(gca,'FontSize',18)

% w = [1 1 1 1 1 1 1 1 1 1 1 1]';
% wnlm = fitnlm(x,y,modelFun,start,'Weight',w);
% 
% figure();
% plot(x,y,'o'); hold on;
% line(xx,predict(wnlm,xx),'color','b');
% hold off;
% grid on;
% xlabel('Time (min)')
% ylabel('T1 (ms)')
% text(20, 400, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', table2array(nlm.Coefficients(:,1))));

%% T1 map values
figure('Position', [100 100 900 900]);
% imagesc(t1_map_3d_nt(:,:,3,1));
mapoi = t1_map_3d_nt(:,:,1,15);
imagesc(mapoi); axis image;
roi_coord = drawpolygon;
roi = createMask(roi_coord);

mean(nonzeros(roi .* mapoi))

%% For T2star
t = [1,4,6,8,10,12,14,16,25,27,30,31]; % min

x2 = [1:15]';
t2star_mGRE = [20.15, 24.88, 29.15, 28.37, nan, 30.37, 31.05, 30.34, 30.39, 32.32, 30.54, 28.14];
t2star_fit = [34.6629951076170;37.9634637671116;26.6064726314687;27.4264529502619;33.9960956064881;29.0026013945881;26.3592967474504;39.2714334404468;27.8109645369952;29.9374908786416;31.3614863759062;25.5883599133104;35.3968000891043;44.4150531573087;27.1834348873005];

x = t';
y = t2star_mGRE;
% This seems to be working
figure();
modelFun = @(b,x) b(1).*exp(b(2).*(x))+b(3);  
start = [-1000; -1; 100];

nlm = fitnlm(x,y,modelFun,start);
xx = linspace(0,40)';
plot(x,y,'o'); hold on;
line(xx,predict(nlm,xx),'linestyle','--','color','k')
hold off;
grid on;
xlabel('Time (min)');
ylabel('T1 (ms)')
text(8, 20, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f, R^2 = %.4f', [table2array(nlm.Coefficients(:,1)); nlm.Rsquared.Adjusted]), 'FontSize', 18);
ylim([18, 40]);

hold on;
plot(x2, t2star_fit, 'o');
t2star_fit_med = medfilt1(t2star_fit,3);
p1 = plot(x2, t2star_fit_med, 'LineWidth', 1.5);

%% 02/15/2022
% Discussed with Randy and doing T1 value progression on 
% Infarct, MVO, remote and blood pool
% Take 21P35 as an example

[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file));
figure();
for i = 1:15
subplot(4,4,i); imagesc(map_to_save.t1_map(:,:,i)); axis image; caxis([200 1000])
end

%% Figure and Draw ROIs
MI = zeros(size(map_to_save.t1_map,1), size(map_to_save.t1_map,2));
remote = zeros(size(map_to_save.t1_map,1), size(map_to_save.t1_map,2));
MVO = zeros(size(map_to_save.t1_map,1), size(map_to_save.t1_map,2));
bloodpool = zeros(size(map_to_save.t1_map,1), size(map_to_save.t1_map,2));
%%    
figure();
imagesc(abs(map_to_save.t1_map(:,:,end))); axis image; caxis([200 1000]);

roi = drawpolygon;
MI = createMask(roi);
% roi = drawpolygon;
% remote = createMask(roi);
% roi = drawpolygon;
% MVO = createMask(roi);
% roi = drawpolygon;
% bloodpool = createMask(roi);

%%
t1_map = map_to_save.t1_map;
MI_t1 = zeros(15,1);
remote_t1 = zeros(15,1);
MVO_t1 = zeros(15,1);
blood_t1 = zeros(15,1);

for i = 1:15
    MI_t1(i) = nanmean(nonzeros(MI .* t1_map(:,:,i)));
    remote_t1(i) = nanmean(nonzeros(remote .* t1_map(:,:,i)));
    MVO_t1(i) = nanmean(nonzeros(MVO .* t1_map(:,:,i)));
    blood_t1(i) = nanmean(nonzeros(bloodpool .* t1_map(:,:,i)));
end

figure();
x = (1:15)';
plot(x, MI_t1); hold on;
plot(x, remote_t1); 
plot(x, MVO_t1); 
plot(x, blood_t1);

legend({'MI', 'Remote', 'MVO', 'Blood Pool'});