function Func_plot_chord_analysis_EpiEndo(MI_Chord_Analysis, tp_dir2)
Mipix_mean_epi = MI_Chord_Analysis(end).Mipix_mean_epi;
Mipix_mean2_epi = MI_Chord_Analysis(end).Mipix_mean2_epi;
Mipix_mean3_epi = MI_Chord_Analysis(end).Mipix_mean3_epi;
Mipix_mean_endo = MI_Chord_Analysis(end).Mipix_mean_endo;
Mipix_mean2_endo = MI_Chord_Analysis(end).Mipix_mean2_endo;
Mipix_mean3_endo = MI_Chord_Analysis(end).Mipix_mean3_endo;

[I,J] = find(Mipix_mean_epi);
Mipix_mean_nnz_epi = [];
Mipix_mean2_nnz_epi = [];
Mipix_mean3_nnz_epi = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz_epi(idx) = Mipix_mean_epi(I(i), J(i));
    Mipix_mean2_nnz_epi(idx) = Mipix_mean2_epi(I(i), J(i));
    Mipix_mean3_nnz_epi(idx) = Mipix_mean3_epi(I(i), J(i));
    idx = idx + 1;
end

% remove NAs
nan_idx_epi = find(~isnan(Mipix_mean2_nnz_epi));
Mipix_mean_nnz_epi = Mipix_mean_nnz_epi(nan_idx_epi);
Mipix_mean2_nnz_epi = Mipix_mean2_nnz_epi(nan_idx_epi);
Mipix_mean3_nnz_epi = Mipix_mean3_nnz_epi(nan_idx_epi);

[I,J] = find(Mipix_mean_endo);
Mipix_mean_nnz_endo = [];
Mipix_mean2_nnz_endo = [];
Mipix_mean3_nnz_endo = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz_endo(idx) = Mipix_mean_endo(I(i), J(i));
    Mipix_mean2_nnz_endo(idx) = Mipix_mean2_endo(I(i), J(i));
    Mipix_mean3_nnz_endo(idx) = Mipix_mean3_endo(I(i), J(i));
    idx = idx + 1;
end

% remove NAs
nan_idx_endo = find(~isnan(Mipix_mean2_nnz_endo));
Mipix_mean_nnz_endo = Mipix_mean_nnz_endo(nan_idx_endo);
Mipix_mean2_nnz_endo = Mipix_mean2_nnz_endo(nan_idx_endo);
Mipix_mean3_nnz_endo = Mipix_mean3_nnz_endo(nan_idx_endo);


figure();
subplot(3,2,1);
imagesc(Mipix_mean_nnz_epi);
caxis([1000 1600]);
colorbar; title('T1 (ms)')
subplot(3,2,3);
imagesc(Mipix_mean2_nnz_epi);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,2,5);
imagesc(Mipix_mean3_nnz_epi);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));

subplot(3,2,2);
imagesc(Mipix_mean_nnz_endo);
caxis([1000 1600]);
colorbar; title('T1 (ms)')
subplot(3,2,4);
imagesc(Mipix_mean2_nnz_endo);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,2,6);
imagesc(Mipix_mean3_nnz_endo);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));
saveas(gcf, cat(2, tp_dir2, 'Chord_Display_EpiEndo.png'));


% Linear regression
x1_lb = 0;
x1_ub = 50;
x2_lb = 0;
x2_ub = 120;
y_lb = 800;
y_ub = 2000;

[x1_endo, b1_endo, X1_endo, yCalc1_endo, yCalc1_for_plot_endo, Rsq1_endo, str1_endo] = ...
    Func_LinearRegression(Mipix_mean2_nnz_endo, Mipix_mean_nnz_endo, x1_lb, x1_ub);
% R2star
[x2_endo, b2_endo, X2_endo, yCalc2_endo, yCalc2_for_plot_endo, Rsq2_endo, str2_endo] = ...
    Func_LinearRegression(Mipix_mean3_nnz_endo, Mipix_mean_nnz_endo, x2_lb, x2_ub);

[x1_epi, b1_epi, X1_epi, yCalc1_epi, yCalc1_for_plot_epi, Rsq1_epi, str1_epi] = ...
    Func_LinearRegression(Mipix_mean2_nnz_epi, Mipix_mean_nnz_epi, x1_lb, x1_ub);
% R2star
[x2_epi, b2_epi, X2_epi, yCalc2_epi, yCalc2_for_plot_epi, Rsq2_epi, str2_epi] = ...
    Func_LinearRegression(Mipix_mean3_nnz_epi, Mipix_mean_nnz_epi, x2_lb, x2_ub);

figure();
subplot(2,2,1);
scatter(Mipix_mean2_nnz_epi, Mipix_mean_nnz_epi, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([y_lb y_ub]); xlim([x1_lb x1_ub]);
title('Epi');
hold on;
plot(x1_epi,yCalc1_for_plot_epi, 'LineWidth', 1.5);
text(x1_ub-48, y_lb+200, str1_epi, 'FontSize', 12);
grid on;

subplot(2,2,3);
scatter(Mipix_mean3_nnz_epi, Mipix_mean_nnz_epi, 48, 'filled');
xlabel('R2star (Hz)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20)
ylim([y_lb y_ub]); xlim([x2_lb x2_ub]);
hold on;
plot(x2_epi,yCalc2_for_plot_epi, 'LineWidth', 1.5);
text(x2_ub-100, y_lb+200, str2_epi, 'FontSize', 12);
grid on;

subplot(2,2,2);
scatter(Mipix_mean2_nnz_endo, Mipix_mean_nnz_endo, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([y_lb y_ub]);xlim([x1_lb x1_ub]);
title('Endo');
hold on;
plot(x1_endo,yCalc1_for_plot_endo, 'LineWidth', 1.5);
text(x1_ub-48, y_lb+200, str1_endo, 'FontSize', 12);
grid on;

subplot(2,2,4);
scatter(Mipix_mean3_nnz_endo, Mipix_mean_nnz_endo, 48, 'filled');
xlabel('R2star (Hz)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([y_lb y_ub]); xlim([x2_lb x2_ub]);
hold on;
plot(x2_endo,yCalc2_for_plot_endo, 'LineWidth', 1.5);
text(x2_ub-100, y_lb+200, str2_endo, 'FontSize', 12);
grid on;

saveas(gcf, cat(2, tp_dir2, 'Chord_Scatter_EpiEndo_LReg.png'));

end