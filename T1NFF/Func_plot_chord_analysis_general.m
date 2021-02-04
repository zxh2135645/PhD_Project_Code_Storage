function Func_plot_chord_analysis_general(MI_Chord_Analysis, tp_dir2)
Mipix_mean = MI_Chord_Analysis(end).Mipix_mean;
Mipix_mean2 = MI_Chord_Analysis(end).Mipix_mean2;
Mipix_mean3 = MI_Chord_Analysis(end).Mipix_mean3;

[I,J] = find(Mipix_mean);
Mipix_mean_nnz = [];
Mipix_mean2_nnz = [];
Mipix_mean3_nnz = [];
idx = 1;
for i = 1:length(I)
    Mipix_mean_nnz(idx) = Mipix_mean(I(i), J(i));
    Mipix_mean2_nnz(idx) = Mipix_mean2(I(i), J(i));
    Mipix_mean3_nnz(idx) = Mipix_mean3(I(i), J(i));
    idx = idx + 1;
end

figure();
subplot(3,1,1);
imagesc(Mipix_mean_nnz);
colorbar; title('T1 (ms)')
caxis([1000 1600]);
subplot(3,1,2);
imagesc(Mipix_mean2_nnz);
colorbar; title('FF (%)')
caxis([0 20]);
subplot(3,1,3);
imagesc(Mipix_mean3_nnz);
colorbar; title('R2star (s^{-1})')
caxis([0 100]);
colormap(brewermap([],'*RdYlBu'));

saveas(gcf, cat(2, tp_dir2, 'Chord_Display.png'));


x1_lb = 0;
x1_ub = 50;
x2_lb = 0;
x2_ub = 120;
y_lb = 800;
y_ub = 2000;
% Linear regression T1
[x1, b1, X1, yCalc1, yCalc1_for_plot, Rsq1, str1] = ...
    Func_LinearRegression(Mipix_mean2_nnz, Mipix_mean_nnz, x1_lb, x1_ub);
% R2star
[x2, b2, X2, yCalc2, yCalc2_for_plot, Rsq2, str2] = ...
    Func_LinearRegression(Mipix_mean3_nnz, Mipix_mean_nnz, x2_lb, x2_ub);

% scatter
figure();
subplot(2,1,1);
scatter(Mipix_mean2_nnz, Mipix_mean_nnz, 48, 'filled');
xlabel('FF (%)'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([y_lb y_ub]); xlim([x1_lb x1_ub]);
hold on;
plot(x1,yCalc1_for_plot, 'LineWidth', 1.5);
text(x1_ub-30, y_lb+200, str1, 'FontSize', 16);
grid on;

subplot(2,1,2);
scatter(Mipix_mean3_nnz, Mipix_mean_nnz, 48, 'filled');
xlabel('R2star (s^{-1})'); ylabel('T1 (ms)');
set(gca, 'FontSize', 20);
ylim([y_lb y_ub]); xlim([x2_lb x2_ub]);
hold on;
plot(x2,yCalc2_for_plot, 'LineWidth', 1.5);
text(x2_ub-60, y_lb+200, str2, 'FontSize', 16);
grid on;

saveas(gcf, cat(2, tp_dir2, 'Chord_Scatter_LReg.png'));
end