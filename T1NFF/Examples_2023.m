clear all;
close all;

%% FF
figure();
roi_in_myo_ff_nan = roi_in_myo_ff;
roi_in_myo_ff_nan(roi_in_myo_ff == 0) = nan;
myo_ff_nan = myo_ff;
myo_ff_nan(myo_ff == 0) = nan;

for slc = 1:size(ff,3)
    ff_slc = ff(:,:,slc);


    ax1 = subplot(2,3,slc);
    imagesc(ff_slc); pbaspect([size(ff_slc, 2), size(ff_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);


    ax2 = axes;
    imagesc(ax2, roi_in_myo_ff_nan(:,:,slc) .* ff_slc, 'AlphaData', myo_ff(:,:,slc));
    pbaspect([size(ff_slc, 2), size(ff_slc, 1) 1]); colormap(ax2, 'jet');
    caxis(ax1, [0 100]); caxis(ax2, [-2 50]); linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';

end

%% T1
figure();
%roi_in_myo_t1_nan = roi_in_myo_t1;
%roi_in_myo_t1_nan = double(roi_rimmed_t1_new);
%roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;
%roi_in_myo_t1_nan(roi_rimmed_t1_new == 0) = nan;

roi_in_myo_t1_nan = double(roi_in_myo_t1);
roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;

myo_t1_nan = myo_t1;
myo_t1_nan(myo_t1 == 0) = nan;
for slc = 1:size(roi_in_myo_t1_nan,3)
    t1_slc = t1(:,:,slc);

    [x,y] = find(remote_in_myo_t1(:,:,slc) == 1);

    t1_array = zeros(length(x),1);
    for i = 1:length(x)
        t1_array(i) = t1(x(i),y(i),slc);
    end

    t1_array(t1_array < 0) = 0;
    mean_t1_remote = mean(t1_array);

    ax1 = subplot(2,3,slc);
    imagesc(t1_slc); pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);

    ax2 = axes;
    imagesc(ax2, myo_t1(:,:,slc), 'AlphaData', myo_t1(:,:,slc));
    pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax2, 'cool');
    ax2.Visible = 'off';

    ax3 = axes;
    imagesc(ax3, roi_in_myo_t1_nan(:,:,slc) .* t1_slc - mean_t1_remote, 'AlphaData', roi_in_myo_t1(:,:,slc));
    pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax3, 'cool');
    caxis(ax1, [0 2000]); caxis(ax2, [0 2]); caxis(ax3, [-500 500]); linkprop([ax1 ax2 ax3], 'Position');
    ax3.Visible = 'off';

end

%% R2star
figure();
roi_in_myo_r2star_nan = roi_in_myo_ff;
roi_in_myo_ff_nan(roi_in_myo_ff == 0) = nan;
myo_ff_nan = myo_ff;
myo_ff_nan(myo_ff == 0) = nan;
for slc = 1:size(r2star,3)
    r2star_slc = r2star(:,:,slc);

    ax1 = subplot(2,2,slc);
    imagesc(r2star_slc); pbaspect([size(r2star_slc, 2), size(r2star_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);


    ax2 = axes;
    imagesc(ax2, roi_in_myo_ff_nan(:,:,slc) .* r2star_slc, 'AlphaData', myo_ff(:,:,slc));
    pbaspect([size(r2star_slc, 2), size(r2star_slc, 1) 1]); colormap(ax2, 'jet');
    caxis(ax1, [0 100]); caxis(ax2, [0 100]); linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';
end

%%  Raw
slc = 2;
figure(); 
imagesc(t1(:,:,slc)); caxis([100 1800]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off; 

%% T1 example
figure();
%roi_in_myo_t1_nan = roi_in_myo_t1;
%roi_in_myo_t1_nan = double(roi_rimmed_t1_new);
%roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;
%roi_in_myo_t1_nan(roi_rimmed_t1_new == 0) = nan;

roi_in_myo_t1_nan = double(roi_in_myo_t1);
roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;

myo_t1_nan = myo_t1;
myo_t1_nan(myo_t1 == 0) = nan;

t1_slc = t1(:,:,slc);

[x,y] = find(roi_in_myo_t1(:,:,slc) == 1);

t1_array = zeros(length(x),1);
for i = 1:length(x)
    t1_array(i) = t1(x(i),y(i),slc);
end

t1_array(t1_array < 0) = 0;
mean_t1_roi = mean(t1_array);

ax1 = axes;
imagesc(t1_slc); pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]);
colormap(ax1, 'gray');
set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);

ax2 = axes;
imagesc(ax2, myo_t1(:,:,slc), 'AlphaData', myo_t1(:,:,slc));
pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax2, 'cool');
ax2.Visible = 'off';

ax3 = axes;
imagesc(ax3, roi_in_myo_t1_nan(:,:,slc) .* t1_slc - mean_t1_roi, 'AlphaData', roi_in_myo_t1(:,:,slc));
pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax3, 'hot');
caxis(ax1, [100 1800]); caxis(ax2, [0 2]); caxis(ax3, [-500 500]); linkprop([ax1 ax2 ax3], 'Position');
ax3.Visible = 'off';
colorbar
%% T2
figure();
roi_in_myo_t2_nan = roi_in_myo_t2;
roi_in_myo_t2_nan(roi_in_myo_t2 == 0) = nan;
myo_t2_nan = myo_t2;
myo_t2_nan(myo_t2 == 0) = nan;
for slc = 1:size(t2,3)
    t2_slc = t2(:,:,slc);

    ax1 = subplot(2,2,slc);
    imagesc(t2_slc); pbaspect([size(t2_slc, 2), size(t2_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);


    ax2 = axes;
    imagesc(ax2, roi_in_myo_t2_nan(:,:,slc) .* t2_slc, 'AlphaData', myo_t2(:,:,slc));
    pbaspect([size(t2_slc, 2), size(t2_slc, 1) 1]); colormap(ax2, 'jet');
    caxis(ax1, [0 100]); caxis(ax2, [10 80]); linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';
end

%% Analysis
slc = 3;
[x,y] = find(roi_in_myo_ff(:,:,slc) == 1);

ff_array = zeros(length(x),1);
for i = 1:length(x)
    ff_array(i) = ff(x(i),y(i),slc);
end

ff_array(ff_array < 0) = 0;

[x,y] = find(roi_in_myo_t1(:,:,slc) == 1);

t1_array = zeros(length(x),1);
for i = 1:length(x)
    t1_array(i) = t1(x(i),y(i),slc);
end

t1_array(t1_array < 0) = 0;
%%
mean(t1_array)
std(t1_array)
mean(ff_array)
std(ff_array)
figure(); histogram(t1_array, 'Normalization','probability'); xlim([800 1800]);
figure(); histogram(ff_array, 'Normalization','probability'); xlim([0 50]);

%% Plot Linear regression
hetero_t1 = [93.4948784	131.5470668	91.00163032	76.80994275	56.58641649	73.68869728	50.58499628	66.03825125	46.06961325	19.16410659	48.13226983	30.41213435	27.94921931	30.87181825	34.33708437	24.12555005	24.78094678	32.31641443	44.70988451	66.47061968	56.32420973	29.67034455	65.84100505	31.37680055	66.09882975	80.84104207	80.26738402	54.83245465	47.18933309	89.61288998	63.37088444	40.54647826	43.68070754	34.74420871	79.25508985	43.33666063	48.80817989	100.6724241	67.10869104	67.54054888	48.25661544	56.34311632	31.34474322	58.01363067	76.55004576	50.38493558	16.17589911	12.51768825	15.36003866];
ff_mean = [15.60447728	14.78351024	12.17616673	8.138109019	6.199221444	6.298566453	6.349978604	3.980585046	4.247593918	4.852949747	4.657240867	5.177199829	1.631045615	2.728173187	2.586949272	2.267199355	1.965339977	3.132784497	3.067722853	3.436115497	3.447189287	4.457039614	3.453607224	1.485206799	10.56805944	11.57662087	12.4761005	6.204056188	6.20985808	7.109428801	8.894293271	5.763136511	8.114717589	4.279872056	2.760749916	3.46983488	8.439058352	12.71780664	3.220006784	9.364701114	7.498762478	3.849234679	3.67828082	2.795359749	3.624650895	2.636651685	1.430597131	3.262081827	4.778726584];

clusters = kmeans([hetero_t1',ff_mean'], 2);


figure();
plot(hetero_t1(clusters == 1), ff_mean(clusters == 1), 'o');
hold on;
plot(hetero_t1(clusters == 2), ff_mean(clusters == 2), 'o');

%%
clusters = zeros(size(hetero_t1'));
clusters(hetero_t1 < 40 & ff_mean < 5.5) = 1;
clusters(hetero_t1 >= 40 & ff_mean < 5.5) = 2;
clusters(hetero_t1 >= 40 & ff_mean >= 5.5) = 3;
figure();
plot(hetero_t1(clusters == 1), ff_mean(clusters == 1), 'ok', 'MarkerSize', 10);
hold on;
plot(hetero_t1(clusters == 2), ff_mean(clusters == 2), 'ok', 'MarkerSize', 10);
plot(hetero_t1(clusters == 3), ff_mean(clusters == 3), 'ok', 'MarkerSize', 10);
yline(5.5);
xline(40);

mdl1 = fitlm(hetero_t1(clusters == 1), ff_mean(clusters == 1))
mdl2 = fitlm(hetero_t1(clusters == 2), ff_mean(clusters == 2))
mdl3 = fitlm(hetero_t1(clusters == 3), ff_mean(clusters == 3))

mdl13 = fitlm(hetero_t1(clusters == 1 | clusters == 3), ff_mean(clusters == 1 | clusters == 3))
%figure(); plot(mdl13)

ff = [0, 5.5 16];
hetero = [0 40 140];
patch([hetero(1) hetero(2) hetero(2) hetero(1)], [ff(1) ff(1) ff(2) ff(2)], [241 194 151]/255, 'FaceAlpha',.5)
patch([hetero(2) hetero(3) hetero(3) hetero(2)], [ff(1) ff(1) ff(2) ff(2)], [199 213 161]/255, 'FaceAlpha',.5)
patch([hetero(2) hetero(3) hetero(3) hetero(2)], [ff(2) ff(2) ff(3) ff(3)], [159 203 219]/255, 'FaceAlpha',.5)
patch([hetero(1) hetero(2) hetero(2) hetero(1)], [ff(2) ff(2) ff(3) ff(3)], [98 141 207]/255, 'FaceAlpha',.5)


%%
%% LGE
figure();
%roi_in_myo_t1_nan = roi_in_myo_t1;
roi_in_myo_lge_nan = double(roi_in_myo_lge);
%roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;
roi_in_myo_lge_nan(roi_in_myo_lge == 0) = nan;

myo_lge_nan = myo_lge;
myo_lge_nan(myo_lge == 0) = nan;
for slc = 1:size(roi_in_myo_lge_nan,3)
    lge_slc = lge(:,:,slc);

    % [x,y] = find(roi_in_myo_lge(:,:,slc) == 1);
    [x,y] = find(roi_in_myo_lge(:,:,slc) == 1);

    lge_array = zeros(length(x),1);
    for i = 1:length(x)
        lge_array(i) = lge(x(i),y(i),slc);
    end

    lge_array(lge_array < 0) = 0;
    mean_lge_remote = mean(lge_array);% + 5*std(lge_array);

    ax1 = subplot(2,2,slc);
    imagesc(lge_slc); pbaspect([size(lge_slc, 2), size(lge_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);

    ax2 = axes;
    imagesc(ax2, myo_lge(:,:,slc), 'AlphaData', myo_lge(:,:,slc));
    pbaspect([size(lge_slc, 2), size(lge_slc, 1) 1]); colormap(ax2, 'cool');
    ax2.Visible = 'off';

    ax3 = axes;
    imagesc(ax3, roi_in_myo_lge_nan(:,:,slc) .* lge_slc - mean_lge_remote, 'AlphaData', roi_in_myo_lge(:,:,slc));
    pbaspect([size(lge_slc, 2), size(lge_slc, 1) 1]); colormap(ax3, 'cool');
    caxis(ax1, [0 4196]); caxis(ax2, [0 2]); caxis(ax3, [-500 500]); linkprop([ax1 ax2 ax3], 'Position');
    ax3.Visible = 'off';

end



