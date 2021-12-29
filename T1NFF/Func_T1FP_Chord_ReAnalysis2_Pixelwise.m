function Func_T1FP_Chord_ReAnalysis2_Pixelwise(Segn, Groove, t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,LR_mdl_fname,chord_values_fname)

t1_cell = cell(1,size(roi_in_myo_t1, 3));
ff_cell = cell(1,size(roi_in_myo_t1, 3));
r2star_cell = cell(1,size(roi_in_myo_t1, 3));
t1_remote_cell = cell(1,size(roi_in_myo_t1, 3));
ff_remote_cell = cell(1,size(roi_in_myo_t1, 3));
r2star_remote_cell = cell(1,size(roi_in_myo_t1, 3));
t1_hemo_cell = cell(1,size(roi_in_myo_t1, 3));
ff_hemo_cell = cell(1,size(roi_in_myo_t1, 3));
r2star_hemo_cell = cell(1,size(roi_in_myo_t1, 3));

for i = 1:size(roi_in_myo_t1, 3)
    img = t1(:,:,i);
    img2 = ff(:,:,i);
    img3 = r2star(:,:,i);
    
    fixed = myo_ff(:,:,i);
    moving = myo_t1(:,:,i);
    
    %img2(img2 > 100) = 100;
    %img2(img2 < 0) = 0;
    
    figure('Position', [100 0 400 800]);
    subplot(3,2,1);
    imagesc(fixed); axis image; title('FF (fixed)'); colormap(brewermap([],'*RdYlBu'));axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,2);
    imagesc(moving); axis image; title('T1 (moving)');axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    se = strel('disk', 1);
    
    I1 = moving; I2 = fixed;
    % Set static and moving image
    S=I2; M=I1;
    
    % resizepercentag
    [movingRegistered,Bx,By,Fx,Fy] = register_images(M,S);
    movingRegistered = movingRegistered > 0.5;
    
    if (size(img, 1) ~= size(img2, 1)) || (size(img, 2) ~= size(img2, 2))
        img = imresize(img,size(img2),'bicubic');
        myo_t1_temp = imresize(myo_t1(:,:,i),size(img2),'bicubic');
        roi_in_myo_t1_temp = imresize(roi_in_myo_t1(:,:,i),size(img2),'bicubic');
        remote_in_myo_t1_temp = imresize(remote_in_myo_t1(:,:,i),size(img2),'bicubic');
    else
        myo_t1_temp = myo_t1(:,:,i);
        roi_in_myo_t1_temp = roi_in_myo_t1(:,:,i);
        remote_in_myo_t1_temp = remote_in_myo_t1(:,:,i);
    end
    
    img = movepixels(img,Bx,By);
    
    subplot(3,2,3); imagesc(img2); axis image; title('FF map'); axis off; caxis([0 50]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,4); imagesc(img); axis image; title('T1 map'); axis off;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    movingRegistered_myo_t1 = movepixels(myo_t1_temp,Bx,By)>0.5;
    movingRegistered_roi_t1 = movepixels(roi_in_myo_t1_temp,Bx,By)>0.5;
    movingRegistered_remote_t1 = movepixels(remote_in_myo_t1_temp,Bx,By)>0.5;
    
    univ_roi = roi_in_myo_ff(:,:,i) & movingRegistered_roi_t1;
    univ_myo = movingRegistered_myo_t1&fixed;
    univ_myo_eroded = imerode(univ_myo, se);
    
    subplot(3,2,5);
    imshowpair(fixed,movingRegistered,'Scaling','joint'); title('Registered');
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    subplot(3,2,6); imagesc(double(univ_myo) + double(univ_roi) + 2*double(movingRegistered_remote_t1)); axis image;
    title(cat(2, 'Slice = ', num2str(i)));
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    saveas(gcf, cat(2, tp_dir2, 'MyocardiumRegistration_demon_Slice', num2str(i), '.png'));
    
    univ_myo = movingRegistered_myo_t1&fixed;
    
    fixed_eroded = imerode(myo_ff(:,:,i), se);
    BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
    center_fixed = imfill(BW_skel, 'hole');
    center_fixed = imopen(center_fixed, se); % Removing spikes
    fixedRegistered_epi = fixed_eroded - center_fixed > 0;
    fixedRegistered_endo = center_fixed + fixed_eroded > 1;
    
    movingRegistered_eroded = imerode(movingRegistered_myo_t1, se);
    BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
    center_moving = imfill(BW_skel, 'hole');
    center_moving = imopen(center_moving, se); % Removing spikes
    movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
    movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
    
    univ_myo_eroded = imerode(univ_myo, se);
    
    remote_mean_t1 = mean(nonzeros(movingRegistered_remote_t1 .* img));
    remote_sd_t1 = std(nonzeros(movingRegistered_remote_t1 .* img));
    thresh = remote_mean_t1 - 2*remote_sd_t1;
    hemo_mask = (img<thresh).*univ_roi;
    
    % For debugging
    img_t1_masked = img .* (univ_roi-hemo_mask) .* movingRegistered_eroded;
    img_ff_masked = img2 .* (univ_roi-hemo_mask) .* fixed_eroded;
    img_r2star_masked = img3 .* (univ_roi-hemo_mask) .* fixed_eroded;
    
    img_t1_remote_masked = img .* movingRegistered_remote_t1 .* movingRegistered_eroded;
    img_ff_remote_masked = img2 .* movingRegistered_remote_t1 .* fixed_eroded;
    img_r2star_remote_masked = img3 .* movingRegistered_remote_t1 .* fixed_eroded;
    
    img_t1_hemo_masked = img .* hemo_mask .* movingRegistered_eroded;
    img_ff_hemo_masked = img2 .* hemo_mask .* fixed_eroded;
    img_r2star_hemo_masked = img3 .* hemo_mask .* fixed_eroded;
    
    %
    if any(img_t1_masked(:))
        [row,col,img_t1_masked_1d] = find(img_t1_masked);
        sz = size(img_ff_masked);
        ind = sub2ind(sz,row,col);
        img_ff_masked_1d = img_ff_masked(ind);
        img_r2star_masked_1d = img_r2star_masked(ind);
        
        t1_cell{i} = img_t1_masked_1d;
        ff_cell{i} = img_ff_masked_1d;
        r2star_cell{i} = img_r2star_masked_1d;
    end
    
    if any(img_t1_remote_masked(:))
        % remote
        [row,col,img_t1_remote_masked_1d] = find(img_t1_remote_masked);
        sz = size(img_ff_remote_masked);
        ind = sub2ind(sz,row,col);
        img_ff_remote_masked_1d = img_ff_remote_masked(ind);
        img_r2star_remote_masked_1d = img_r2star_remote_masked(ind);
       
        t1_remote_cell{i} = img_t1_remote_masked_1d;
        ff_remote_cell{i} = img_ff_remote_masked_1d;
        r2star_remote_cell{i} = img_r2star_remote_masked_1d;
    end
    
    if any(img_t1_hemo_masked(:))
        [row,col,img_t1_hemo_masked_1d] = find(img_t1_hemo_masked);
        sz = size(img_ff_hemo_masked);
        ind = sub2ind(sz,row,col);
        img_ff_hemo_masked_1d = img_ff_hemo_masked(ind);
        img_r2star_hemo_masked_1d = img_r2star_hemo_masked(ind);
        
        t1_hemo_cell{i} = img_t1_hemo_masked_1d;
        ff_hemo_cell{i} = img_ff_hemo_masked_1d;
        r2star_hemo_cell{i} = img_r2star_hemo_masked_1d;
    end
end

%% Plot one linear regression
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mean_t1_array_nz = [];
mean_ff_array_nz = [];
mean_r2star_array_nz = [];

mean_t1_hemo_array_nz = [];
mean_ff_hemo_array_nz = [];
mean_r2star_hemo_array_nz = [];

mean_t1_array_nz_remote = [];
mean_ff_array_nz_remote = [];
mean_r2star_array_nz_remote = [];

for i = 1:length(t1_remote_cell)
    mean_t1_array_nz = [mean_t1_array_nz; t1_cell{i}];
    mean_ff_array_nz = [mean_ff_array_nz; ff_cell{i}];
    mean_r2star_array_nz = [mean_r2star_array_nz; r2star_cell{i}];    
    
    mean_t1_hemo_array_nz = [mean_t1_hemo_array_nz; t1_hemo_cell{i}];
    mean_ff_hemo_array_nz = [mean_ff_hemo_array_nz; ff_hemo_cell{i}];    
    mean_r2star_hemo_array_nz = [mean_r2star_hemo_array_nz; r2star_hemo_cell{i}];
    
    mean_t1_array_nz_remote = [mean_t1_array_nz_remote; t1_remote_cell{i}];
    mean_ff_array_nz_remote = [mean_ff_array_nz_remote; ff_remote_cell{i}];
    mean_r2star_array_nz_remote = [mean_r2star_array_nz_remote; r2star_remote_cell{i}];
end

mean_ff_array_nz_new = mean_ff_array_nz;
mean_t1_array_nz_new = mean_t1_array_nz;
mean_r2star_array_nz_new = mean_r2star_array_nz;

mean_t1_hemo_array_nz_new = mean_t1_hemo_array_nz;
mean_r2star_hemo_array_nz_new = mean_r2star_hemo_array_nz;
mean_ff_hemo_array_nz_new = mean_ff_hemo_array_nz;

mean_ff_array_nz_new(mean_ff_array_nz<0 | mean_ff_array_nz>100) = [];
mean_t1_array_nz_new(mean_ff_array_nz<0 | mean_ff_array_nz>100) = [];
mean_r2star_array_nz_new(mean_ff_array_nz<0 | mean_ff_array_nz>100) = [];

mean_ff_hemo_array_nz_new(mean_r2star_hemo_array_nz<0) = [];
mean_t1_hemo_array_nz_new(mean_r2star_hemo_array_nz<0) = [];
mean_r2star_hemo_array_nz_new(mean_r2star_hemo_array_nz<0) = [];

mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
mean_r2star_array_nz_new_remote = mean_r2star_array_nz_remote;

mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<0 | mean_ff_array_nz_remote>100) = [];
mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<0 | mean_ff_array_nz_remote>100) = [];
mean_r2star_array_nz_new_remote(mean_ff_array_nz_remote<0 | mean_ff_array_nz_remote>100) = [];

plot_save = cat(2, tp_dir2, 'LinearReg/');
if ~exist(plot_save, 'dir')
    mkdir(plot_save);
end

mdl = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);

if any(mean_r2star_hemo_array_nz_new)
    mdl2 = fitlm(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new);
    
    figure('Position', [100 100 600 800]);
    title(cat(2, name, ' ', time_point))
    scatter(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_r2star_hemo_array_nz_new .* mdl2.Coefficients.Estimate(2) + mdl2.Coefficients.Estimate(1);
    hold on;
    plot(mean_r2star_hemo_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_r2star_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    xlabel('R2star (s^{-1})');
    ylabel('T1 (ms)');
    yl = ylim;
    xl = xlim;
    text(0.5*xl(2), yl(1)+200, cat(2,'Y = ', num2str(mdl2.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl2.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl2.Rsquared.Ordinary,3)), 'FontSize', 12)
    saveas(gcf, cat(2, plot_save, 'r2starVST1_1line_demon_in_hemo_pixelwise.png'));
else
    mdl2 = struct;
    
    figure('Position', [100 100 600 800]);
    title(cat(2, name, ' ', time_point))
    scatter(mean_r2star_hemo_array_nz_new, mean_t1_hemo_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = [];
    hold on;
    plot(mean_r2star_hemo_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_r2star_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    xlabel('R2star (s^{-1})');
    ylabel('T1 (ms)');
    yl = ylim;
    xl = xlim;
    text(0.5*xl(2), 0.95*yl(2), 'NA', 'FontSize', 12)
    saveas(gcf, cat(2, plot_save, 'r2starVST1_1line_demon_in_hemo_pixelwise.png'));
end

chord_value_results_pixelwise = struct;
chord_value_results_pixelwise.mean_t1_array = mean_t1_array_nz_new;
chord_value_results_pixelwise.mean_ff_array = mean_ff_array_nz_new;
chord_value_results_pixelwise.mean_r2star_array = mean_r2star_array_nz_new;
chord_value_results_pixelwise.mean_t1_array_remote = mean_t1_array_nz_new_remote;
chord_value_results_pixelwise.mean_ff_array_remote = mean_ff_array_nz_new_remote;
chord_value_results_pixelwise.mean_r2star_array_remote = mean_r2star_array_nz_new_remote;
chord_value_results_pixelwise.mean_t1_hemo_array = mean_t1_hemo_array_nz;
chord_value_results_pixelwise.mean_ff_hemo_array = mean_ff_hemo_array_nz;
chord_value_results_pixelwise.mean_r2star_hemo_array = mean_r2star_hemo_array_nz;

save(chord_values_fname, '-struct', 'chord_value_results_pixelwise');

figure('Position', [100 100 600 800]);
title(cat(2, name, ' ', time_point))
rows = ceil(size(roi_in_myo_t1, 3) / 2) + 1;
subplot(rows,1,rows);
scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_nz_new .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
xlabel('FF (%)');
ylabel('T1 (ms)');
yl = ylim;
xl = xlim;
text(0.5*xl(2), yl(1)+200, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

mdl_general_ffvst1 = mdl;
mdl_general_r2vst1_in_hemo = mdl2;

mdl_slc_cell = cell(size(roi_in_myo_t1, 3), 1);
for i = 1:size(roi_in_myo_t1, 3)
    mean_ff_array_nz = ff_cell{i};
    mean_t1_array_nz = t1_cell{i};
    mean_ff_array_nz_new = mean_ff_array_nz;
    mean_t1_array_nz_new = mean_t1_array_nz;
    mean_ff_array_nz_new(mean_ff_array_nz<0 | mean_ff_array_nz>100) = [];
    mean_t1_array_nz_new(mean_ff_array_nz<0 | mean_ff_array_nz>100) = [];
    
    mean_ff_array_nz_remote = ff_remote_cell{i};
    mean_t1_array_nz_remote = t1_remote_cell{i};
    mean_ff_array_nz_new_remote = mean_ff_array_nz_remote;
    mean_t1_array_nz_new_remote = mean_t1_array_nz_remote;
    mean_ff_array_nz_new_remote(mean_ff_array_nz_remote<0 | mean_ff_array_nz_remote>100) = [];
    mean_t1_array_nz_new_remote(mean_ff_array_nz_remote<0 | mean_ff_array_nz_remote>100) = [];
    
    mdl_slc = fitlm(mean_ff_array_nz_new, mean_t1_array_nz_new);
    
    subplot(rows,2,i);
    scatter(mean_ff_array_nz_new, mean_t1_array_nz_new, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
    Y = mean_ff_array_nz_new .* mdl_slc.Coefficients.Estimate(2) + mdl_slc.Coefficients.Estimate(1);
    hold on;
    plot(mean_ff_array_nz_new, Y, 'k', 'LineWidth', 1);
    scatter(mean_ff_array_nz_new_remote, mean_t1_array_nz_new_remote, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
    text(0.2*xl(2), yl(1)+100, cat(2,'Y = ', num2str(mdl_slc.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_slc.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_slc.Rsquared.Ordinary,3)), 'FontSize', 10)
    
    xlim(xl); ylim(yl);
    title(cat(2, 'Slice ', num2str(i)));
    xlabel('FF (%)');
    ylabel('T1 (ms)');
    
    mdl_slc_cell{i} = mdl_slc;
end

saveas(gcf, cat(2, plot_save, 'T1vsFF_1line_demon_pixelwise.png'));

mdl_results_pixelwise = struct;
mdl_results_pixelwise.mdl_general_ffvst1 = mdl;
mdl_results_pixelwise.mdl_general_r2vst1_in_hemo = mdl2;
mdl_results_pixelwise.mdl_slc = mdl_slc_cell;

save(LR_mdl_fname, '-struct', 'mdl_results_pixelwise');



end