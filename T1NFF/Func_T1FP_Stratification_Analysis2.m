function Func_T1FP_Stratification_Analysis2(t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,strat_fname,status)
strat = 7;
strat_cell = cell(strat, size(roi_in_myo_t1, 3));

for i = 1:size(roi_in_myo_t1, 3)
    if status(i) == 1
        img = t1(:,:,i);
        img2 = ff(:,:,i);
        img3 = r2star(:,:,i);
        
        moving = myo_t1(:,:,i);
        fixed = myo_ff(:,:,i);
        
        img2(img2 > 100) = 100;
        img2(img2 < 0) = 0;
        
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
        
        movingRegistered_myo_t1 = movepixels(myo_t1_temp,Bx,By)>0.5;
        movingRegistered_roi_t1 = movepixels(roi_in_myo_t1_temp,Bx,By)>0.5;
        movingRegistered_remote_t1 = movepixels(remote_in_myo_t1_temp,Bx,By)>0.5;
        
        univ_roi = roi_in_myo_ff(:,:,i) & movingRegistered_roi_t1;
        univ_myo = movingRegistered_myo_t1&fixed;
        
        
        %saveas(gcf, cat(2, tp_dir2, 'MyocardiumRegistration_demon_Slice', num2str(i), '.png'));
        fixed_eroded = imerode(myo_ff(:,:,i), se);
        BW_skel = bwmorph(fixed_eroded, 'skel', Inf);
%         center_fixed = imfill(BW_skel, 'hole');
%         center_fixed = imopen(center_fixed, se); % Removing spikes
%         fixedRegistered_epi = fixed_eroded - center_fixed > 0;
%         fixedRegistered_endo = center_fixed + fixed_eroded > 1;
        
        movingRegistered_eroded = imerode(movingRegistered_myo_t1, se);
        BW_skel = bwmorph(movingRegistered_eroded, 'skel', Inf);
%         center_moving = imfill(BW_skel, 'hole');
%         center_moving = imopen(center_moving, se); % Removing spikes
%         movingRegistered_epi = movingRegistered_eroded - center_moving > 0;
%         movingRegistered_endo = center_moving + movingRegistered_eroded > 1;
        
        univ_myo_eroded = imerode(univ_myo, se);
        
        h = figure('Position', [100 0 400 600]);
        subplot(3,2,1);
        temp = univ_myo_eroded .* img;
        temp(temp == 0) = nan;
        imagesc(temp); axis image; title('T1 (moving)'); colormap(brewermap([],'*RdYlBu'));axis off;
        caxis([700 1500]);
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        
        subplot(3,2,2);
        temp = univ_myo_eroded .* img2;
        temp(temp == 0) = nan;
        imagesc(temp); axis image; title('FF (fixed)'); colormap(brewermap([],'*RdYlBu'));axis off;
        caxis([0 50]);
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        
        ff_cat = zeros(size(img,1),size(img,2),strat);
        for n = 1:size(ff_cat, 3)
            if n ~= size(ff_cat, 3)
                ff_cat(:,:,n) = (fixed_eroded .* img2)>(n-1)*5 & (fixed_eroded .* img2)<=n*5;
            else
                ff_cat(:,:,n) = (fixed_eroded .* img2)>(n-1)*5;
            end
        end
        
        t1_cat = movingRegistered_eroded .* img .* ff_cat .* univ_roi;
        t1_array_mean = zeros(size(ff_cat, 3), 1);
        t1_array_sd = zeros(size(ff_cat, 3), 1);
        for n = 1:size(ff_cat, 3)
            t1_array_mean(n) = mean(nonzeros(t1_cat(:,:,n)));
            t1_array_sd(n) = std(nonzeros(t1_cat(:,:,n)));
            
            strat_cell{n,i} = nonzeros(t1_cat(:,:,n));
        end
        
        axis tight manual % this ensures that getframe() returns a consistent size
        %filename = 'testAnimated.gif';
        x = (0:(size(ff_cat,3)-1))*5+2.5;
        subplot(2,2,3); errorbar(x, t1_array_mean, t1_array_sd); xlim([0 size(ff_cat,3)*5]);
        xlabel('FF (%)'); ylabel('T1 (ms)'); title('FF vs T1');
        hold on;
        
        filename = cat(2, tp_dir2, 'FF_Stratify_Slice', num2str(i), '.gif');
        
        for n = 1:size(ff_cat, 3)
            subplot(2,2,3);
            h1 = xline((n-1)*5+2.5);
            drawnow 
            subplot(2,2,4); imagesc(t1_cat(:,:,n)); axis image; axis off; caxis([700 1500]);
            title(cat(2, 'Slice = ', num2str(i)));
            pause(1);
            
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            % Write to the GIF File
            if n == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
            set(h1, 'LineStyle', 'none');
        end

    end
    
    save(strat_fname, 'strat_cell');
end