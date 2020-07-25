function [ax] = plotQCimages(infarct_mask, original_image, myocardium_mask, label)
% Plot a 2 x 2 plot of 2D image
mid_img_value = median(median(nonzeros(original_image)));

    switch label
        
        case 1
            ax = subplot(2,2,1);
            img_to_plot = original_image + mid_img_value * infarct_mask;
            imagesc(img_to_plot);
            axis image;
        case 2
            ax = subplot(2,2,2);
            img_to_plot = original_image + mid_img_value * infarct_mask;
            imagesc(img_to_plot);
            axis image;
        case 3
            ax = subplot(2,2,3);
            img_to_plot = original_image;
            imagesc(img_to_plot);
            axis image;
        case 4
            ax = subplot(2,2,4);
            img_to_plot = original_image .* myocardium_mask;
            imagesc(img_to_plot);
            axis image;
    end
end