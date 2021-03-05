function [snr] = SNR(images_input)
    imshow(images_input);
    disp('draw SI');
    roi = impoly;
    roi_mask = createMask(roi);
    maskSI = images_input.* roi_mask;
    maskSI = maskSI(:);
    SI_indices = maskSI~=0;
    SI = mean(maskSI(SI_indices));
    disp('SI');disp(SI);
    SI_sd = std(maskSI(SI_indices));
    disp('SI_sd');disp(SI_sd);

    disp('draw noise');
    roi = impoly;
    roi_mask = createMask(roi);
    maskNoise = images_input.*roi_mask;
    maskNoise = maskNoise(:);
    Noise_indices = maskNoise~=0;
    Noise = mean(maskNoise(Noise_indices));
    disp('Noise');disp(Noise);
    Noise_sd = std(maskNoise(Noise_indices));
    disp('Noise_sd');disp(Noise_sd);

    snr = SI/Noise_sd;
    disp('SNR');disp(snr);
end


