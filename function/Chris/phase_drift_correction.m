function iField_drift_corrected = phase_drift_correction(iField, TE, voxel_size)
    phase_diff = angle(iField(:,:,:,2)./iField(:,:,:,1));
    iMag= sqrt(sum(abs(iField(:,:,:,1:end)).^2,4));
    iMag = iMag./max(iMag(:));
    Mask = autoMask(iMag,voxel_size);
    p = 2;
    [phase_diff_unwrap,~,~] = phase_unwrap_3d_UNIC_Chris_1002(phase_diff,p,iMag,voxel_size,Mask);
    phase_diff_unwrap_dmean = phase_diff_unwrap - round(mean(phase_diff_unwrap(Mask~=0)/2/pi))*2*pi;
    phase_diff_unwrap_dmean_te1 = phase_diff_unwrap_dmean./(TE(2) - TE(1))*TE(1);

    offsets = iField(:,:,:,1)./exp(1j*phase_diff_unwrap_dmean_te1);
    offsets = offsets./abs(offsets); % complex phase offset
    offsets(isnan(offsets)) = 0;

    % offsets(:,:,:,1) = smooth3(offsets(:,:,:,1),'box',round(5)*2+1); 
    % %       offsets(:,:,:,1,chan) = smooth3(offsets(:,:,:,1,chan),'box',round(2./vox)*2+1); 
    % offsets(:,:,:,1) = offsets(:,:,:,1)./abs(offsets(:,:,:,1));
    % 
    % offsets(:,:,:,1) = imgaussfilt3(real(offsets(:,:,:,1)),6) + 1j*imgaussfilt3(imag(offsets(:,:,:,1)),6);
    % offsets(:,:,:,1,chan) = offsets(:,:,:,1)./abs(offsets(:,:,:,1));

    iField_drift_corrected = iField./offsets;
end