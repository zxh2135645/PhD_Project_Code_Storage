function Func_T1FP_ROI_Analysis(t1, ff, r2star, myo_t1, myo_ff, roi_in_myo_t1, roi_in_myo_ff, roi_in_myo_r2star, remote_in_myo_t1, remote_in_myo_ff, remote_in_myo_r2star,tp_dir2,name,time_point,strat_fname,status)

for i = 1:size(roi_in_myo_t1, 3)
    
    t1_roi_masked = roi_in_myo_t1(:,:,i) .* t1(:,:,i);
    ff_roi_masked = roi_in_myo_ff(:,:,i) .* ff(:,:,i);
    r2star_roi_masked = roi_in_myo_r2star(:,:,i) .* r2star(:,:,i);
    
    t1_remote_masked = remote_in_myo_t1(:,:,i) .* t1(:,:,i);
    ff_remote_masked = remote_in_myo_ff(:,:,i) .* ff(:,:,i);
    r2star_remote_masked = remote_in_myo_r2star(:,:,i) .* r2star(:,:,i);
    
end

end