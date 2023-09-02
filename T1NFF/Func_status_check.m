function [roi_in_myo_new, remote_in_myo_new, roi_new, remote_new, orig_new, myo_new] = ...
    Func_status_check(status_check, n, name, tp_count, roi_in_myo, remote_in_myo,...
    roi, remote, orig, myo)
slc_count = 1;
roi_in_myo_new = [];
remote_in_myo_new = [];
roi_new = [];
remote_new = [];
orig_new = [];
myo_new = [];

switch name
    
    case {'Ryn', 'Gobi', '18D16'}
        
        for i = 1:sum(status_check(n).status(tp_count,:))
            if  status_check(n).status(tp_count, i) == 1
                roi_in_myo_new(:,:,slc_count) = roi_in_myo(:,:,i);
                remote_in_myo_new(:,:,slc_count) = remote_in_myo(:,:,i);
                roi_new(:,:,slc_count) = roi(:,:,i);
                remote_new(:,:,slc_count) = remote(:,:,i);
                orig_new(:,:,slc_count) = orig(:,:,i);
                myo_new(:,:,slc_count) = myo(:,:,i);
                slc_count = slc_count + 1;
            end
        end
        
    otherwise
        for i = 1:length(status_check(n).status_final)
            if  status_check(n).status_final(i) == 1
                roi_in_myo_new(:,:,slc_count) = roi_in_myo(:,:,i);
                remote_in_myo_new(:,:,slc_count) = remote_in_myo(:,:,i);
                roi_new(:,:,slc_count) = roi(:,:,i);
                remote_new(:,:,slc_count) = remote(:,:,i);
                orig_new(:,:,slc_count) = orig(:,:,i);
                myo_new(:,:,slc_count) = myo(:,:,i);
                slc_count = slc_count + 1;
            end
        end
        
end

end