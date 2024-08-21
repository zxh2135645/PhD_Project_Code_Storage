clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% pre_QualControl.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% data_storage_rim.mat
% Time_Evolution_rim.png
% Time_Evolution.png
% for_analysis.mat
% for_analysis_rim.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../function/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
%time_points = {'0D_baseline','1D', '7D', '28D', '8WK', '6MO', '9MO', '1YR', '15YR'};
time_points = {'7D', '8WK', '12WK', '14WK' '4MO', '6MO', '9MO', '1YR', '15YR'};

dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    };

sequence_label = {'T1', 'T2star', 'LGE'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));

%% Before analysis, parse pre_QualControl
status_check = struct;
for n = 1:length(Names)
    name = Names{n};
    status_check(n).Name = name;
    status_check(n).status = [];
    status_check(n).status_final = [];
    tp_count = 1;
    for tp = 1:length(time_points)

        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');

        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
                slc_loc = cat(2, 'Slice', num2str(i));
                if pre_QualControl(n).status(end-tp+1).(slc_loc) == 1
                    status_check(n).status(tp_count, i) = 1;
                elseif pre_QualControl(n).status(end-tp+1).(slc_loc) == 0
                    status_check(n).status(tp_count, i) = 0;
                end
            end

            tp_count = tp_count + 1;
        end
    end

    for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
        if i <= size(status_check(n).status,2)
            status_check(n).status_final(1,i) = all(status_check(n).status(:,i));
            % The timepoints of status_check goes from end to beginning
        end
    end
end

%% Get the matrics of T1, FF and R2star ( exclude edges of MI region)
%for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality
label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis_rim = struct;
data_storage_rim = struct;

time_points = {'6MO', '9MO', '1YR', '15YR'};

for n = 1:length(Names)
%for n = 12:12
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    for_analysis_rim(n).Name = name;
    for_analysis_rim(n).metrics = struct;

    tp_count = 1;
    data_storage_rim(n).Name = name;
    data_storage_rim(n).data = struct;

    metrics(n).name = name;
    metrics(n).TimePoints = struct;

    for tp = 1:length(time_points)
    %for tp = 4:4
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            myo_glob = glob(cat(2, tp_dir, label_t1, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{6}, '/*'));
            
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
            
            [slc_array_t1, idx_reordered] = sort(slc_array);
            
            roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
            roi_t1 = roi_in_myo_t1 .* vol_img_3D;
            remote_t1 = remote_in_myo_t1 .* vol_img_3D;
            t1 = vol_img_3D;
            myo_t1 = mask_myocardium_3D;
            
            roi_in_myo_t1 = roi_in_myo_t1(:,:,idx_reordered);
            remote_in_myo_t1 = remote_in_myo_t1(:,:,idx_reordered);
            roi_t1 = roi_t1(:,:,idx_reordered);
            remote_t1 = remote_t1(:,:,idx_reordered);
            t1 = t1(:,:,idx_reordered);
            myo_t1 = myo_t1(:,:,idx_reordered);
            
            % Pre-QC
            [roi_in_myo_t1_new, remote_in_myo_t1_new, roi_t1_new, remote_t1_new, t1_new, myo_t1_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_t1, remote_in_myo_t1,...
                roi_t1, remote_t1, t1, myo_t1);
            
            % remove edges of MI region
            roi_edg_t1_new = zeros(size(roi_in_myo_t1_new));
            for i = 1:size(roi_in_myo_t1_new, 3)
                roi_edg_t1_new(:,:,i)  = edge(squeeze(roi_in_myo_t1_new(:,:,i)),'Canny');
            end
            roi_rimmed_t1_new = (roi_in_myo_t1_new - roi_edg_t1_new)>0;
            [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_rimmed_t1_new);
            [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1_new);
            
            % FF
            myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
            
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat'));
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
            slc_array_ff = slc_array;
            
            ff_map = cell(1, length(glob_names));
            for f = 1:length(ff_map)
                ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
            end
            
            % convert ff_map to matrix
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{2}.fwmc_ff, 2), length(ff_map));
            for f = 1:length(ff_map)
               ff(:,:,f) = ff_map{f}.fwmc_ff; 
            end
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            ff = ff(:,:,idx_reordered);
            myo_ff = myo_ff(:,:,idx_reordered);
            remote_ff = remote_ff(:,:,idx_reordered);
            roi_ff = roi_ff(:,:,idx_reordered);
            remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
            roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            
            % Pre-QC
            [roi_in_myo_ff_new, remote_in_myo_ff_new, roi_ff_new, remote_ff_new, ff_new, myo_ff_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_ff, remote_in_myo_ff,...
                roi_ff, remote_ff, ff, myo_ff);
            
            % remove edges of MI region
            roi_edg_ff_new = zeros(size(roi_in_myo_ff_new));
            for i = 1:size(roi_in_myo_ff_new, 3)
                roi_edg_ff_new(:,:,i)  = edge(squeeze(roi_in_myo_ff_new(:,:,i)),'Canny');
            end
            roi_rimmed_ff_new = (roi_in_myo_ff_new - roi_edg_ff_new)>0;
            roi_ff_new = roi_ff_new .* roi_rimmed_ff_new;


            
            [row_roi, col_roi, v_roi] = find(roi_rimmed_ff_new);
            [row_remote, col_remote, v_remote] = find(remote_in_myo_ff_new);
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
            for f = 1:length(r2star_map)
               r2star(:,:,f) = r2star_map{f}.fwmc_r2star; 
            end
            
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;

            r2star = r2star(:,:,idx_reordered);
            myo_r2star = myo_r2star(:,:,idx_reordered);
            remote_r2star = remote_r2star(:,:,idx_reordered);
            roi_r2star = roi_r2star(:,:,idx_reordered);
            remote_in_myo_r2star = remote_in_myo_r2star(:,:,idx_reordered);
            roi_in_myo_r2star = roi_in_myo_r2star(:,:,idx_reordered);
            
            % Pre-QC
            [roi_in_myo_r2star_new, remote_in_myo_r2star_new, roi_r2star_new, remote_r2star_new, r2star_new, myo_r2star_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_r2star, remote_in_myo_r2star,...
                roi_r2star, remote_r2star, r2star, myo_r2star);
            
            % remove edges of MI region
            roi_edg_r2star_new = zeros(size(roi_in_myo_r2star_new));
            for i = 1:size(roi_in_myo_r2star_new, 3)
                roi_edg_r2star_new(:,:,i)  = edge(squeeze(roi_in_myo_r2star_new(:,:,i)),'Canny');
            end
            roi_rimmed_r2star_new = (roi_in_myo_r2star_new - roi_edg_r2star_new)>0;
            roi_r2star_new = roi_rimmed_r2star_new .* roi_r2star_new;

            [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_rimmed_r2star_new);
            [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star_new);

            
            roi_ff_array = zeros(1, length(row_roi));
            remote_ff_array = zeros(1, length(row_remote));
            roi_r2star_array = zeros(1, length(row_roi_r2star));
            remote_r2star_array = zeros(1, length(row_remote_r2star));
            roi_t1_array = zeros(1, length(row_roi_t1));
            remote_t1_array = zeros(1, length(row_remote_t1));
            
            for fff = 1:length(row_roi)
               roi_ff_array(fff) = roi_ff_new(row_roi(fff), col_roi(fff));
            end
            
            for fff = 1:length(row_remote)
                remote_ff_array(fff) = remote_ff_new(row_remote(fff), col_remote(fff));
            end
            
            for fff = 1:length(row_roi_r2star)
                roi_r2star_array(fff) = roi_r2star_new(row_roi_r2star(fff), col_roi_r2star(fff));
            end
            
            for fff = 1:length(row_remote_r2star)
                remote_r2star_array(fff) = remote_r2star_new(row_remote_r2star(fff), col_remote_r2star(fff));
            end
            
            for fff = 1:length(row_roi_t1)
                roi_t1_array(fff) = roi_t1_new(row_roi_t1(fff), col_roi_t1(fff));
            end
            
            for fff = 1:length(row_remote_t1)
                remote_t1_array(fff) = remote_t1_new(row_remote_t1(fff), col_remote_t1(fff));
            end
            
            roi_ff_array(roi_ff_array < 0) = 0;
            roi_ff_array(roi_ff_array > 100) = 100;
            remote_ff_array(remote_ff_array < 0) = 0;
            remote_ff_array(remote_ff_array > 100) = 100;
            
            roi_r2star_array(roi_r2star_array > 100) = 100;
            remote_r2star_array(remote_r2star_array > 100) = 100;

            
             % XZ 09/20/2023
            roi_ff_new(roi_ff_new < 0) = 0;
            roi_ff_new(roi_ff_new > 100) = 100;
            roi_r2star_new(roi_r2star_new > 100) = 100;

            thresh = mean(remote_r2star_array) + 2*std(remote_r2star_array);
            fib_perc = zeros(1, size(roi_ff_new, 3));
            mi_perc = zeros(1, size(roi_ff_new, 3));
            mi_pix = zeros(1, size(roi_ff_new, 3));
            ff_pix = zeros(1, size(roi_ff_new, 3));
            r2star_pix = zeros(1, size(roi_ff_new, 3));
            union_pix = zeros(1, size(roi_ff_new, 3));
            intercept_pix = zeros(1, size(roi_ff_new, 3));

            if size(roi_ff_new, 3) <= 3
                roi_ff_thresh = roi_ff_new(:) > 6;
                ff_pix = sum(roi_ff_thresh);

                roi_r2star_thresh = roi_r2star_new(:) > thresh;
                r2star_pix = sum(roi_r2star_thresh);

                roi_union = roi_ff_thresh | roi_r2star_thresh;
                union_pix = sum(roi_union);

                roi_intercept = roi_ff_thresh & roi_r2star_thresh;
                intercept_pix = sum(roi_intercept);

                fib_perc = (1 - sum(roi_union(:)) / sum(roi_rimmed_ff_new(:)))*100;

                mi_perc = sum(roi_in_myo_ff(:)) / sum(myo_ff_new(:));

                mi_pix = sum(roi_in_myo_ff(:));

            else
                for xx = 1:(size(roi_ff_new, 3)-2)
                    roi_ff_thresh = roi_ff_new(:,:,xx:xx+2) > 6;
                    ff_pix(xx) = sum(sum(roi_ff_thresh(:)));

                    roi_r2star_thresh = roi_r2star_new(:,:,xx:xx+2) > thresh;
                    r2star_pix(xx) = sum(sum(roi_r2star_thresh(:)));

                    roi_union = roi_ff_thresh | roi_r2star_thresh;
                    union_pix(xx) = sum(roi_union(:));

                    roi_intercept = roi_ff_thresh & roi_r2star_thresh;
                    intercept_pix(xx) = sum(roi_intercept(:));

                    fib_perc(xx) = (1 - sum(roi_union(:)) / sum(sum(sum(roi_rimmed_ff_new(:,:,xx:xx+2)))))*100;

                    mi_perc(xx) = sum(sum(sum(roi_in_myo_ff(:,:,xx:xx+2)))) / sum(sum(sum(myo_ff_new(:,:,xx:xx+2))));

                    mi_pix(xx) = sum(sum(sum(roi_in_myo_ff(:,:,xx:xx+2))));
                end

            end

            

            ff_new(ff_new < 0) = 1e-8;
            ff_new(ff_new > 100) = 100;
            r2star_new(r2star_new<0)=1e-8;
            r2star_new(r2star_new>200)=200;


            metrics(n).TimePoints(tp_count).time_point = time_point;

            metrics(n).TimePoints(tp_count).SliceAnalysis = struct;
            metrics(n).TimePoints(tp_count).HeteroAnalysis = struct;

            temp_roi =  t1_new .* roi_rimmed_t1_new;
            temp_remote = remote_t1_new;

            num_bins = 56;
            if size(roi_ff_new, 3) <= 3

                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(t1_new .* roi_rimmed_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(t1_new .* remote_in_myo_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(t1_new .* roi_rimmed_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(t1_new .* remote_in_myo_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(t1_new .* roi_rimmed_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(t1_new .* remote_in_myo_t1_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(t1_new .* roi_rimmed_t1_new))-3;
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(t1_new .* remote_in_myo_t1_new))-3;

                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(ff_new .* roi_rimmed_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(ff_new .* remote_in_myo_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(ff_new .* roi_rimmed_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(ff_new .* remote_in_myo_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi = skewness(nonzeros(ff_new .* roi_rimmed_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_remote = skewness(nonzeros(ff_new .* remote_in_myo_ff_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi = kurtosis(nonzeros(ff_new .* roi_rimmed_ff_new))-3;
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_remote = kurtosis(nonzeros(ff_new .* remote_in_myo_ff_new))-3;

                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi = mean(nonzeros(r2star_new .* roi_rimmed_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_remote = mean(nonzeros(r2star_new .* remote_in_myo_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi = std(nonzeros(r2star_new .* roi_rimmed_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_remote = std(nonzeros(r2star_new .* remote_in_myo_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi = skewness(nonzeros(r2star_new .* roi_rimmed_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_remote = skewness(nonzeros(r2star_new .* remote_in_myo_r2star_new));
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi = kurtosis(nonzeros(r2star_new .* roi_rimmed_r2star_new))-3;
                metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_remote = kurtosis(nonzeros(r2star_new .* remote_in_myo_r2star_new))-3;

                bipolar = temp_roi(:) - mean(nonzeros(temp_remote(:)));
                weighted_map = double(bipolar < 0) .* 2 .* roi_rimmed_t1_new(:) .* bipolar + double(bipolar >= 0) .* roi_rimmed_t1_new(:) .* bipolar;

                lb = 2*(400 - mean(nonzeros(temp_remote(:))));
                ub = 1800 - mean(nonzeros(temp_remote(:)));
                weighted_map(weighted_map < lb) = lb;
                weighted_map(weighted_map > ub) = ub;

                %temp_norm_roi = uint8((temp_roi(:,:,slc) - 800)./ (1800-800) * 256);
                %temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);

                temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
                temp_norm_remote = uint8((temp_remote(:) - 800)./ (1800-800) * 256);

                non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);
                t1_norm_roi_slc = temp_norm_roi;
                t1_norm_roi_slc(t1_norm_roi_slc ==  non_roi_value) = [];
                figure(); imhist(t1_norm_roi_slc, num_bins);
                p_roi = imhist(t1_norm_roi_slc, num_bins);
                nonZeros = find(p_roi);
                len = length((nonZeros));
                pNonZeros = zeros(1,len);

                for i = 1:len
                    pNonZeros(i) = p_roi(nonZeros(i));
                end

                % normalize pNonZeros so that sum(p) is one.
                pNonZeros = pNonZeros ./ sum(p_roi);
                E_roi = -sum(pNonZeros.*log2(pNonZeros));


                t1_norm_remote_slc = temp_norm_remote;
                t1_norm_remote_slc(t1_norm_remote_slc == 0) = [];
                figure(); imhist(t1_norm_remote_slc, num_bins);
                p_remote = imhist(t1_norm_remote_slc, num_bins);

                nonZeros = find(p_remote);
                len = length((nonZeros));
                pNonZeros = zeros(1,len);

                for i = 1:len
                    pNonZeros(i) = p_remote(nonZeros(i));
                end

                % normalize pNonZeros so that sum(p) is one.
                pNonZeros = pNonZeros ./ sum(p_remote);
                E_remote = -sum(pNonZeros.*log2(pNonZeros));


                metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;


            else

                for slc = 1:(size(roi_in_myo_ff,3)-2)
                    %for slc = 3:3

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(t1_new(:,:,slc:slc+2) .* roi_rimmed_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(t1_new(:,:,slc:slc+2) .* remote_in_myo_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(t1_new(:,:,slc:slc+2) .* roi_rimmed_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(t1_new(:,:,slc:slc+2) .* remote_in_myo_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(t1_new(:,:,slc:slc+2) .* roi_rimmed_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(t1_new(:,:,slc:slc+2) .* remote_in_myo_t1_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(t1_new(:,:,slc:slc+2) .* roi_rimmed_t1_new(:,:,slc:slc+2)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(t1_new(:,:,slc:slc+2) .* remote_in_myo_t1_new(:,:,slc:slc+2)))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(ff_new(:,:,slc:slc+2) .* roi_rimmed_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(ff_new(:,:,slc:slc+2) .* remote_in_myo_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(ff_new(:,:,slc:slc+2) .* roi_rimmed_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(ff_new(:,:,slc:slc+2) .* remote_in_myo_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi = skewness(nonzeros(ff_new(:,:,slc:slc+2) .* roi_rimmed_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_remote = skewness(nonzeros(ff_new(:,:,slc:slc+2) .* remote_in_myo_ff_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi = kurtosis(nonzeros(ff_new(:,:,slc:slc+2) .* roi_rimmed_ff_new(:,:,slc:slc+2)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_remote = kurtosis(nonzeros(ff_new(:,:,slc:slc+2) .* remote_in_myo_ff_new(:,:,slc:slc+2)))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi = mean(nonzeros(r2star_new(:,:,slc:slc+2) .* roi_rimmed_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_remote = mean(nonzeros(r2star_new(:,:,slc:slc+2) .* remote_in_myo_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi = std(nonzeros(r2star_new(:,:,slc:slc+2) .* roi_rimmed_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_remote = std(nonzeros(r2star_new(:,:,slc:slc+2) .* remote_in_myo_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi = skewness(nonzeros(r2star_new(:,:,slc:slc+2) .* roi_rimmed_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_remote = skewness(nonzeros(r2star_new(:,:,slc:slc+2) .* remote_in_myo_r2star_new(:,:,slc:slc+2)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi = kurtosis(nonzeros(r2star_new(:,:,slc:slc+2) .* roi_rimmed_r2star_new(:,:,slc:slc+2)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_remote = kurtosis(nonzeros(r2star_new(:,:,slc:slc+2) .* remote_in_myo_r2star_new(:,:,slc:slc+2)))-3;



                    % Hard-coded for T1 mapping

                    bipolar = temp_roi(:,:,slc:slc+2) - mean(nonzeros(temp_remote(:,:,slc:slc+2)));
                    weighted_map = double(bipolar < 0) .* 2 .* roi_rimmed_t1_new(:,:,slc:slc+2) .* bipolar + double(bipolar >= 0) .* roi_rimmed_t1_new(:,:,slc:slc+2) .* bipolar;

                    lb = 2*(400 - mean(nonzeros(temp_remote(:,:,slc:slc+2))));
                    ub = 1800 - mean(nonzeros(temp_remote(:,:,slc:slc+2)));
                    weighted_map(weighted_map < lb) = lb;
                    weighted_map(weighted_map > ub) = ub;

                    %temp_norm_roi = uint8((temp_roi(:,:,slc) - 800)./ (1800-800) * 256);
                    %temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);

                    temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
                    temp_norm_remote = uint8((temp_remote(:,:,slc:slc+2) - 800)./ (1800-800) * 256);

                    non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);
                    t1_norm_roi_slc = temp_norm_roi;
                    t1_norm_roi_slc(t1_norm_roi_slc ==  non_roi_value) = [];
                    figure(); imhist(t1_norm_roi_slc, num_bins);
                    p_roi = imhist(t1_norm_roi_slc, num_bins);
                    nonZeros = find(p_roi);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_roi(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_roi);
                    E_roi = -sum(pNonZeros.*log2(pNonZeros));


                    t1_norm_remote_slc = temp_norm_remote;
                    t1_norm_remote_slc(t1_norm_remote_slc == 0) = [];
                    figure(); imhist(t1_norm_remote_slc, num_bins);
                    p_remote = imhist(t1_norm_remote_slc, num_bins);

                    nonZeros = find(p_remote);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_remote(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_remote);
                    E_remote = -sum(pNonZeros.*log2(pNonZeros));


                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                end
            end


            name_tp = cat(2, name, '_', time_point);

            for_analysis_rim(n).metrics(tp).time_point = time_point;
            for_analysis_rim(n).metrics(tp).mean_roi_t1 = mean(nonzeros(roi_t1_array));
            for_analysis_rim(n).metrics(tp).sd_roi_t1 = std(nonzeros(roi_t1_array));
            for_analysis_rim(n).metrics(tp).mean_remote_t1 = mean(nonzeros(remote_t1_array));
            for_analysis_rim(n).metrics(tp).sd_remote_t1 = std(nonzeros(remote_t1_array));

            for_analysis_rim(n).metrics(tp).mean_roi_ff = mean(roi_ff_array);
            for_analysis_rim(n).metrics(tp).sd_roi_ff = std(roi_ff_array);
            for_analysis_rim(n).metrics(tp).mean_remote_ff = mean(remote_ff_array);
            for_analysis_rim(n).metrics(tp).sd_remote_ff = std(remote_ff_array);
            
            for_analysis_rim(n).metrics(tp).mean_roi_r2star = mean(nonzeros(roi_r2star_array));
            for_analysis_rim(n).metrics(tp).sd_roi_r2star = std(nonzeros(roi_r2star_array));
            for_analysis_rim(n).metrics(tp).mean_remote_r2star = mean(nonzeros(remote_r2star_array));
            for_analysis_rim(n).metrics(tp).sd_remote_r2star = std(nonzeros(remote_r2star_array));

                        
            data_storage_rim(n).data(tp).time_point = time_point;
            data_storage_rim(n).data(tp).roi_ff_array = roi_ff_array;
            data_storage_rim(n).data(tp).remote_ff_array = remote_ff_array;
            data_storage_rim(n).data(tp).roi_r2star_array = roi_r2star_array;
            data_storage_rim(n).data(tp).remote_r2star_array = remote_r2star_array;
            data_storage_rim(n).data(tp).roi_t1_array = roi_t1_array;
            data_storage_rim(n).data(tp).remote_t1_array = remote_t1_array;

            data_storage_rim(n).data(tp).fib_perc = fib_perc;
            data_storage_rim(n).data(tp).mi_perc = mi_perc;
            data_storage_rim(n).data(tp).mi_pix = mi_pix;
            
            data_storage_rim(n).data(tp).ff_pix = ff_pix;
            data_storage_rim(n).data(tp).r2star_pix = r2star_pix;
            data_storage_rim(n).data(tp).union_pix = union_pix;
            data_storage_rim(n).data(tp).intercept_pix = intercept_pix;

            tp_count = tp_count + 1;
        end
    end
    close all;
end
