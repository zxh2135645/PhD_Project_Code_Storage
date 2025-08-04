clear all;
close all;

addpath('../function/');
addpath('../GUI_n_Analysis_for_LRT/');
addpath('../T1NFF/');

% ================================ Identify your major folder, in this case
% ROC_analysis/
base_dir = uigetdir; % more generic -> % ROC_analysis
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = ExtractNames(folder_glob);

% time_points = {'D5', 'D7'};
% time_points = {'D6', 'D8'};
time_points = {'WK8', 'WK8+2'};
methods = {'LRT', 'ConvCMR'};

OutputPath = GetFullPath(cat(2, base_dir, '/Results/'));
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end


sequence_label = {'LGE', 'EGE', 'T2_MULTIECHO', 'LGE-MVO'};
sequence_label_sub = {'EGE1', 'EGE2'};
%% LGE
name_check = 'LATTE_chronic';
starting_point = find(strcmp(name_check, Names),1);


for n = starting_point:starting_point
    % for n = 2:2
    name = Names{n};
    real_name_temp = strsplit(name, '_');
    real_name = real_name_temp{1};

    SubjectPath = GetFullPath(cat(2, OutputPath, '/', name, '/'));
    if ~exist(SubjectPath, 'dir')
        mkdir(SubjectPath);
    end

    load(cat(2, SubjectPath, 'Heuristics.mat'));
    shift_idx = Heuristics.(name).shift_idx;
    idx_mag = Heuristics.(name).idx_mag;
    matching_cell = Heuristics.(name).matching_cell;

    for tp = 1:length(time_points)
        %for tp = 2:2

        time_point = time_points{end-tp+1};

        contour_dir_lrt_check = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{1}, '/', sequence_label{1}, '/');
        if ~exist(cat(2, contour_dir_lrt_check, sequence_label{1}, '_vol_img_3D.mat'))

            disp(['skip: ', time_point])

        else
            TimePtPath = GetFullPath(cat(2, SubjectPath, '/', time_point, '/'));
            if ~exist(TimePtPath, 'dir')
                mkdir(TimePtPath);
            end

            LRTPath = GetFullPath(cat(2, TimePtPath, '/', methods{1}, '/'));
            if ~exist(LRTPath, 'dir')
                mkdir(LRTPath);
            end

            CMRPath = GetFullPath(cat(2, TimePtPath, '/', methods{2}, '/'));
            if ~exist(CMRPath, 'dir')
                mkdir(CMRPath);
            end

            AHA_segs = struct;

            for ll = 1:length(sequence_label)
            % for ll = 2:2
                label = sequence_label{ll};

                if strcmp(label, 'LGE-MVO')
                    contour_dir_lrt = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{1}, '/', sequence_label{1}, '/');
                else
                    contour_dir_lrt = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{1}, '/', label, '/');
                end

                if (isempty(matching_cell) && strcmp(label, 'T2_MULTIECHO'))
                    t2s_six_segments_mean_lrt = [];
                    Segmentpix_lrt = [];
                    stats_lrt = [];
                    Mask_index_lrt = [];

                    t2s_six_segments_mean_cmr = [];
                    Segmentpix_cmr = [];
                    stats_cmr = [];
                    Mask_index_cmr = [];

                    AHA_segs(ll).t2s_six_segments_mean_lrt = t2s_six_segments_mean_lrt;
                    AHA_segs(ll).t2s_six_segments_mean_cmr = t2s_six_segments_mean_cmr;

                    AHA_segs(ll).Segmentpix_cmr = Segmentpix_cmr;
                    AHA_segs(ll).stats_cmr = stats_cmr;
                    AHA_segs(ll).Mask_index_cmr = Mask_index_cmr;
                    AHA_segs(ll).Segmentpix_lrt = Segmentpix_lrt;
                    AHA_segs(ll).stats_lrt = stats_lrt;
                    AHA_segs(ll).Mask_index_lrt = Mask_index_lrt;
                    continue; % only works for T2_MULTIECHO
                elseif (strcmp(time_point, 'WK8') || strcmp(time_point, 'WK8+2')) && strcmp(label, 'EGE')
                    t2s_six_segments_mean_lrt_early = [];
                    Segmentpix_lrt_early = [];
                    stats_lrt_early = [];
                    Mask_index_lrt_early = [];
                    t2s_six_segments_mean_lrt_pseudo = [];
                    Segmentpix_lrt_pseudo = [];
                    stats_lrt_pseudo = [];
                    Mask_index_lrt_pseudo = [];

                    t2s_six_segments_mean_cmr_early = [];
                    Segmentpix_cmr_early = [];
                    stats_cmr_early = [];
                    Mask_index_cmr_early = [];
                    t2s_six_segments_mean_cmr_pseudo = [];
                    Segmentpix_cmr_pseudo = [];
                    stats_cmr_pseudo = [];
                    Mask_index_cmr_pseudo = [];

                    AHA_segs(ll).t2s_six_segments_mean_lrt_early = t2s_six_segments_mean_lrt_early;
                    AHA_segs(ll).t2s_six_segments_mean_cmr_early = t2s_six_segments_mean_cmr_early;

                    AHA_segs(ll).Segmentpix_cmr_early = Segmentpix_cmr_early;
                    AHA_segs(ll).stats_cmr_early = stats_cmr_early;
                    AHA_segs(ll).Mask_index_cmr_early = Mask_index_cmr_early;
                    AHA_segs(ll).Segmentpix_lrt_early = Segmentpix_lrt_early;
                    AHA_segs(ll).stats_lrt_early = stats_lrt_early;
                    AHA_segs(ll).Mask_index_lrt_early = Mask_index_lrt_early;

                    AHA_segs(ll).t2s_six_segments_mean_lrt_pseudo = t2s_six_segments_mean_lrt_pseudo;
                    AHA_segs(ll).t2s_six_segments_mean_cmr_pseudo = t2s_six_segments_mean_cmr_pseudo;

                    AHA_segs(ll).Segmentpix_cmr_pseudo = Segmentpix_cmr_pseudo;
                    AHA_segs(ll).stats_cmr_pseudo = stats_cmr_pseudo;
                    AHA_segs(ll).Mask_index_cmr_pseudo = Mask_index_cmr_pseudo;
                    AHA_segs(ll).Segmentpix_lrt_pseudo = Segmentpix_lrt_pseudo;
                    AHA_segs(ll).stats_lrt_pseudo = stats_lrt_pseudo;
                    AHA_segs(ll).Mask_index_lrt_pseudo = Mask_index_lrt_pseudo;
                    continue;
                elseif (strcmp(time_point, 'WK8') || strcmp(time_point, 'WK8+2')) && strcmp(label, 'LGE-MVO')
                    t2s_six_segments_mean_lgemvo_lrt = [];
                    Segmentpix_lgemvo_lrt = [];
                    stats_lgemvo_lrt = [];
                    Mask_index_lgemvo_lrt = [];

                    t2s_six_segments_mean_lgemvo_cmr = [];
                    Segmentpix_lgemvo_cmr = [];
                    stats_lgemvo_cmr = [];
                    Mask_index_lgemvo_cmr = [];

                    AHA_segs(ll).t2s_six_segments_mean_lgemvo_lrt = t2s_six_segments_mean_lgemvo_lrt;
                    AHA_segs(ll).t2s_six_segments_mean_lgemvo_cmr = t2s_six_segments_mean_lgemvo_cmr;

                    AHA_segs(ll).Segmentpix_lgemvo_cmr = Segmentpix_lgemvo_cmr;
                    AHA_segs(ll).stats_lgemvo_cmr = stats_lgemvo_cmr;
                    AHA_segs(ll).Mask_index_lgemvo_cmr = Mask_index_lgemvo_cmr;
                    AHA_segs(ll).Segmentpix_lgemvo_lrt = Segmentpix_lgemvo_lrt;
                    AHA_segs(ll).stats_lgemvo_lrt = stats_lgemvo_lrt;
                    AHA_segs(ll).Mask_index_lgemvo_lrt = Mask_index_lgemvo_lrt;
                    continue; % only works for T2_MULTIECHO
                end

                if strcmp(label, 'LGE-MVO')
                    LGE_LRT = load(cat(2, contour_dir_lrt, sequence_label{1}, '_vol_img_3D.mat'));
                    sliceloc_LRT = load(cat(2, contour_dir_lrt, sequence_label{1}, '_SliceLoc.mat'));
                else
                    LGE_LRT = load(cat(2, contour_dir_lrt, label, '_vol_img_3D.mat'));
                    sliceloc_LRT = load(cat(2, contour_dir_lrt, label, '_SliceLoc.mat'));
                end

                myo_LRT = load(cat(2, contour_dir_lrt, '/Myocardium/mask_myocardium.mat'));
                freeroi_LRT = load(cat(2, contour_dir_lrt, '/freeROI/freeROI.mat'));
                heart_LRT = load(cat(2, contour_dir_lrt, '/Heart/mask_heart.mat'));
                blood_LRT = load(cat(2, contour_dir_lrt, '/BloodPool/mask_blood.mat'));
                noReflow_LRT = load(cat(2, contour_dir_lrt, '/noReflowArea/noReflow.mat'));
                excludeArea_LRT = load(cat(2, contour_dir_lrt, '/excludeArea/excludeArea.mat'));
                myoRef_LRT = load(cat(2, contour_dir_lrt, '/MyoReference/myoRef.mat'));

                % LRT
                vol_img_3D_lrt = LGE_LRT.vol_img_3D;
                slc_array_lrt = sliceloc_LRT.slc_array;
                mask_myocardium_3D_lrt = myo_LRT.mask_myocardium_3D;
                freeROIMask_3D_lrt = freeroi_LRT.freeROIMask_3D;
                mask_heart_3D_lrt = heart_LRT.mask_heart_3D;
                mask_blood_3D_lrt = blood_LRT.mask_blood_3D;
                noReflowMask_3D_lrt = noReflow_LRT.noReflowMask_3D;
                excludeMask_3D_lrt = excludeArea_LRT.excludeMask_3D;
                myoRefMask_3D_lrt = myoRef_LRT.myoRefMask_3D;


                slc_size = size(vol_img_3D_lrt, 3);
                if strcmp(label, 'EGE') 
                    [idx_lrt] = find(any(reshape(noReflowMask_3D_lrt, [], slc_size)));
                elseif strcmp(label, 'LGE') || strcmp(label, 'LGE-MVO')
                    [idx_lrt] = find(any(reshape(freeROIMask_3D_lrt, [], slc_size)));
                elseif strcmp(label, 'T2_MULTIECHO')
                    [idx_lrt] = find(any(reshape(myoRefMask_3D_lrt, [], slc_size)));
                end

                vol_img_3D_lrt = vol_img_3D_lrt(:,:,idx_lrt);
                slc_array_lrt = slc_array_lrt(idx_lrt);
                mask_myocardium_3D_lrt = mask_myocardium_3D_lrt(:,:,idx_lrt);
                freeROIMask_3D_lrt = freeROIMask_3D_lrt(:,:,idx_lrt);
                mask_heart_3D_lrt = mask_heart_3D_lrt(:,:,idx_lrt);
                mask_blood_3D_lrt = mask_blood_3D_lrt(:,:,idx_lrt);
                noReflowMask_3D_lrt = noReflowMask_3D_lrt(:,:,idx_lrt);
                excludeMask_3D_lrt = excludeMask_3D_lrt(:,:,idx_lrt);
                myoRefMask_3D_lrt = myoRefMask_3D_lrt(:,:,idx_lrt);


                if strcmp(label, 'T2_MULTIECHO')
                    LGE_LRT_te = load(cat(2, contour_dir_lrt, label, '_vol_img_3D_te.mat'));
                    freeroi_LRT_te = load(cat(2, contour_dir_lrt, '/freeROI/freeROI_te.mat'));
                    excludeArea_LRT_te = load(cat(2, contour_dir_lrt, '/excludeArea/excludeMask_te.mat'));
                    myoRef_LRT_te = load(cat(2, contour_dir_lrt, '/MyoReference/myoRef_te.mat'));

                    vol_img_3D_te_lrt = LGE_LRT_te.vol_img_3D_te;
                    freeROIMask_3D_te_lrt = freeroi_LRT_te.freeROIMask_3D_te;
                    excludeMask_3D_te_lrt = excludeArea_LRT_te.excludeMask_3D_te;
                    myoRefMask_3D_te_lrt = myoRef_LRT_te.myoRef_3D_te;

                    vol_img_3D_te_lrt = vol_img_3D_te_lrt(:,:,idx_lrt);
                    freeROIMask_3D_te_lrt = freeROIMask_3D_te_lrt(:,:,idx_lrt);
                    excludeMask_3D_te_lrt = excludeMask_3D_te_lrt(:,:,idx_lrt);
                    myoRefMask_3D_te_lrt = myoRefMask_3D_te_lrt(:,:,idx_lrt);
                end

                %%
                if strcmp(label, 'EGE')
                    % vol_img_3D_lrt_early = vol_img_3D_lrt(:,:,2);
                    % mask_myocardium_3D_lrt_early = mask_myocardium_3D_lrt(:,:,2);
                    % noReflowMask_3D_lrt_early = noReflowMask_3D_lrt(:,:,2);
                    % myoRefMask_3D_lrt_early = myoRefMask_3D_lrt(:,:,2);
                    % mask_blood_3D_lrt_early = mask_blood_3D_lrt(:,:,2);
                    % 
                    % vol_img_3D_lrt_pseudo = vol_img_3D_lrt(:,:,13);
                    % mask_myocardium_3D_lrt_pseudo = mask_myocardium_3D_lrt(:,:,13);
                    % noReflowMask_3D_lrt_pseudo = noReflowMask_3D_lrt(:,:,13);
                    % myoRefMask_3D_lrt_pseudo = myoRefMask_3D_lrt(:,:,13);
                    % mask_blood_3D_lrt_pseudo = mask_blood_3D_lrt(:,:,13);

                    % thresh = mean(nonzeros(vol_img_3D_lrt_early .* myoRefMask_3D_lrt_early)) + ...
                    %     5*std(nonzeros(vol_img_3D_lrt_early .* myoRefMask_3D_lrt_early));
                    % roi_early = vol_img_3D_lrt_early >= thresh;
                    % noReflow_lrt_early = roi_early .* mask_myocardium_3D_lrt_early .* noReflowMask_3D_lrt_early;
                    % 
                    % thresh = mean(nonzeros(vol_img_3D_lrt_pseudo .* myoRefMask_3D_lrt_pseudo)) + ...
                    %     5*std(nonzeros(vol_img_3D_lrt_pseudo .* myoRefMask_3D_lrt_pseudo));
                    % roi_pseudo = vol_img_3D_lrt_pseudo >= thresh;
                    % noReflow_lrt_pseudo = roi_pseudo .* mask_myocardium_3D_lrt_pseudo .* noReflowMask_3D_lrt_pseudo;


                    vol_img_3D_lrt_pseudo = vol_img_3D_lrt(:,:,1);
                    x = zeros(size(vol_img_3D_lrt_pseudo, 3),1);
                    y = zeros(size(vol_img_3D_lrt_pseudo, 3),1);

                    InsertionPtLRTPath = cat(2, LRTPath, '/InsertionPts/');
                    if ~exist(InsertionPtLRTPath, 'dir')
                        mkdir(InsertionPtLRTPath);
                    end

                    InsertionPt_LRT_f = cat(2, InsertionPtLRTPath, '/', 'pts_EGE.mat');

                    if ~exist(InsertionPt_LRT_f)
                        for slc = 1:size(vol_img_3D_lrt_pseudo, 3)
                            figure(104);
                            imagesc(vol_img_3D_lrt_pseudo); colormap gray; axis image;
                            title('Draw Insertion point (LRT-EGE)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 104
                        end
                        ref_pts_lrt.x = x;
                        ref_pts_lrt.y = y;
                        save(InsertionPt_LRT_f, 'ref_pts_lrt');
                    else
                        load(InsertionPt_LRT_f);
                    end

                    clear BaseGroove_lrt
                    for slc = 1:size(vol_img_3D_lrt_pseudo,3)
                        C = regionprops(mask_blood_3D_lrt(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_lrt.x(slc) - x_centroid, ref_pts_lrt.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_lrt(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_lrt(slc) = -Groove_temp -180;
                        else
                            BaseGroove_lrt(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end


                    roi_lrt = zeros(size(vol_img_3D_lrt));
                    noReflow_lrt = zeros(size(vol_img_3D_lrt));
                    % thresh_array = [932 1861 1745 1627 1520 1418 1319 1257 1214 1119 1116 1099 1044 1010 976];
                    % CHILI D8? CINNAMON????
                    if strcmp(name, 'XXX')
                        thresh_array = [932 1861 1745 1627 1520 1418 1319 1257 1214 1119 1116 1099 1044 1010 976];
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            thresh = mean(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc))) + ...
                                5*std(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc)));
                            % disp(thresh)
                            thresh = thresh_array(slc);
                            roi_lrt(:,:,slc) = vol_img_3D_lrt(:,:,slc) <= thresh;
                            noReflow_lrt(:,:,slc) = roi_lrt(:,:,slc) .* mask_myocardium_3D_lrt(:,:,slc) .* noReflowMask_3D_lrt(:,:,slc) .* ~excludeMask_3D_lrt(:,:,slc);
                        end
                    else
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            thresh = mean(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc))) + ...
                                5*std(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc)));
                            % disp(thresh)
                            % thresh = thresh_array(slc);
                            roi_lrt(:,:,slc) = vol_img_3D_lrt(:,:,slc) <= thresh;
                            noReflow_lrt(:,:,slc) = roi_lrt(:,:,slc) .* mask_myocardium_3D_lrt(:,:,slc) .* noReflowMask_3D_lrt(:,:,slc) .* ~excludeMask_3D_lrt(:,:,slc);
                        end
                    end

                    mvo_size_array = sum(reshape(noReflow_lrt, [], size(noReflow_lrt,3)),1);
                    [mvo_size_array_sorted, idx] = sort(mvo_size_array);
                    idx_pseudo = find(mvo_size_array_sorted(3) == mvo_size_array);
                    if length(idx_pseudo) > 1
                        idx_pseudo = idx_pseudo(end-1);
                    end
                    idx_early = find(mvo_size_array_sorted(14) == mvo_size_array);
                    if length(idx_early) > 1
                        idx_early = idx_early(1);
                    end

                    [Segmentpix_lrt_early, stats_lrt_early, Mask_index_lrt_early] =AHASegmentation_qsm_invivo(noReflow_lrt(:,:,idx_early), mask_myocardium_3D_lrt(:,:,idx_early),6,BaseGroove_lrt);
                    t2s_six_segments_mean_lrt_early = squeeze(stats_lrt_early(2,:,:));
                    t2s_six_segments_mean_lrt_early = t2s_six_segments_mean_lrt_early(:);

                    [Segmentpix_lrt_pseudo, stats_lrt_pseudo, Mask_index_lrt_pseudo] =AHASegmentation_qsm_invivo(noReflow_lrt(:,:,idx_pseudo), mask_myocardium_3D_lrt(:,:,idx_pseudo),6,BaseGroove_lrt);
                    t2s_six_segments_mean_lrt_pseudo = squeeze(stats_lrt_pseudo(2,:,:));
                    t2s_six_segments_mean_lrt_pseudo = t2s_six_segments_mean_lrt_pseudo(:);

                elseif strcmp(label, 'LGE') || strcmp(label, 'LGE-MVO')
                    %shift_idx = [5 4 3 2 1 6];
                    %shift_idx = [4 3 2 1 10 9 8 7 6 5];
                    vol_img_3D_lrt = vol_img_3D_lrt(:,:,shift_idx);
                    slc_array_lrt = slc_array_lrt(shift_idx);
                    mask_myocardium_3D_lrt = mask_myocardium_3D_lrt(:,:,shift_idx);
                    freeROIMask_3D_lrt = freeROIMask_3D_lrt(:,:,shift_idx);
                    mask_heart_3D_lrt = mask_heart_3D_lrt(:,:,shift_idx);
                    mask_blood_3D_lrt = mask_blood_3D_lrt(:,:,shift_idx);
                    noReflowMask_3D_lrt = noReflowMask_3D_lrt(:,:,shift_idx);
                    excludeMask_3D_lrt = excludeMask_3D_lrt(:,:,shift_idx);
                    myoRefMask_3D_lrt = myoRefMask_3D_lrt(:,:,shift_idx);

                    x = zeros(size(vol_img_3D_lrt, 3),1);
                    y = zeros(size(vol_img_3D_lrt, 3),1);

                    InsertionPtLRTPath = cat(2, LRTPath, '/InsertionPts/');
                    if ~exist(InsertionPtLRTPath, 'dir')
                        mkdir(InsertionPtLRTPath);
                    end

                    InsertionPt_LRT_f = cat(2, InsertionPtLRTPath, '/', 'pts.mat');

                    if ~exist(InsertionPt_LRT_f)
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            figure(104);
                            imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (LRT)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 104
                        end
                        ref_pts_lrt.x = x;
                        ref_pts_lrt.y = y;
                        save(InsertionPt_LRT_f, 'ref_pts_lrt');
                    else
                        load(InsertionPt_LRT_f);
                    end

                    clear BaseGroove_lrt
                    for slc = 1:size(vol_img_3D_lrt,3)
                        C = regionprops(mask_blood_3D_lrt(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_lrt.x(slc) - x_centroid, ref_pts_lrt.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_lrt(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_lrt(slc) = -Groove_temp -180;
                        else
                            BaseGroove_lrt(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end

                    [Segmentpix_lrt, stats_lrt, Mask_index_lrt] =AHASegmentation_qsm_invivo(freeROIMask_3D_lrt, mask_myocardium_3D_lrt,6,BaseGroove_lrt);
                    t2s_six_segments_mean_lrt = squeeze(stats_lrt(2,:,:));
                    t2s_six_segments_mean_lrt = t2s_six_segments_mean_lrt(:);


                    roi_lrt = zeros(size(vol_img_3D_lrt));
                    noReflow_lrt = zeros(size(vol_img_3D_lrt));

                    for slc = 1:size(vol_img_3D_lrt, 3)
                        thresh = mean(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc))) + ...
                            5*std(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc)));
                        % disp(thresh)
                        % thresh = thresh_array(slc);
                        roi_lrt(:,:,slc) = vol_img_3D_lrt(:,:,slc) <= thresh;
                        noReflow_lrt(:,:,slc) = roi_lrt(:,:,slc) .* mask_myocardium_3D_lrt(:,:,slc) .* noReflowMask_3D_lrt(:,:,slc) .* ~excludeMask_3D_lrt(:,:,slc);
                    end

                    [Segmentpix_lgemvo_lrt, stats_lgemvo_lrt, Mask_index_lgemvo_lrt] =AHASegmentation_qsm_invivo(noReflow_lrt, mask_myocardium_3D_lrt,6,BaseGroove_lrt);
                    t2s_six_segments_mean_lgemvo_lrt = squeeze(stats_lgemvo_lrt(2,:,:));
                    t2s_six_segments_mean_lgemvo_lrt = t2s_six_segments_mean_lgemvo_lrt(:);

                elseif strcmp(label, 'T2_MULTIECHO')
                    slc_size = size(mask_myocardium_3D_lrt,3);
                    % matching_cell = [5, 0];
                    shift_idx_te = matching_cell(1);
                    % shift_idx = 4;
                    vol_img_3D_lrt = circshift(vol_img_3D_lrt,shift_idx_te,3);
                    slc_array_lrt = circshift(slc_array_lrt, shift_idx_te);
                    mask_myocardium_3D_lrt = circshift(mask_myocardium_3D_lrt, shift_idx_te, 3);
                    freeROIMask_3D_lrt = circshift(freeROIMask_3D_lrt,shift_idx_te, 3);
                    mask_heart_3D_lrt = circshift(mask_heart_3D_lrt,shift_idx_te, 3);
                    mask_blood_3D_lrt = circshift(mask_blood_3D_lrt,shift_idx_te, 3);
                    noReflowMask_3D_lrt = circshift(noReflowMask_3D_lrt,shift_idx_te, 3);
                    excludeMask_3D_lrt = circshift(excludeMask_3D_lrt,shift_idx_te, 3);
                    myoRefMask_3D_lrt = circshift(myoRefMask_3D_lrt,shift_idx_te, 3);

                    vol_img_3D_te_lrt = circshift(vol_img_3D_te_lrt,shift_idx_te, 3);
                    freeROIMask_3D_te_lrt = circshift(freeROIMask_3D_te_lrt,shift_idx_te, 3);
                    excludeMask_3D_te_lrt = circshift(excludeMask_3D_te_lrt,shift_idx_te, 3);
                    myoRefMask_3D_te_lrt = circshift(myoRefMask_3D_te_lrt,shift_idx_te, 3);


                    exclude_within_myo = ~excludeMask_3D_lrt.*mask_myocardium_3D_lrt;
                    [idx_lrt2] = find(any(reshape(exclude_within_myo, [], slc_size)));

                    vol_img_3D_lrt = vol_img_3D_lrt(:,:,idx_lrt2);
                    slc_array_lrt = slc_array_lrt(idx_lrt2);
                    mask_myocardium_3D_lrt = mask_myocardium_3D_lrt(:,:,idx_lrt2);
                    freeROIMask_3D_lrt = freeROIMask_3D_lrt(:,:,idx_lrt2);
                    mask_heart_3D_lrt = mask_heart_3D_lrt(:,:,idx_lrt2);
                    mask_blood_3D_lrt = mask_blood_3D_lrt(:,:,idx_lrt2);
                    noReflowMask_3D_lrt = noReflowMask_3D_lrt(:,:,idx_lrt2);
                    excludeMask_3D_lrt = excludeMask_3D_lrt(:,:,idx_lrt2);
                    myoRefMask_3D_lrt = myoRefMask_3D_lrt(:,:,idx_lrt2);
                    vol_img_3D_te_lrt = vol_img_3D_te_lrt(:,:,idx_lrt2);
                    freeROIMask_3D_te_lrt = freeROIMask_3D_te_lrt(:,:,idx_lrt2);
                    excludeMask_3D_te_lrt = excludeMask_3D_te_lrt(:,:,idx_lrt2);
                    myoRefMask_3D_te_lrt = myoRefMask_3D_te_lrt(:,:,idx_lrt2);

                    reverse_label = matching_cell(2);
                    if reverse_label == 1
                        vol_img_3D_lrt = flip(vol_img_3D_lrt , 3);
                        slc_array_lrt = flip(slc_array_lrt);
                        mask_myocardium_3D_lrt = flip(mask_myocardium_3D_lrt, 3);
                        freeROIMask_3D_lrt = flip(freeROIMask_3D_lrt, 3);
                        mask_heart_3D_lrt = flip(mask_heart_3D_lrt, 3);
                        mask_blood_3D_lrt = flip(mask_blood_3D_lrt, 3);
                        noReflowMask_3D_lrt = flip(noReflowMask_3D_lrt, 3);
                        excludeMask_3D_lrt = flip(excludeMask_3D_lrt,3);
                        myoRefMask_3D_lrt = flip(myoRefMask_3D_lrt, 3);

                        vol_img_3D_te_lrt = flip(vol_img_3D_te_lrt, 3);
                        freeROIMask_3D_te_lrt = flip(freeROIMask_3D_te_lrt, 3);
                        excludeMask_3D_te_lrt = flip(excludeMask_3D_te_lrt, 3);
                        myoRefMask_3D_te_lrt = flip(myoRefMask_3D_te_lrt, 3);

                    end

                    hemo_lrt = zeros(size(vol_img_3D_lrt));
                    for slc = 1:size(vol_img_3D_lrt, 3)
                        thresh = mean(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc))) - ...
                            2*std(nonzeros(vol_img_3D_lrt(:,:,slc) .* myoRefMask_3D_lrt(:,:,slc)));
                        roi_temp = vol_img_3D_lrt(:,:,slc) < thresh;
                        hemo_lrt(:,:,slc) = roi_temp .* mask_myocardium_3D_lrt(:,:,slc) .* ~excludeMask_3D_lrt(:,:,slc);
                    end

                    x = zeros(size(vol_img_3D_lrt, 3),1);
                    y = zeros(size(vol_img_3D_lrt, 3),1);

                    InsertionPtLRTPath = cat(2, LRTPath, '/InsertionPts/');
                    if ~exist(InsertionPtLRTPath, 'dir')
                        mkdir(InsertionPtLRTPath);
                    end

                    InsertionPt_LRT_f = cat(2, InsertionPtLRTPath, '/', 'pts_te.mat');

                    if ~exist(InsertionPt_LRT_f)
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            figure(107);
                            imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (LRT-MULTIECHO)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 107
                        end
                        ref_pts_lrt.x = x;
                        ref_pts_lrt.y = y;
                        save(InsertionPt_LRT_f, 'ref_pts_lrt');
                    else
                        load(InsertionPt_LRT_f);
                    end

                    clear BaseGroove_lrt
                    for slc = 1:size(vol_img_3D_lrt,3)
                        C = regionprops(mask_blood_3D_lrt(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_lrt.x(slc) - x_centroid, ref_pts_lrt.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_lrt(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_lrt(slc) = -Groove_temp -180;
                        else
                            BaseGroove_lrt(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end

                    [Segmentpix_lrt, stats_lrt, Mask_index_lrt] =AHASegmentation_qsm_invivo(hemo_lrt, mask_myocardium_3D_lrt,6,BaseGroove_lrt);
                    t2s_six_segments_mean_lrt = squeeze(stats_lrt(2,:,:));
                    t2s_six_segments_mean_lrt = t2s_six_segments_mean_lrt(:);

                end

                %% CMR
                clear vol_img_3D_cmr
                if strcmp(label, 'EGE')
                    CMR = struct;
                    for lll = 1:length(sequence_label_sub)
                        label_sub = sequence_label_sub{lll};
                        time_point_theother = time_points{3 - find(strcmp(time_point, time_points))};
                        contour_dir_convcmr = cat(2, base_dir, '/ContourData/', name, '/', time_point_theother, '/', methods{2}, '/', label_sub, '/');

                        CMR(lll).vol_img_3D = load(cat(2, contour_dir_convcmr, label_sub, '_vol_img_3D.mat'));
                        CMR(lll).slc_array = load(cat(2, contour_dir_convcmr, label_sub,'_SliceLoc.mat'));
                        CMR(lll).mask_myocardium_3D = load(cat(2, contour_dir_convcmr, '/Myocardium/mask_myocardium.mat'));
                        CMR(lll).freeROIMask_3D = load(cat(2, contour_dir_convcmr, '/freeROI/freeROI.mat'));
                        CMR(lll).mask_heart_3D = load(cat(2, contour_dir_convcmr, '/Heart/mask_heart.mat'));
                        CMR(lll).mask_blood_3D = load(cat(2, contour_dir_convcmr, '/BloodPool/mask_blood.mat'));
                        CMR(lll).excludeMask_3D = load(cat(2, contour_dir_convcmr, '/excludeArea/excludeArea.mat'));
                        CMR(lll).myoRefMask_3D = load(cat(2, contour_dir_convcmr, '/MyoReference/myoRef.mat'));
                        CMR(lll).noReflowMask_3D = load(cat(2, contour_dir_convcmr, '/noReflowArea/noReflow.mat'));
                    end
                elseif strcmp(label, 'T2_MULTIECHO')
                    CMR = struct;
                    contour_dir_convcmr = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{2}, '/', label, '/');
                    CMR(1).vol_img_3D = load(cat(2, contour_dir_convcmr, label, '_vol_img_3D.mat'));
                    CMR(1).vol_img_3D_te = load(cat(2, contour_dir_convcmr, label, '_vol_img_3D_te.mat'));
                    CMR(1).slc_array = load(cat(2, contour_dir_convcmr, label, '_SliceLoc.mat'));
                    CMR(1).Index = load(cat(2, contour_dir_convcmr, label, '_index.mat'));
                    CMR(1).mask_myocardium_3D = load(cat(2, contour_dir_convcmr, '/Myocardium/mask_myocardium.mat'));
                    CMR(1).freeROIMask_3D = load(cat(2, contour_dir_convcmr, '/freeROI/freeROI.mat'));
                    CMR(1).mask_heart_3D = load(cat(2, contour_dir_convcmr, '/Heart/mask_heart.mat'));
                    CMR(1).mask_blood_3D = load(cat(2, contour_dir_convcmr, '/BloodPool/mask_blood.mat'));
                    CMR(1).excludeMask_3D = load(cat(2, contour_dir_convcmr, '/excludeArea/excludeArea.mat'));
                    CMR(1).myoRefMask_3D = load(cat(2, contour_dir_convcmr, '/MyoReference/myoRef.mat'));
                    CMR(1).noReflowMask_3D = load(cat(2, contour_dir_convcmr, '/noReflowArea/noReflow.mat'));

                    CMR(1).freeROIMask_3D_te = load(cat(2, contour_dir_convcmr, '/freeROI/freeROI_te.mat'));
                    CMR(1).myoRefMask_3D_te = load(cat(2, contour_dir_convcmr, '/MyoReference/myoRef_te.mat'));
                    CMR(1).excludeMask_3D_te = load(cat(2, contour_dir_convcmr, '/excludeArea/excludeMask_te.mat'));
                elseif strcmp(label, 'LGE')
                    CMR = struct;
                    contour_dir_convcmr = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{2}, '/', label, '/');
                    CMR(1).vol_img_3D = load(cat(2, contour_dir_convcmr, label, '_vol_img_3D.mat'));
                    CMR(1).slc_array = load(cat(2, contour_dir_convcmr, label, '_SliceLoc.mat'));
                    CMR(1).mask_myocardium_3D = load(cat(2, contour_dir_convcmr, '/Myocardium/mask_myocardium.mat'));
                    CMR(1).freeROIMask_3D = load(cat(2, contour_dir_convcmr, '/freeROI/freeROI.mat'));
                    CMR(1).mask_heart_3D = load(cat(2, contour_dir_convcmr, '/Heart/mask_heart.mat'));
                    CMR(1).mask_blood_3D = load(cat(2, contour_dir_convcmr, '/BloodPool/mask_blood.mat'));
                    CMR(1).excludeMask_3D = load(cat(2, contour_dir_convcmr, '/excludeArea/excludeArea.mat'));
                    CMR(1).myoRefMask_3D = load(cat(2, contour_dir_convcmr, '/MyoReference/myoRef.mat'));
                    CMR(1).noReflowMask_3D = load(cat(2, contour_dir_convcmr, '/noReflowArea/noReflow.mat'));
                
                elseif strcmp(label, 'LGE-MVO')
                    CMR = struct;
                    contour_dir_convcmr = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{2}, '/', sequence_label{1}, '/');
                    CMR(1).vol_img_3D = load(cat(2, contour_dir_convcmr, sequence_label{1}, '_vol_img_3D.mat'));
                    CMR(1).slc_array = load(cat(2, contour_dir_convcmr, sequence_label{1}, '_SliceLoc.mat'));
                    CMR(1).mask_myocardium_3D = load(cat(2, contour_dir_convcmr, '/Myocardium/mask_myocardium.mat'));
                    CMR(1).freeROIMask_3D = load(cat(2, contour_dir_convcmr, '/freeROI/freeROI.mat'));
                    CMR(1).mask_heart_3D = load(cat(2, contour_dir_convcmr, '/Heart/mask_heart.mat'));
                    CMR(1).mask_blood_3D = load(cat(2, contour_dir_convcmr, '/BloodPool/mask_blood.mat'));
                    CMR(1).excludeMask_3D = load(cat(2, contour_dir_convcmr, '/excludeArea/excludeArea.mat'));
                    CMR(1).myoRefMask_3D = load(cat(2, contour_dir_convcmr, '/MyoReference/myoRef.mat'));
                    CMR(1).noReflowMask_3D = load(cat(2, contour_dir_convcmr, '/noReflowArea/noReflow.mat'));
                end


                % This only works when EGE and LGE are in same matrix size, but
                % it's not always the case!!!
                %%
                clear vol_img_3D_cmr slc_array_cmr mask_myocardium_3D_cmr freeROIMask_3D_cmr mask_heart_3D_cmr mask_blood_3D_cmr excludeMask_3D_cmr myoRefMask_3D_cmr noReflowMask_3D_cmr
                clear vol_img_3D_cmr_te glob_names_cmr echo_idx_te_cmr

                % if it is 1, it means it will be cell
                cc = 1;

                vol_img_3D_cmr = CMR(cc).vol_img_3D.vol_img_3D;
                slc_array_cmr = CMR(cc).slc_array.slc_array;
                mask_myocardium_3D_cmr = CMR(cc).mask_myocardium_3D.mask_myocardium_3D;
                freeROIMask_3D_cmr = CMR(cc).freeROIMask_3D.freeROIMask_3D;
                mask_heart_3D_cmr = CMR(cc).mask_heart_3D.mask_heart_3D;
                mask_blood_3D_cmr = CMR(cc).mask_blood_3D.mask_blood_3D;
                excludeMask_3D_cmr = CMR(cc).excludeMask_3D.excludeMask_3D;
                myoRefMask_3D_cmr = CMR(cc).myoRefMask_3D.myoRefMask_3D;
                noReflowMask_3D_cmr = CMR(cc).noReflowMask_3D.noReflowMask_3D;

                if exist('vol_img_3D_cmr', 'var') && length(CMR) > 1
                    size_check = size(CMR(2).vol_img_3D.vol_img_3D,1)  ~= size(CMR(1).vol_img_3D.vol_img_3D,1) | size(CMR(2).vol_img_3D.vol_img_3D,2) ~= size(CMR(1).vol_img_3D.vol_img_3D,2);
                else
                    size_check = 0;
                end


                for cc = 1:length(CMR)
                    if cc > 1
                        if size_check
                            slc_array_cmr = cat(2, slc_array_cmr, CMR(cc).slc_array.slc_array);

                            vol_img_3D_cmr_cell = {};
                            mask_myocardium_3D_cmr_cell = {};
                            freeROIMask_3D_cmr_cell = {};
                            mask_heart_3D_cmr_cell = {};
                            mask_blood_3D_cmr_cell = {};
                            excludeMask_3D_cmr_cell = {};
                            myoRefMask_3D_cmr_cell = {};
                            noReflowMask_3D_cmr_cell = {};

                            vol_img_3D_cmr_cell{1} = vol_img_3D_cmr;
                            mask_myocardium_3D_cmr_cell{1} = mask_myocardium_3D_cmr;
                            freeROIMask_3D_cmr_cell{1} = freeROIMask_3D_cmr;
                            mask_heart_3D_cmr_cell{1} = freeROIMask_3D_cmr;
                            mask_blood_3D_cmr_cell{1} = mask_blood_3D_cmr;
                            excludeMask_3D_cmr_cell{1} = excludeMask_3D_cmr;
                            myoRefMask_3D_cmr_cell{1} = myoRefMask_3D_cmr;
                            noReflowMask_3D_cmr_cell{1} = noReflowMask_3D_cmr;

                            vol_img_3D_cmr_cell{cc} = CMR(cc).vol_img_3D.vol_img_3D;
                            mask_myocardium_3D_cmr_cell{cc} = CMR(cc).mask_myocardium_3D.mask_myocardium_3D;
                            freeROIMask_3D_cmr_cell{cc} = CMR(cc).freeROIMask_3D.freeROIMask_3D;
                            mask_heart_3D_cmr_cell{cc} = CMR(cc).mask_heart_3D.mask_heart_3D;
                            mask_blood_3D_cmr_cell{cc} = CMR(cc).mask_blood_3D.mask_blood_3D;
                            excludeMask_3D_cmr_cell{cc} = CMR(cc).excludeMask_3D.excludeMask_3D;
                            myoRefMask_3D_cmr_cell{cc} = CMR(cc).myoRefMask_3D.myoRefMask_3D;
                            noReflowMask_3D_cmr_cell{cc} = CMR(cc).noReflowMask_3D.noReflowMask_3D;
                        else
                            vol_img_3D_cmr = cat(3, vol_img_3D_cmr, CMR(cc).vol_img_3D.vol_img_3D);
                            slc_array_cmr = cat(2, slc_array_cmr, CMR(cc).slc_array.slc_array);
                            mask_myocardium_3D_cmr = cat(3, mask_myocardium_3D_cmr, CMR(cc).mask_myocardium_3D.mask_myocardium_3D);
                            freeROIMask_3D_cmr = cat(3, freeROIMask_3D_cmr, CMR(cc).freeROIMask_3D.freeROIMask_3D);
                            mask_heart_3D_cmr = cat(3, mask_heart_3D_cmr, CMR(cc).mask_heart_3D.mask_heart_3D);
                            mask_blood_3D_cmr = cat(3, mask_blood_3D_cmr, CMR(cc).mask_blood_3D.mask_blood_3D);
                            excludeMask_3D_cmr = cat(3, excludeMask_3D_cmr, CMR(cc).excludeMask_3D.excludeMask_3D);
                            myoRefMask_3D_cmr = cat(3, myoRefMask_3D_cmr, CMR(cc).myoRefMask_3D.myoRefMask_3D);
                            noReflowMask_3D_cmr = cat(3, noReflowMask_3D_cmr, CMR(cc).noReflowMask_3D.noReflowMask_3D);
                        end
                    else

                        if strcmp(label, 'T2_MULTIECHO')

                            freeROIMask_3D_te_cmr = CMR(cc).freeROIMask_3D_te.freeROIMask_3D_te;
                            myoRefMask_3D_te_cmr = CMR(cc).myoRefMask_3D_te.myoRef_3D_te;
                            excludeMask_3D_te_cmr = CMR(cc).excludeMask_3D_te.excludeMask_3D_te;

                            vol_img_3D_cmr_te = CMR(cc).vol_img_3D_te.vol_img_3D_te;
                            Index_cmr = CMR(cc).Index;
                            glob_names_cmr = Index_cmr.glob_names;
                            echo_idx_te_cmr = Index_cmr.echo_idx_te;
                        end
                    end
                end



                %% for cell (CMR EGE analysis)
                if strcmp(label, 'EGE')

                    if size_check == 1 % This only happens in EGE analysis
                        vol_img_3D_cmr_early = vol_img_3D_cmr_cell{1};
                        mask_myocardium_3D_cmr_early = mask_myocardium_3D_cmr_cell{1};
                        freeROIMask_3D_cmr_early = freeROIMask_3D_cmr_cell{1};
                        mask_blood_3D_cmr_early = mask_blood_3D_cmr_cell{1};
                        excludeMask_3D_cmr_early = excludeMask_3D_cmr_cell{1};
                        myoRefMask_3D_cmr_early = myoRefMask_3D_cmr_cell{1};
                        noReflowMask_3D_cmr_early = noReflowMask_3D_cmr_cell{1};

                        vol_img_3D_cmr_pseudo = vol_img_3D_cmr_cell{2};
                        mask_myocardium_3D_cmr_pseudo = mask_myocardium_3D_cmr_cell{2};
                        noReflowMask_3D_cmr_pseudo = noReflowMask_3D_cmr_cell{2};
                        freeROIMask_3D_cmr_pseudo = freeROIMask_3D_cmr_cell{2};
                        mask_blood_3D_cmr_pseudo = mask_blood_3D_cmr_cell{2};
                        excludeMask_3D_cmr_pseudo = excludeMask_3D_cmr_cell{2};
                        myoRefMask_3D_cmr_pseudo = myoRefMask_3D_cmr_cell{2};
                    else
                        vol_img_3D_cmr_early = vol_img_3D_cmr(:,:,1);
                        mask_myocardium_3D_cmr_early = mask_myocardium_3D_cmr(:,:,1);
                        freeROIMask_3D_cmr_early = freeROIMask_3D_cmr(:,:,1);
                        mask_blood_3D_cmr_early = mask_blood_3D_cmr(:,:,1);
                        excludeMask_3D_cmr_early = excludeMask_3D_cmr(:,:,1);
                        myoRefMask_3D_cmr_early = myoRefMask_3D_cmr(:,:,1);
                        noReflowMask_3D_cmr_early = noReflowMask_3D_cmr(:,:,1);

                        vol_img_3D_cmr_pseudo = vol_img_3D_cmr(:,:,2);
                        mask_myocardium_3D_cmr_pseudo = mask_myocardium_3D_cmr(:,:,2);
                        freeROIMask_3D_cmr_pseudo = freeROIMask_3D_cmr(:,:,2);
                        mask_blood_3D_cmr_pseudo = mask_blood_3D_cmr(:,:,2);
                        excludeMask_3D_cmr_pseudo = excludeMask_3D_cmr(:,:,2);
                        myoRefMask_3D_cmr_pseudo = myoRefMask_3D_cmr(:,:,2);
                        noReflowMask_3D_cmr_pseudo = noReflowMask_3D_cmr(:,:,2);
                    end



                    %% TODO
                    x = zeros(size(vol_img_3D_cmr_early, 3),1);
                    y = zeros(size(vol_img_3D_cmr_early, 3),1);

                    InsertionPtCMRPath = cat(2, CMRPath, '/InsertionPts/');
                    if ~exist(InsertionPtCMRPath, 'dir')
                        mkdir(InsertionPtCMRPath);
                    end

                    InsertionPt_CMR_f = cat(2, InsertionPtCMRPath, '/', 'pts_early.mat');

                    if ~exist(InsertionPt_CMR_f)
                        for slc = 1:size(vol_img_3D_cmr_early, 3)
                            figure(105);
                            imagesc(vol_img_3D_cmr_early(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (Early)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 105
                        end
                        ref_pts_cmr_early.x = x;
                        ref_pts_cmr_early.y = y;
                        save(InsertionPt_CMR_f, 'ref_pts_cmr_early');
                    else
                        load(InsertionPt_CMR_f);
                    end

                    clear BaseGroove_cmr_early
                    for slc = 1:size(vol_img_3D_cmr_early,3)
                        C = regionprops(mask_blood_3D_cmr_early(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_cmr_early.x(slc) - x_centroid, ref_pts_cmr_early.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_cmr_early(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_cmr_early(slc) = -Groove_temp -180;
                        else
                            BaseGroove_cmr_early(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end


                    x = zeros(size(vol_img_3D_cmr_pseudo, 3),1);
                    y = zeros(size(vol_img_3D_cmr_pseudo, 3),1);

                    InsertionPtCMRPath = cat(2, CMRPath, '/InsertionPts/');
                    if ~exist(InsertionPtCMRPath, 'dir')
                        mkdir(InsertionPtCMRPath);
                    end

                    InsertionPt_CMR_f = cat(2, InsertionPtCMRPath, '/', 'pts_pseudo.mat');

                    if ~exist(InsertionPt_CMR_f)
                        for slc = 1:size(vol_img_3D_cmr_pseudo, 3)
                            figure(106);
                            imagesc(vol_img_3D_cmr_early(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (Pseudo)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 106
                        end
                        ref_pts_cmr_pseudo.x = x;
                        ref_pts_cmr_pseudo.y = y;
                        save(InsertionPt_CMR_f, 'ref_pts_cmr_pseudo');
                    else
                        load(InsertionPt_CMR_f);
                    end


                    clear BaseGroove_cmr_pseudo
                    for slc = 1:size(vol_img_3D_cmr_pseudo,3)
                        C = regionprops(mask_blood_3D_cmr_pseudo(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_cmr_pseudo.x(slc) - x_centroid, ref_pts_cmr_pseudo.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_cmr_pseudo(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_cmr_pseudo(slc) = -Groove_temp -180;
                        else
                            BaseGroove_cmr_pseudo(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end


                    %% AHA analysis
                    noReflow_early = zeros(size(vol_img_3D_cmr_early));
                    for slc = 1:size(vol_img_3D_cmr_early,3)
                        remote_temp = nonzeros(vol_img_3D_cmr_early(:,:,slc) .* myoRefMask_3D_cmr_early(:,:,slc));
                        thresh = mean(remote_temp) + 5 * std(remote_temp);
                        roi_temp = vol_img_3D_cmr_early(:,:,slc) <= thresh;
                        noReflow_early(:,:,slc) = roi_temp .* mask_myocardium_3D_cmr_early(:,:,slc) .* noReflowMask_3D_cmr_early(:,:,slc);
                    end

                    [Segmentpix_cmr_early, stats_cmr_early, Mask_index_cmr_early] =AHASegmentation_qsm_invivo(noReflow_early, mask_myocardium_3D_cmr_early,6,BaseGroove_cmr_early);
                    t2s_six_segments_mean_cmr_early = squeeze(stats_cmr_early(2,:,:));
                    t2s_six_segments_mean_cmr_early = t2s_six_segments_mean_cmr_early(:);


                    noReflow_pseudo = zeros(size(vol_img_3D_cmr_pseudo));
                    for slc = 1:size(vol_img_3D_cmr_pseudo,3)
                        remote_temp = nonzeros(vol_img_3D_cmr_pseudo(:,:,slc) .* myoRefMask_3D_cmr_pseudo(:,:,slc));
                        thresh = mean(remote_temp) + 5 * std(remote_temp);
                        roi_temp = vol_img_3D_cmr_pseudo(:,:,slc) <= thresh;
                        noReflow_pseudo(:,:,slc) = roi_temp .* mask_myocardium_3D_cmr_pseudo(:,:,slc) .* noReflowMask_3D_cmr_pseudo(:,:,slc);
                    end

                    [Segmentpix_cmr_pseudo, stats_cmr_pseudo, Mask_index_cmr_pseudo] =AHASegmentation_qsm_invivo(noReflow_pseudo, mask_myocardium_3D_cmr_pseudo,6,BaseGroove_cmr_pseudo);
                    t2s_six_segments_mean_cmr_pseudo = squeeze(stats_cmr_pseudo(2,:,:));
                    t2s_six_segments_mean_cmr_pseudo = t2s_six_segments_mean_cmr_pseudo(:);

                elseif strcmp(label, 'LGE') || strcmp(label, 'LGE-MVO')
                    %% Draw insertion point (CMR - LGE)

                    slc_size = size(vol_img_3D_cmr, 3);
                    % idx_mag = 1:2:slc_size;
                    % idx_mag = [5 3 1 20 18 16 14 12 10 8 6];
                    [idx_cmr] = find(any(reshape(freeROIMask_3D_cmr, [], slc_size)));
                    idx_cmr_mag_temp = intersect(idx_mag, idx_cmr);
                    idx_cmr_mag = zeros(1, length(idx_mag));
                    for xx = 1:length(idx_cmr_mag_temp)
                        ind = find(idx_cmr_mag_temp(xx) == idx_mag);
                        idx_cmr_mag(ind) = idx_cmr_mag_temp(xx);
                    end
                    idx_cmr_mag = idx_cmr_mag(~idx_cmr_mag == 0);

                    % idx_cmr_mag = [3 1 20 18 16 14 12 10 8 6];
                    vol_img_3D_cmr = vol_img_3D_cmr(:,:,idx_cmr_mag);
                    slc_array_cmr = slc_array_cmr(idx_cmr_mag);
                    mask_myocardium_3D_cmr = mask_myocardium_3D_cmr(:,:,idx_cmr_mag);
                    freeROIMask_3D_cmr = freeROIMask_3D_cmr(:,:,idx_cmr_mag);
                    mask_heart_3D_cmr = mask_heart_3D_cmr(:,:,idx_cmr_mag);
                    mask_blood_3D_cmr = mask_blood_3D_cmr(:,:,idx_cmr_mag);
                    noReflowMask_3D_cmr = noReflowMask_3D_cmr(:,:,idx_cmr_mag);
                    myoRefMask_3D_cmr = myoRefMask_3D_cmr(:,:,idx_cmr_mag);
                    excludeMask_3D_cmr = excludeMask_3D_cmr(:,:,idx_cmr_mag);


                    x = zeros(size(vol_img_3D_cmr, 3),1);
                    y = zeros(size(vol_img_3D_cmr, 3),1);

                    InsertionPtCMRPath = cat(2, CMRPath, '/InsertionPts/');
                    if ~exist(InsertionPtCMRPath, 'dir')
                        mkdir(InsertionPtCMRPath);
                    end

                    InsertionPt_CMR_f = cat(2, InsertionPtCMRPath, '/', 'pts.mat');

                    if ~exist(InsertionPt_CMR_f)
                        for slc = 1:size(vol_img_3D_cmr, 3)
                            figure(104);
                            imagesc(vol_img_3D_cmr(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (CMR)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 104
                        end
                        ref_pts_cmr.x = x;
                        ref_pts_cmr.y = y;
                        save(InsertionPt_CMR_f, 'ref_pts_cmr');
                    else
                        load(InsertionPt_CMR_f);
                    end

                    clear BaseGroove_cmr
                    for slc = 1:size(vol_img_3D_cmr,3)
                        C = regionprops(mask_blood_3D_cmr(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_cmr.x(slc) - x_centroid, ref_pts_cmr.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_cmr(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_cmr(slc) = -Groove_temp -180;
                        else
                            BaseGroove_cmr(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end

                    [Segmentpix_cmr, stats_cmr, Mask_index_cmr] =AHASegmentation_qsm_invivo(freeROIMask_3D_cmr, mask_myocardium_3D_cmr,6,BaseGroove_cmr);
                    t2s_six_segments_mean_cmr = squeeze(stats_cmr(2,:,:));
                    t2s_six_segments_mean_cmr = t2s_six_segments_mean_cmr(:);


                    noReflow_cmr = zeros(size(vol_img_3D_cmr));
                    for slc = 1:size(vol_img_3D_cmr,3)
                        remote_temp = nonzeros(vol_img_3D_cmr(:,:,slc) .* myoRefMask_3D_cmr(:,:,slc));
                        thresh = mean(remote_temp) + 5 * std(remote_temp);
                        roi_temp = vol_img_3D_cmr(:,:,slc) <= thresh;
                        noReflow_cmr(:,:,slc) = roi_temp .* mask_myocardium_3D_cmr(:,:,slc) .* noReflowMask_3D_cmr(:,:,slc);
                    end

                    [Segmentpix_lgemvo_cmr, stats_lgemvo_cmr, Mask_index_lgemvo_cmr] =AHASegmentation_qsm_invivo(noReflow_cmr, mask_myocardium_3D_cmr,6,BaseGroove_cmr);
                    t2s_six_segments_mean_lgemvo_cmr = squeeze(stats_lgemvo_cmr(2,:,:));
                    t2s_six_segments_mean_lgemvo_cmr = t2s_six_segments_mean_lgemvo_cmr(:);


                elseif strcmp(label, 'T2_MULTIECHO')

                    slc_size = size(vol_img_3D_cmr, 3);
                    [idx_cmr] = find(any(reshape(freeROIMask_3D_te_cmr, [], slc_size)));

                    vol_img_3D_cmr = vol_img_3D_cmr(:,:,idx_cmr);
                    slc_array_cmr = slc_array_cmr(idx_cmr);
                    mask_myocardium_3D_cmr = mask_myocardium_3D_cmr(:,:,idx_cmr);
                    freeROIMask_3D_cmr = freeROIMask_3D_cmr(:,:,idx_cmr);
                    mask_heart_3D_cmr = mask_heart_3D_cmr(:,:,idx_cmr);
                    mask_blood_3D_cmr = mask_blood_3D_cmr(:,:,idx_cmr);
                    excludeMask_3D_cmr = excludeMask_3D_cmr(:,:,idx_cmr);
                    myoRefMask_3D_cmr = myoRefMask_3D_cmr(:,:,idx_cmr);
                    noReflowMask_3D_cmr = noReflowMask_3D_cmr(:,:,idx_cmr);
                    vol_img_3D_cmr_te = vol_img_3D_cmr_te(:,:,idx_cmr);

                    freeROIMask_3D_te_cmr = freeROIMask_3D_te_cmr(:,:,idx_cmr);
                    myoRefMask_3D_te_cmr = myoRefMask_3D_te_cmr(:,:,idx_cmr);
                    excludeMask_3D_te_cmr = excludeMask_3D_te_cmr(:,:,idx_cmr);


                    x = zeros(size(vol_img_3D_cmr_te, 3),1);
                    y = zeros(size(vol_img_3D_cmr_te, 3),1);

                    InsertionPtCMRPath = cat(2, CMRPath, '/InsertionPts/');
                    if ~exist(InsertionPtCMRPath, 'dir')
                        mkdir(InsertionPtCMRPath);
                    end

                    InsertionPt_CMR_f = cat(2, InsertionPtCMRPath, '/', 'pts_te.mat');

                    if ~exist(InsertionPt_CMR_f)
                        for slc = 1:size(vol_img_3D_cmr_te, 3)
                            figure(104);
                            imagesc(vol_img_3D_cmr_te(:,:,slc)); colormap gray; axis image;
                            title('Draw Insertion point (CMR)')
                            [x(slc),y(slc)] = getpts(gca);
                            close figure 104
                        end
                        ref_pts_cmr.x = x;
                        ref_pts_cmr.y = y;
                        save(InsertionPt_CMR_f, 'ref_pts_cmr');
                    else
                        load(InsertionPt_CMR_f);
                    end

                    clear BaseGroove_cmr
                    for slc = 1:size(vol_img_3D_cmr,3)
                        C = regionprops(mask_blood_3D_cmr(:,:,slc));
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);

                        Groove_temp = atan2(ref_pts_cmr.x(slc) - x_centroid, ref_pts_cmr.y(slc) - y_centroid) * 180 / pi;
                        if Groove_temp > 90
                            BaseGroove_cmr(slc) = -Groove_temp + 90;
                        elseif Groove_temp < -90
                            BaseGroove_cmr(slc) = -Groove_temp -180;
                        else
                            BaseGroove_cmr(slc) = Groove_temp;
                        end

                        %disp(ref_pts_cmr.x(slc) - x_centroid)
                        %disp(ref_pts_cmr.y(slc) - y_centroid)
                    end

                    hemo_cmr = zeros(size(mask_myocardium_3D_cmr));
                    thresh_array = zeros(1, size(mask_myocardium_3D_cmr, 3));

                    for slc = 1:size(mask_myocardium_3D_cmr, 3)
                        thresh = mean(nonzeros(myoRefMask_3D_te_cmr(:,:,slc) .* vol_img_3D_cmr_te(:,:,slc))) - 2 * std(nonzeros(myoRefMask_3D_te_cmr(:,:,slc) .* vol_img_3D_cmr_te(:,:,slc)));
                        thresh_array(slc) = thresh;
                        roi_temp = vol_img_3D_cmr_te(:,:,slc) < thresh;
                        hemo_cmr(:,:,slc) = roi_temp .* mask_myocardium_3D_cmr(:,:,slc) .* freeROIMask_3D_te_cmr(:,:,slc) .* ~excludeMask_3D_te_cmr(:,:,slc);
                    end

                    idx_exclude = find(sum(reshape(hemo_cmr, [], size(hemo_cmr,3))) ~= 0);
                    hemo_cmr = hemo_cmr(:,:,idx_exclude);
                    mask_myocardium_3D_cmr = mask_myocardium_3D_cmr(:,:,idx_exclude);
                    BaseGroove_cmr = BaseGroove_cmr(idx_exclude);
                    vol_img_3D_cmr_te = vol_img_3D_cmr_te(:,:,idx_exclude);

                    [Segmentpix_cmr, stats_cmr, Mask_index_cmr] =AHASegmentation_qsm_invivo(hemo_cmr, mask_myocardium_3D_cmr,6,BaseGroove_cmr);
                    t2s_six_segments_mean_cmr = squeeze(stats_cmr(2,:,:));
                    t2s_six_segments_mean_cmr = t2s_six_segments_mean_cmr(:);


                end




                % %% ROC
                % % 1. 5% for Ground Truth
                % posclass = 1;
                % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr>0.05,t2s_six_segments_mean,posclass)
                % figure();
                % plot(X,Y, 'LineWidth', 2);
                % title('Threshold = 0.05')
                %
                % % 2. 1% for Ground Truth
                % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr>0.01,t2s_six_segments_mean,posclass)
                % figure();
                % plot(X,Y, 'LineWidth', 2);
                % title('Threshold = 0.01')
                %
                % % 3. 10% for Ground Truth
                % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr>0.1,t2s_six_segments_mean,posclass)
                % figure();
                % plot(X,Y, 'LineWidth', 2);
                % title('Threshold = 0.10')
                %
                %
                % %%
                % % 2. 1% for Ground Truth
                % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr_early>0.01,t2s_six_segments_mean_lrt_early,posclass)
                % figure();
                % plot(X,Y, 'LineWidth', 2);
                % title('Threshold = 0.01')
                %
                % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr_pseudo>0.01,t2s_six_segments_mean_lrt_early,posclass)
                % figure();
                % plot(X,Y, 'LineWidth', 2);
                % title('Threshold = 0.01')
                % %% 4. 20% for Ground Truth
                % % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr>0.2,t2s_six_segments_mean,posclass)
                % % figure();
                % % plot(X,Y, 'LineWidth', 2);
                % % title('Threshold = 0.20')
                % %
                % % % 5. 30% for Ground Truth
                % % [X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr>0.3,t2s_six_segments_mean,posclass)
                % % figure();
                % % plot(X,Y, 'LineWidth', 2);
                % % title('Threshold = 0.30')

                % save as struct
                if ~(strcmp(label, 'EGE') || strcmp(label, 'LGE-MVO'))
                    AHA_segs(ll).t2s_six_segments_mean_lrt = t2s_six_segments_mean_lrt;
                    AHA_segs(ll).t2s_six_segments_mean_cmr = t2s_six_segments_mean_cmr;

                    AHA_segs(ll).Segmentpix_cmr = Segmentpix_cmr;
                    AHA_segs(ll).stats_cmr = stats_cmr;
                    AHA_segs(ll).Mask_index_cmr = Mask_index_cmr;
                    AHA_segs(ll).Segmentpix_lrt = Segmentpix_lrt;
                    AHA_segs(ll).stats_lrt = stats_lrt;
                    AHA_segs(ll).Mask_index_lrt = Mask_index_lrt;

                elseif strcmp(label, 'EGE')
                    AHA_segs(ll).t2s_six_segments_mean_lrt_pseudo = t2s_six_segments_mean_lrt_pseudo;
                    AHA_segs(ll).t2s_six_segments_mean_lrt_early = t2s_six_segments_mean_lrt_early;
                    AHA_segs(ll).t2s_six_segments_mean_cmr_pseudo = t2s_six_segments_mean_cmr_pseudo;
                    AHA_segs(ll).t2s_six_segments_mean_cmr_early = t2s_six_segments_mean_cmr_early;

                    AHA_segs(ll).Segmentpix_cmr_pseudo = Segmentpix_cmr_pseudo;
                    AHA_segs(ll).stats_cmr_pseudo = stats_cmr_pseudo;
                    AHA_segs(ll).Mask_index_cmr_pseudo = Mask_index_cmr_pseudo;
                    AHA_segs(ll).Segmentpix_lrt_pseudo = Segmentpix_lrt_pseudo;
                    AHA_segs(ll).stats_lrt_pseudo = stats_lrt_pseudo;
                    AHA_segs(ll).Mask_index_lrt_pseudo = Mask_index_lrt_pseudo;

                    AHA_segs(ll).Segmentpix_cmr_early = Segmentpix_cmr_early;
                    AHA_segs(ll).stats_cmr_early = stats_cmr_early;
                    AHA_segs(ll).Mask_index_cmr_early = Mask_index_cmr_early;
                    AHA_segs(ll).Segmentpix_lrt_early = Segmentpix_lrt_early;
                    AHA_segs(ll).stats_lrt_early = stats_lrt_early;
                    AHA_segs(ll).Mask_index_lrt_early = Mask_index_lrt_early;

                elseif strcmp(label, 'LGE-MVO')
                    AHA_segs(ll).t2s_six_segments_mean_lgemvo_lrt = t2s_six_segments_mean_lgemvo_lrt;
                    AHA_segs(ll).t2s_six_segments_mean_lgemvo_cmr = t2s_six_segments_mean_lgemvo_cmr;

                    AHA_segs(ll).Segmentpix_lgemvo_cmr = Segmentpix_lgemvo_cmr;
                    AHA_segs(ll).stats_lgemvo_cmr = stats_lgemvo_cmr;
                    AHA_segs(ll).Mask_index_lgemvo_cmr = Mask_index_lgemvo_cmr;
                    AHA_segs(ll).Segmentpix_lgemvo_lrt = Segmentpix_lgemvo_lrt;
                    AHA_segs(ll).stats_lgemvo_lrt = stats_lgemvo_lrt;
                    AHA_segs(ll).Mask_index_lgemvo_lrt = Mask_index_lgemvo_lrt;
                end



            end


            save_path = cat(2, TimePtPath, 'AHA_segs.mat');
            if length(AHA_segs) > (length(sequence_label) - 1)
                save(save_path, 'AHA_segs');
            end

        end

    end
end















%% Image Display
figure();
for slc = 1:size(vol_img_3D_lrt, 3)
    subplot(2,3,slc)
    imagesc(vol_img_3D_lrt(:,:,slc).*mask_myocardium_3D_lrt(:,:,slc)); colormap gray; axis image;
end

figure();
for slc = 1:size(vol_img_3D_cmr_te, 3)
    subplot(2,3,slc)
    imagesc(vol_img_3D_cmr_te(:,:,slc).*mask_myocardium_3D_cmr(:,:,slc)); colormap gray; axis image;
end

%%
figure();
for slc = 1:size(Mask_index_lrt, 3)
    subplot(2,3,slc)
    imagesc(Mask_index_lrt(:,:,slc)); colormap gray; axis image;
end

figure();
for slc = 1:size(Mask_index_cmr, 3)
    subplot(2,3,slc)
    imagesc(Mask_index_cmr(:,:,slc)); colormap gray; axis image;
end

%%
%temp = circshift(vol_img_3D_lrt,5,3);
figure();
for slc = 1:size(vol_img_3D_cmr_te, 3)
    subplot(4,4,slc)

    imagesc(vol_img_3D_cmr_te(:,:,slc)); colormap gray; axis image;

    subplot(4,4,slc+8)
    if slc <= size(vol_img_3D_lrt, 3)
        imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
    end
end


%% ROC
% 1. 5% for Ground Truth
posclass = 1;
[X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr_pseudo>0.05,t2s_six_segments_mean_lrt_early,posclass)
figure();
plot(X,Y, 'LineWidth', 2);
title('Threshold = 0.05')

% 2. 1% for Ground Truth
[X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr_pseudo>0.01,t2s_six_segments_mean_lrt_early,posclass)
figure();
plot(X,Y, 'LineWidth', 2);
title('Threshold = 0.01')

% 3. 10% for Ground Truth
[X,Y,T,AUC,OPTROCPT] = perfcurve(t2s_six_segments_mean_cmr_pseudo>0.1,t2s_six_segments_mean_lrt_early,posclass)
figure();
plot(X,Y, 'LineWidth', 2);
title('Threshold = 0.10')
%%
figure();
for slc = 1:size(vol_img_3D_cmr_te, 3)
    subplot(4,4,slc)
    imagesc(hemo_cmr(:,:,slc).*mask_myocardium_3D_cmr(:,:,slc)); colormap gray; axis image;

end

%%
shift_idx = [4 3 2 1 10 9 8 7 6 5];
temp = vol_img_3D_lrt(:,:,shift_idx);
figure();
for slc = 1:size(vol_img_3D_lrt, 3)
    subplot(4,4,slc)
    imagesc(temp(:,:,slc)); colormap gray; axis image;

end

%%
figure();
for slc = 1:size(vol_img_3D, 3)
    subplot(4,4,slc)
    imagesc(vol_img_3D(:,:,slc)); colormap gray; axis image;
end

%%
figure();
for slc = 1:size(vol_img_3D_lrt, 3)
    subplot(4,4,slc)
    imagesc(excludeMask_3D(:,:,slc)); colormap gray; axis image;
end
%%
% mask_myocardium_3D_lrt
% freeROIMask_3D_lrt
% mask_heart_3D_lrt
% mask_blood_3D_lrt
% noReflowMask_3D_lrt;
% excludeMask_3D_lrt;
% myoRefMask_3D_lrt;
