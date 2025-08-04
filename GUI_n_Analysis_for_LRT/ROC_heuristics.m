addpath('../function/');
addpath('../GUI_n_Analysis_for_LRT/');
addpath('../T1NFF/');

% ================================ Identify your major folder, in this case
% ROC_analysis/
base_dir = uigetdir; % more generic -> % ROC_analysis
% base_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
folder_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = ExtractNames(folder_glob);

time_points = {'D6', 'D8'};
time_points = {'D5', 'D7'};
time_points = {'WK8', 'WK8+2'}
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
Heuristics = struct;

for n = starting_point:starting_point
    name = Names{n};
    real_name_temp = strsplit(name, '_');
    real_name = real_name_temp{1};


    SubjectPath = GetFullPath(cat(2, OutputPath, '/', name, '/'));
    if ~exist(SubjectPath, 'dir')
        mkdir(SubjectPath);
    end

    Heuristics.(name) = struct;
    for tp = 1:length(time_points)
        %for tp = 2:2

        time_point = time_points{end-tp+1};

        contour_dir_lrt_check = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{1}, '/', sequence_label{1}, '/');
        if ~exist(cat(2, contour_dir_lrt_check, sequence_label{1}, '_vol_img_3D.mat'))

            disp(['skip: ', time_point])

        else

            Heuristics.(name).time_point = time_point;


            for ll = 1:length(sequence_label)
                % for ll = 2:2
                label = sequence_label{ll};

                if (strcmp(time_point, 'WK8') || strcmp(time_point, 'WK8+2')) && (strcmp(label, 'EGE') || strcmp(label, 'LGE-MVO'))
                    

                    continue;

                else
                    contour_dir_lrt = cat(2, base_dir, '/ContourData/', name, '/', time_point, '/', methods{1}, '/', label, '/');
                    LRT = load(cat(2, contour_dir_lrt, label, '_vol_img_3D.mat'));
                    myoRef_LRT = load(cat(2, contour_dir_lrt, '/MyoReference/myoRef.mat'));
                    noReflow_LRT = load(cat(2, contour_dir_lrt, '/noReflowArea/noReflow.mat'));
                    freeroi_LRT = load(cat(2, contour_dir_lrt, '/freeROI/freeROI.mat'));

                    vol_img_3D_lrt = LRT.vol_img_3D;
                    freeROIMask_3D_lrt = freeroi_LRT.freeROIMask_3D;
                    noReflowMask_3D_lrt = noReflow_LRT.noReflowMask_3D;
                    myoRefMask_3D_lrt = myoRef_LRT.myoRefMask_3D;

                    slc_size = size(vol_img_3D_lrt, 3);
                    if strcmp(label, 'EGE')
                        [idx_lrt] = find(any(reshape(noReflowMask_3D_lrt, [], slc_size)));
                    elseif strcmp(label, 'LGE')
                        [idx_lrt] = find(any(reshape(freeROIMask_3D_lrt, [], slc_size)));
                    elseif strcmp(label, 'T2_MULTIECHO')
                        [idx_lrt] = find(any(reshape(myoRefMask_3D_lrt, [], slc_size)));
                    end
                    vol_img_3D_lrt = vol_img_3D_lrt(:,:,idx_lrt);



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
                    else
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
                    end

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



                    if strcmp(label, 'LGE')

                        slc_size_lrt = size(vol_img_3D_lrt, 3);
                        n = ceil(sqrt(slc_size_lrt));
                        figure();
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            subplot(n,n,slc)
                            imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
                        end
                        disp(' Please determine the shift_idx. apex -> base ');
                        pause;

                        slc_size_cmr = size(vol_img_3D_cmr, 3);
                        n = ceil(sqrt(slc_size_cmr));
                        figure();
                        for slc = 1:size(vol_img_3D_cmr, 3)
                            subplot(n,n,slc)
                            imagesc(vol_img_3D_cmr(:,:,slc)); colormap gray; axis image;
                        end
                        disp(' Please determine the idx_mag. apex -> base ');
                        pause;


                    elseif strcmp(label, 'EGE')
                        slc_size_lrt = size(vol_img_3D_lrt, 3);
                        n = ceil(sqrt(slc_size_lrt));
                        figure();
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            subplot(n,n,slc)
                            imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
                        end
                        disp(' Make sure this is from early to late enhancement (LRT) ');
                        pause;



                    elseif strcmp(label, 'T2_MULTIECHO')
                        slc_size_lrt_te = size(vol_img_3D_lrt, 3);
                        n = ceil(sqrt(slc_size_lrt_te));
                        figure();
                        for slc = 1:size(vol_img_3D_lrt, 3)
                            subplot(n,n,slc)
                            imagesc(vol_img_3D_lrt(:,:,slc)); colormap gray; axis image;
                        end
                        disp(' Please determine the matching_cell [shift_idx (shift toward right), reverse_label]. apex -> base');
                        pause;


                        slc_size_cmr_te = size(vol_img_3D_cmr, 3);
                        n = ceil(sqrt(slc_size_cmr_te));
                        figure();
                        for slc = 1:size(vol_img_3D_cmr, 3)
                            subplot(n,n,slc)
                            imagesc(vol_img_3D_cmr(:,:,slc)); colormap gray; axis image;
                        end
                        disp(' Make sure this is Apex to Base ');

                    end
                end
            end
        end
    end
end

%%
% LISBON_acute
shift_idx = [5 4 3 2 1 6];
idx_mag = 1:2:slc_size_cmr;
matching_cell = [5, 0];

% CARLOS_acute
shift_idx = [5 6 7 8 9 10 1 2 3 4];
idx_mag = [6 8 10 12 14 16 18 20 1 3 5];
matching_cell = [5, 0];

% CHILI_acute
shift_idx = [2 3 4 5 6 1];
idx_mag = [3 5 7 9 11 13 15 1];
matching_cell = [5, 0];

% CINNAMON_acute
shift_idx = [2 3 4 5 1];
idx_mag = 1:2:slc_size_cmr;
matching_cell = [];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% DAVE_acute
shift_idx = [1 2 3 4 5];
idx_mag = [2 4 6 8 10 12 14 1];
matching_cell = [0, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% GINER_acute
shift_idx = [1 2 3 4 5];
idx_mag = [3 5 7 9 11 13 15 1];
matching_cell = [4, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% NUTMEG_acute
shift_idx = [2 3 4 5 6 1];
idx_mag = 1:2:slc_size_cmr;
matching_cell = [6, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% PAPRIKA_acute
shift_idx = [1 2 3 4 5];
idx_mag = [3 5 7 9 11 13 15 1];
matching_cell = [0, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% PARIS_acute
shift_idx = [1 2 3 4 5 6];
idx_mag = [3 5 7 9 11 13 15 17 1];
matching_cell = [];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% JESSE_acute
shift_idx = [5 6 7 8 1 2 3 4];
idx_mag = [3 5 7 9 11 13 15 1];
matching_cell = [5, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% SOFIA_acute
shift_idx = [5 4 3 2 1];
idx_mag = [3 5 7 9 11 13 15 17 1];
matching_cell = [];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% LISBON_chronic
shift_idx = [3 4 5 6 7 1 2];
idx_mag = [3 5 7 9 11 13 15 1];
matching_cell = [6, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% SOFIA_chronic
shift_idx = [4 3 2 1];
idx_mag = [5 7 9 11 13 15 17 19 1 3];
matching_cell = [];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% PARIS_chronic
shift_idx = [1 2 3 4 5];
idx_mag = [5 7 9 11 13 15 17 19 1 3];
matching_cell = [];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

% JESSE_chronic
shift_idx = [2 3 4 5 6 7 1];
idx_mag = [1 3 5 7 9 11 13];
matching_cell = [5, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

%% 
% LATTE_chronic
shift_idx = [1 2 3 4];
idx_mag = [2 3 4 5 6 7 8 1];
matching_cell = [0, 0];
Heuristics.(name).shift_idx = shift_idx;
Heuristics.(name).idx_mag = idx_mag;
Heuristics.(name).matching_cell = matching_cell;

save_path = cat(2, SubjectPath, 'Heuristics.mat');
if ~exist(save_path, 'file')
    save(save_path, 'Heuristics');
end





