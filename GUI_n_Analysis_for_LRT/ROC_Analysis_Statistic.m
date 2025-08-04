clear all;
close all;

addpath('../function/');
addpath('../GUI_n_Analysis_for_LRT/');
addpath('../T1NFF/');

base_dir = uigetdir; % more generic -> % ROC_analysis
folder_glob = glob(cat(2, base_dir, '/Results/*'));
OutputPath = GetFullPath(cat(2, base_dir, '/Results/'));
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end

Names = ExtractNames(folder_glob);

time_points = {'D5', 'D6', 'D7', 'D8', 'WK8', 'WK8+2'};
% time_points_chronic = {'WK8', 'WK8+2'};
%%
%name_check = 'SOFIA_acute';
%starting_point = find(strcmp(name_check, Names),1);

AHA_segs_lrt_LGE_array = [];
AHA_segs_cmr_LGE_array = [];
AHA_segs_lrt_EGE1_array = [];
AHA_segs_cmr_EGE1_array = [];
AHA_segs_lrt_EGE2_array = [];
AHA_segs_cmr_EGE2_array = [];
AHA_segs_lrt_T2MULTIECHO_array = [];
AHA_segs_cmr_T2MULTIECHO_array = [];
AHA_segs_lrt_LGEMVO_array = [];
AHA_segs_cmr_LGEMVO_array = [];

AHA_segs_lrt_LGE_array_acute = [];
AHA_segs_cmr_LGE_array_acute = [];
AHA_segs_lrt_EGE1_array_acute = [];
AHA_segs_cmr_EGE1_array_acute = [];
AHA_segs_lrt_EGE2_array_acute = [];
AHA_segs_cmr_EGE2_array_acute = [];
AHA_segs_lrt_T2MULTIECHO_array_acute = [];
AHA_segs_cmr_T2MULTIECHO_array_acute = [];
AHA_segs_lrt_LGEMVO_array_acute = [];
AHA_segs_cmr_LGEMVO_array_acute = [];

AHA_segs_lrt_LGE_array_chronic = [];
AHA_segs_cmr_LGE_array_chronic = [];
AHA_segs_lrt_EGE1_array_chronic = [];
AHA_segs_cmr_EGE1_array_chronic = [];
AHA_segs_lrt_EGE2_array_chronic = [];
AHA_segs_cmr_EGE2_array_chronic = [];
AHA_segs_lrt_T2MULTIECHO_array_chronic = [];
AHA_segs_cmr_T2MULTIECHO_array_chronic = [];
AHA_segs_lrt_LGEMVO_array_chronic = [];
AHA_segs_cmr_LGEMVO_array_chronic = [];

for n = 1:length(Names)
% for n = 1:9
    
    name = Names{n};
    disp(name)
    real_name_temp = strsplit(name, '_');
    real_name = real_name_temp{1};

    SubjectPath = GetFullPath(cat(2, OutputPath, '/', name, '/'));
    if ~exist(SubjectPath, 'dir')
        mkdir(SubjectPath);
    end

    for tp = 1:length(time_points)
        time_point = time_points{end-tp+1};
        AHA_segs_check = cat(2, SubjectPath, time_point, '/AHA_segs.mat');
        if ~exist(AHA_segs_check)

            disp(['skip: ', time_point])

        else
            load(AHA_segs_check);
            AHA_segs_lrt_LGE = AHA_segs(1).t2s_six_segments_mean_lrt;
            AHA_segs_cmr_LGE = AHA_segs(1).t2s_six_segments_mean_cmr;


            if strcmp(name, 'JESSE_acute')
                AHA_segs_lrt_LGE(1:6) = circshift(AHA_segs_lrt_LGE(1:6), 5);
                AHA_segs_lrt_LGE(7:12) = circshift(AHA_segs_lrt_LGE(7:12), 1);
                AHA_segs_lrt_LGE(13:18) = circshift(AHA_segs_lrt_LGE(13:18), 1);
                AHA_segs_lrt_LGE(19:24) = circshift(AHA_segs_lrt_LGE(19:24), 1);
                AHA_segs_cmr_LGE(37) = 0;
            elseif strcmp(name, 'LISBON_chronic')
                AHA_segs_lrt_LGE(25:30) = [];
            elseif strcmp(name, 'SOFIA_chronic')
                AHA_segs_lrt_LGE(1:6) = circshift(AHA_segs_lrt_LGE(1:6), 1);
            elseif strcmp(name, 'JESSE_chronic')
                AHA_segs_lrt_LGE(1:6) = circshift(AHA_segs_lrt_LGE(1:6), 2);
                AHA_segs_lrt_LGE(10) = 0;
            elseif strcmp(name, 'LATTE_chronic')
                AHA_segs_lrt_LGE(1:6) = circshift(AHA_segs_lrt_LGE(1:6), 1);
                AHA_segs_lrt_LGE(7:12) = circshift(AHA_segs_lrt_LGE(7:12), 1);
                % AHA_segs_lrt_LGE(3) = 0;
            end


            AHA_segs_lrt_EGE1 = AHA_segs(2).t2s_six_segments_mean_lrt_early;
            AHA_segs_lrt_EGE2 = AHA_segs(2).t2s_six_segments_mean_lrt_pseudo;
            
            if strcmp(name, 'CARLOS_acute')
                AHA_segs_lrt_EGE1 = circshift(AHA_segs(2).t2s_six_segments_mean_lrt_early, 1);
            elseif strcmp(name, 'PAPRIKA_acute')
                AHA_segs_lrt_EGE1 = circshift(AHA_segs(2).t2s_six_segments_mean_lrt_early, 5);
            end

            if strcmp(name, 'SOFIA_acute')
                AHA_segs_lrt_EGE2(1:6) = [0;0;0;0;0;0];
            elseif strcmp(name, 'CARLOS_acute')
                AHA_segs_lrt_EGE2 = circshift(AHA_segs(2).t2s_six_segments_mean_lrt_pseudo, 1);
            elseif strcmp(name, 'GINGER_acute')
                AHA_segs_lrt_EGE2 = circshift(AHA_segs(2).t2s_six_segments_mean_lrt_pseudo, 5);
            elseif strcmp(name, 'PARIS_acute')
                AHA_segs_lrt_EGE2(1) = 0;
            end

            AHA_segs_cmr_EGE1 = AHA_segs(2).t2s_six_segments_mean_cmr_early;
            AHA_segs_cmr_EGE2 = AHA_segs(2).t2s_six_segments_mean_cmr_pseudo;
            AHA_segs_lrt_T2MULTIECHO = AHA_segs(3).t2s_six_segments_mean_lrt;
            AHA_segs_cmr_T2MULTIECHO = AHA_segs(3).t2s_six_segments_mean_cmr;

            
            if strcmp(name, 'PAPRIKA_acute')
                AHA_segs_lrt_T2MULTIECHO(1:6) = circshift(AHA_segs_lrt_T2MULTIECHO(1:6), 1);
                AHA_segs_cmr_T2MULTIECHO(12) = 0;
            elseif strcmp(name, 'JESSE_acute')
                AHA_segs_lrt_T2MULTIECHO(1:6) = circshift(AHA_segs_lrt_T2MULTIECHO(1:6), 5);
            elseif strcmp(name, 'CARLOS_acute')
                AHA_segs_lrt_T2MULTIECHO(7:12) = circshift(AHA_segs_lrt_T2MULTIECHO(7:12), 5);
                AHA_segs_cmr_T2MULTIECHO(27) = 0;
            elseif strcmp(name, 'LISBON_chronic')
                AHA_segs_lrt_T2MULTIECHO(31) = AHA_segs_lrt_T2MULTIECHO(32);
                AHA_segs_lrt_T2MULTIECHO(32) = 0;
                AHA_segs_lrt_T2MULTIECHO(19:24) = [];
            elseif strcmp(name, 'JESSE_chronic')
                AHA_segs_lrt_T2MULTIECHO(1:6) = circshift(AHA_segs_lrt_T2MULTIECHO(1:6), 2);
            elseif strcmp(name, 'LATTE_chronic')
                AHA_segs_cmr_T2MULTIECHO(16) = 0;
            end



            if length(AHA_segs) == 4
                AHA_segs_lrt_LGEMVO = AHA_segs(4).t2s_six_segments_mean_lgemvo_lrt;
                AHA_segs_cmr_LGEMVO = AHA_segs(4).t2s_six_segments_mean_lgemvo_cmr;

                if strcmp(name, 'JESSE_acute')
                    AHA_segs_lrt_LGEMVO(1:6) = circshift(AHA_segs_lrt_LGEMVO(1:6), 5);
                    AHA_segs_lrt_LGEMVO(7:12) = circshift(AHA_segs_lrt_LGEMVO(7:12), 1);
                    AHA_segs_lrt_LGEMVO(13:18) = circshift(AHA_segs_lrt_LGEMVO(13:18), 1);
                    AHA_segs_lrt_LGEMVO(19:24) = circshift(AHA_segs_lrt_LGEMVO(19:24), 1);
                    AHA_segs_cmr_LGE(25) = 0;
                    AHA_segs_cmr_LGE(31) = 0;
                else

                end

                
                AHA_segs_lrt_LGE_array = [AHA_segs_lrt_LGE_array; AHA_segs_lrt_LGE];
                AHA_segs_cmr_LGE_array = [AHA_segs_cmr_LGE_array; AHA_segs_cmr_LGE];

                AHA_segs_lrt_EGE1_array = [AHA_segs_lrt_EGE1_array; AHA_segs_lrt_EGE1];
                AHA_segs_cmr_EGE1_array = [AHA_segs_cmr_EGE1_array; AHA_segs_cmr_EGE1];
                AHA_segs_lrt_EGE2_array = [AHA_segs_lrt_EGE2_array; AHA_segs_lrt_EGE2];
                AHA_segs_cmr_EGE2_array = [AHA_segs_cmr_EGE2_array; AHA_segs_cmr_EGE2];

                AHA_segs_lrt_T2MULTIECHO_array = [AHA_segs_lrt_T2MULTIECHO_array; AHA_segs_lrt_T2MULTIECHO];
                AHA_segs_cmr_T2MULTIECHO_array = [AHA_segs_cmr_T2MULTIECHO_array; AHA_segs_cmr_T2MULTIECHO];

                AHA_segs_lrt_LGEMVO_array = [AHA_segs_lrt_LGEMVO_array; AHA_segs_lrt_LGEMVO];
                AHA_segs_cmr_LGEMVO_array = [AHA_segs_cmr_LGEMVO_array; AHA_segs_cmr_LGEMVO];

                if tp > 2
                    AHA_segs_lrt_LGE_array_acute = [AHA_segs_lrt_LGE_array_acute; AHA_segs_lrt_LGE];
                    AHA_segs_cmr_LGE_array_acute = [AHA_segs_cmr_LGE_array_acute; AHA_segs_cmr_LGE];

                    AHA_segs_lrt_EGE1_array_acute = [AHA_segs_lrt_EGE1_array_acute; AHA_segs_lrt_EGE1];
                    AHA_segs_cmr_EGE1_array_acute = [AHA_segs_cmr_EGE1_array_acute; AHA_segs_cmr_EGE1];
                    AHA_segs_lrt_EGE2_array_acute = [AHA_segs_lrt_EGE2_array_acute; AHA_segs_lrt_EGE2];
                    AHA_segs_cmr_EGE2_array_acute = [AHA_segs_cmr_EGE2_array_acute; AHA_segs_cmr_EGE2];

                    AHA_segs_lrt_T2MULTIECHO_array_acute = [AHA_segs_lrt_T2MULTIECHO_array_acute; AHA_segs_lrt_T2MULTIECHO];
                    AHA_segs_cmr_T2MULTIECHO_array_acute = [AHA_segs_cmr_T2MULTIECHO_array_acute; AHA_segs_cmr_T2MULTIECHO];

                    AHA_segs_lrt_LGEMVO_array_acute = [AHA_segs_lrt_LGEMVO_array_acute; AHA_segs_lrt_LGEMVO];
                    AHA_segs_cmr_LGEMVO_array_acute = [AHA_segs_cmr_LGEMVO_array_acute; AHA_segs_cmr_LGEMVO];

                elseif tp <= 2
                    AHA_segs_lrt_LGE_array_chronic = [AHA_segs_lrt_LGE_array_chronic; AHA_segs_lrt_LGE];
                    AHA_segs_cmr_LGE_array_chronic = [AHA_segs_cmr_LGE_array_chronic; AHA_segs_cmr_LGE];

                    AHA_segs_lrt_T2MULTIECHO_array_chronic = [AHA_segs_lrt_T2MULTIECHO_array_chronic; AHA_segs_lrt_T2MULTIECHO];
                    AHA_segs_cmr_T2MULTIECHO_array_chronic = [AHA_segs_cmr_T2MULTIECHO_array_chronic; AHA_segs_cmr_T2MULTIECHO];

                end
            end

        end
    end
end

%% ROC
filename = 'ROC_Total'
thresh_array = [0.01, 0.05, 0.10];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[0 600 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array>thresh,AHA_segs_lrt_LGE_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end
saveas(gcf, fullfile(OutputPath, [filename, '_MI.png']))

figure('Position',[0 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_EGE1_array>thresh,AHA_segs_lrt_EGE1_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end
saveas(gcf, fullfile(OutputPath, [filename, '_MVO1.png']))


figure('Position',[500 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_EGE2_array>thresh,AHA_segs_lrt_EGE2_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

saveas(gcf, fullfile(OutputPath, [filename, '_MVO2.png']))

figure('Position',[0 200 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array>thresh,AHA_segs_lrt_T2MULTIECHO_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

saveas(gcf, fullfile(OutputPath, [filename, '_MVO3.png']))



figure('Position',[1000 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGEMVO_array>thresh,AHA_segs_lrt_LGEMVO_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end
saveas(gcf, fullfile(OutputPath, [filename, '_HEMO.png']))

%% Acute
thresh_array = [0.01, 0.05, 0.10];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[0 600 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array_acute>thresh,AHA_segs_lrt_LGE_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

figure('Position',[0 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_EGE1_array_acute>thresh,AHA_segs_lrt_EGE1_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

figure('Position',[500 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_EGE2_array_acute>thresh,AHA_segs_lrt_EGE2_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

figure('Position',[0 200 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_acute>thresh,AHA_segs_lrt_T2MULTIECHO_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end


figure('Position',[1000 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGEMVO_array_acute>thresh,AHA_segs_lrt_LGEMVO_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

%% Chronic
thresh_array = [0.01, 0.05, 0.10];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[0 600 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array_chronic>thresh,AHA_segs_lrt_LGE_array_chronic,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end



figure('Position',[0 400 500 125]);
for i = 1:length(thresh_array)
    subplot(1,3,i);
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_chronic>thresh,AHA_segs_lrt_T2MULTIECHO_array_chronic,posclass);
    plot(X, Y, 'LineWidth', 2);
    title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    text(0.3, 0.5, txt);
end

%% HOLD ON
thresh_array = [0.01, 0.05, 0.10];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[0 600 500 125]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    subplot(1,3,i);
    hold on;
    
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array>thresh,AHA_segs_lrt_LGE_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array_acute>thresh,AHA_segs_lrt_LGE_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array_chronic>thresh,AHA_segs_lrt_LGE_array_chronic,posclass);
    plot(X, Y, 'LineWidth', 2);

    title(cat(2, 'Threshold = ', num2str(thresh)));
    %txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    %text(0.3, 0.5, txt);
end

%% HOLD ON
thresh_array = [0.01, 0.05, 0.10];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[0 600 500 125]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    subplot(1,3,i);
    hold on;
    
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array>thresh,AHA_segs_lrt_T2MULTIECHO_array,posclass);
    plot(X, Y, 'LineWidth', 2);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_acute>thresh,AHA_segs_lrt_T2MULTIECHO_array_acute,posclass);
    plot(X, Y, 'LineWidth', 2);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_chronic>thresh,AHA_segs_lrt_T2MULTIECHO_array_chronic,posclass);
    plot(X, Y, 'LineWidth', 2);

    title(cat(2, 'Threshold = ', num2str(thresh)));
    %txt = cat(2, 'AUC = ', num2str(AUC, formatSpec));
    %text(0.3, 0.5, txt);
end

%% Produce quality plots
filename = 'ROC_Total'
thresh_array = [0.01];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[100 100 650 600]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array>thresh,AHA_segs_lrt_LGE_array,posclass);
    plot(X, Y, '-.', 'LineWidth', 4, 'Color', 'k');
    se = getAUC_CI(AHA_segs_cmr_LGE_array, AHA_segs_lrt_LGE_array, AUC);
    hold on;
    [X_acute,Y_acute,T_acute,AUC_acute,OPTROCPT_acute] = perfcurve(AHA_segs_cmr_LGE_array_acute>thresh,AHA_segs_lrt_LGE_array_acute,posclass);
    plot(X_acute, Y_acute, 'LineWidth', 4, 'Color', [179 0 0]/255);
    se_acute = getAUC_CI(AHA_segs_cmr_LGE_array_acute, AHA_segs_lrt_LGE_array_acute, AUC_acute);

    [X_chronic,Y_chronic,T_chronic,AUC_chronic,OPTROCPT_chronic] = perfcurve(AHA_segs_cmr_LGE_array_chronic>thresh,AHA_segs_lrt_LGE_array_chronic,posclass);
    plot(X_chronic, Y_chronic, 'LineWidth', 4, 'Color', [4 90 141]/255);
    se_chronic = getAUC_CI(AHA_segs_cmr_LGE_array_chronic, AHA_segs_lrt_LGE_array_chronic, AUC_chronic);

    %title(cat(2, 'Threshold = ', num2str(thresh)));

    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec), ', SE = ', num2str(1.96*se, formatSpec));
    text(0.3, 0.5, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_acute, formatSpec), ', SE = ', num2str(1.96*se_acute, formatSpec));
    text(0.3, 0.4, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_chronic, formatSpec), ', SE = ', num2str(1.96*se_chronic, formatSpec));
    text(0.3, 0.3, txt);
end
%saveas(gcf, fullfile(OutputPath, [filename, '_MI.png']))


%% HEMO
filename = 'ROC_Total'
thresh_array = [0.01];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[100 100 650 600]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array>thresh,AHA_segs_lrt_T2MULTIECHO_array,posclass);
    plot(X, Y, '-.', 'LineWidth', 4, 'Color', 'k');
    se = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array, AHA_segs_lrt_T2MULTIECHO_array, AUC);

    hold on;
    [X_acute,Y_acute,T_acute,AUC_acute,OPTROCPT_acute] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_acute>thresh,AHA_segs_lrt_T2MULTIECHO_array_acute,posclass);
    plot(X_acute, Y_acute, 'LineWidth', 4, 'Color', [179 0 0]/255);
    se_acute = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array_acute, AHA_segs_lrt_T2MULTIECHO_array_acute, AUC_acute);

    [X_chronic,Y_chronic,T_chronic,AUC_chronic,OPTROCPT_chronic] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_chronic>thresh,AHA_segs_lrt_T2MULTIECHO_array_chronic,posclass);
    plot(X_chronic, Y_chronic, 'LineWidth', 4, 'Color', [4 90 141]/255);
    se_chronic = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array_chronic, AHA_segs_lrt_T2MULTIECHO_array_chronic, AUC_chronic);

    %title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec), ', SE = ', num2str(1.96*se, formatSpec));
    text(0.3, 0.5, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_acute, formatSpec), ', SE = ', num2str(1.96*se_acute, formatSpec));
    text(0.3, 0.4, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_chronic, formatSpec), ', SE = ', num2str(1.96*se_chronic, formatSpec));
    text(0.3, 0.3, txt);
end

%% MVO
filename = 'ROC_MVO'
thresh_array = [0.01];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[100 100 650 600]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    % [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGEMVO_array>thresh,AHA_segs_lrt_LGEMVO_array,posclass);
    % plot(X, Y, '-.', 'LineWidth', 6, 'Color', 'k');
    % se = getAUC_CI(AHA_segs_cmr_LGEMVO_array, AHA_segs_lrt_LGEMVO_array, AUC);

    hold on;
    [X_acute,Y_acute,T_acute,AUC_acute,OPTROCPT_acute] = perfcurve(AHA_segs_cmr_EGE1_array_acute>thresh,AHA_segs_lrt_EGE1_array_acute,posclass);
    plot(X_acute, Y_acute, 'LineWidth', 6, 'Color', [179 0 0]/255);
    se_ege1 = getAUC_CI(AHA_segs_cmr_EGE1_array_acute, AHA_segs_lrt_EGE1_array_acute, AUC_acute);

    [X_chronic,Y_chronic,T_chronic,AUC_chronic,OPTROCPT_chronic] = perfcurve(AHA_segs_cmr_EGE2_array_acute>thresh,AHA_segs_lrt_EGE2_array_acute,posclass);
    plot(X_chronic, Y_chronic, 'LineWidth', 6, 'Color', [4 90 141]/255);
    se_ege2 = getAUC_CI(AHA_segs_cmr_EGE2_array_acute, AHA_segs_lrt_EGE2_array_acute, AUC_chronic);

    %title(cat(2, 'Threshold = ', num2str(thresh)));
    txt = cat(2, 'AUC = ', num2str(AUC, formatSpec),', SE = ', num2str(1.96*se, formatSpec));
    text(0.3, 0.5, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_acute, formatSpec),', SE = ', num2str(1.96*se_ege1, formatSpec));
    text(0.3, 0.4, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_chronic, formatSpec),', SE = ', num2str(1.96*se_ege2, formatSpec));
    text(0.3, 0.3, txt);
end

%% Produce quality plots （No Total)
filename = 'ROC_Total'
thresh_array = [0.01];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[100 100 650 600]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    % [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_LGE_array>thresh,AHA_segs_lrt_LGE_array,posclass);
    % plot(X, Y, '-.', 'LineWidth', 4, 'Color', 'k');
    % se = getAUC_CI(AHA_segs_cmr_LGE_array, AHA_segs_lrt_LGE_array, AUC);
    hold on;
    [X_acute,Y_acute,T_acute,AUC_acute,OPTROCPT_acute] = perfcurve(AHA_segs_cmr_LGE_array_acute>thresh,AHA_segs_lrt_LGE_array_acute,posclass);
    plot(X_acute, Y_acute, 'LineWidth', 4, 'Color', [179 0 0]/255);
    se_acute = getAUC_CI(AHA_segs_cmr_LGE_array_acute, AHA_segs_lrt_LGE_array_acute, AUC_acute);

    [X_chronic,Y_chronic,T_chronic,AUC_chronic,OPTROCPT_chronic] = perfcurve(AHA_segs_cmr_LGE_array_chronic>thresh,AHA_segs_lrt_LGE_array_chronic,posclass);
    plot(X_chronic, Y_chronic, 'LineWidth', 4, 'Color', [4 90 141]/255);
    se_chronic = getAUC_CI(AHA_segs_cmr_LGE_array_chronic, AHA_segs_lrt_LGE_array_chronic, AUC_chronic);

    %title(cat(2, 'Threshold = ', num2str(thresh)));

    % txt = cat(2, 'AUC = ', num2str(AUC, formatSpec), ', SE = ', num2str(1.96*se, formatSpec));
    % text(0.3, 0.5, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_acute, formatSpec), ', SE = ', num2str(1.96*se_acute, formatSpec));
    text(0.3, 0.4, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_chronic, formatSpec), ', SE = ', num2str(1.96*se_chronic, formatSpec));
    text(0.3, 0.3, txt);
end
%saveas(gcf, fullfile(OutputPath, [filename, '_MI.png']))

%% HEMO
filename = 'ROC_Total'
thresh_array = [0.01];
formatSpec = '%.4f';
% 1. 5% for Ground Truth
posclass = 1;
figure('Position',[100 100 650 600]);
for i = 1:length(thresh_array)
    thresh = thresh_array(i);
    % [X,Y,T,AUC,OPTROCPT] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array>thresh,AHA_segs_lrt_T2MULTIECHO_array,posclass);
    % plot(X, Y, '-.', 'LineWidth', 4, 'Color', 'k');
    % se = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array, AHA_segs_lrt_T2MULTIECHO_array, AUC);

    hold on;
    [X_acute,Y_acute,T_acute,AUC_acute,OPTROCPT_acute] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_acute>thresh,AHA_segs_lrt_T2MULTIECHO_array_acute,posclass);
    plot(X_acute, Y_acute, 'LineWidth', 4, 'Color', [179 0 0]/255);
    se_acute = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array_acute, AHA_segs_lrt_T2MULTIECHO_array_acute, AUC_acute);

    [X_chronic,Y_chronic,T_chronic,AUC_chronic,OPTROCPT_chronic] = perfcurve(AHA_segs_cmr_T2MULTIECHO_array_chronic>thresh,AHA_segs_lrt_T2MULTIECHO_array_chronic,posclass);
    plot(X_chronic, Y_chronic, 'LineWidth', 4, 'Color', [4 90 141]/255);
    se_chronic = getAUC_CI(AHA_segs_cmr_T2MULTIECHO_array_chronic, AHA_segs_lrt_T2MULTIECHO_array_chronic, AUC_chronic);

    %title(cat(2, 'Threshold = ', num2str(thresh)));
    % txt = cat(2, 'AUC = ', num2str(AUC, formatSpec), ', SE = ', num2str(1.96*se, formatSpec));
    % text(0.3, 0.5, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_acute, formatSpec), ', SE = ', num2str(1.96*se_acute, formatSpec));
    text(0.3, 0.4, txt);
    txt = cat(2, 'AUC = ', num2str(AUC_chronic, formatSpec), ', SE = ', num2str(1.96*se_chronic, formatSpec));
    text(0.3, 0.3, txt);
end

%%
addpath('../function/AUC_CI/');
s = randn(50,1) + 1; n = randn(50,1); % simulate binormal model, "signal" and "noise"
% Estimate the AUC and calculate bootstrapped 95% confidence intervals (bias-corrected and accelerated)
[A,Aci] = auc(format_by_class(AHA_segs_cmr_LGE_array>0.01,AHA_segs_lrt_LGE_array),0.05,'boot',1000,'type','bca');


function se = getAUC_CI(cmr_array, lrt_array, AUC)

N1 = sum(cmr_array > 0.01);
N2 = length(cmr_array) - N1;

Q1 = AUC / (2-AUC);
Q2 = 2*AUC^2 / (1 + AUC);
se = sqrt((AUC*(1-AUC) + (N1-1)*(Q1 - AUC^2) + (N2-1)*(Q2-AUC^2)) / (N1*N2));

end