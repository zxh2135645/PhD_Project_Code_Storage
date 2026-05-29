% Load data first
% subject_name_cell = {'18P90', '18P93', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48', '20P40'};
% avg_num_cell = {'Avg0016', 'Invivo'};

% rootDir = '/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/data';
%
% folders = dir(rootDir);
% folders = folders([folders.isdir]);
% folders = folders(~ismember({folders.name}, {'.','..'}));

rootDir = '/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/data';

subject_name_cell = {'18P90', '18P93', '20P10_Exvivo7', '20P11_Exvivo6', ...
    '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48', '20P40'};

avg_num_cell = {'Avg0016', 'Invivo'};

% number of resolution sets
nRes.Avg0016 = 28;
nRes.Invivo  = 20;

% initialize
scores = struct();
scores.Avg0016 = cell(nRes.Avg0016, 1);
scores.Invivo  = cell(nRes.Invivo, 1);

for r = 1:nRes.Avg0016
    scores.Avg0016{r} = [];
end
for r = 1:nRes.Invivo
    scores.Invivo{r} = [];
end

labels_all = [];
patientID_all = [];
subject_loaded = {};

pid = 0;

for i = 1:numel(subject_name_cell)

    subjectName = subject_name_cell{i};
    subjectPath = fullfile(rootDir, subjectName);

    if ~isfolder(subjectPath)
        warning('Missing folder: %s', subjectName);
        continue;
    end

    % -------- load Avg0016 first, because row 1 will define labels --------
    avgFile = fullfile(subjectPath, 'aha_analysis_Avg0016.mat');
    if ~isfile(avgFile)
        warning('Missing file: %s', avgFile);
        continue;
    end

    S_avg = load(avgFile);
    if ~isfield(S_avg, 'aha_analysis') || ~isfield(S_avg.aha_analysis, 'perc_array_mi')
        warning('Invalid structure in %s', avgFile);
        continue;
    end

    perc_avg = S_avg.aha_analysis.perc_array_mi;   % expected: 28 x M

    if size(perc_avg,1) < nRes.Avg0016
        warning('%s: Avg0016 has only %d rows, expected at least %d.', ...
            subjectName, size(perc_avg,1), nRes.Avg0016);
        continue;
    end

    % define label from first resolution row of Avg0016
    gt_row = perc_avg(1,:).';   % M x 1

    % ===== choose how to make binary labels =====
    % Example: infarct if > 0
    thresh = 0.1;
    labels_use = gt_row > thresh;

    % if you want another threshold, replace above, e.g.
    % labels_use = gt_row >= 50;

    M = numel(labels_use);

    pid = pid + 1;
    subject_loaded{end+1,1} = subjectName;

    labels_all    = [labels_all; labels_use];
    patientID_all = [patientID_all; repmat(pid, M, 1)];

    % store all 28 Avg0016 resolutions as scores
    for r = 1:nRes.Avg0016
        thisScore = perc_avg(r,:).';
        if numel(thisScore) ~= M
            warning('%s: Avg0016 row %d has inconsistent segment number.', subjectName, r);
            continue;
        end
        scores.Avg0016{r} = [scores.Avg0016{r}; thisScore];
    end

    % -------- load Invivo --------
    invFile = fullfile(subjectPath, 'aha_analysis_Invivo.mat');
    if ~isfile(invFile)
        warning('Missing file: %s', invFile);
        continue;
    end

    S_inv = load(invFile);
    if ~isfield(S_inv, 'aha_analysis2') || ~isfield(S_inv.aha_analysis2, 'perc_array_mi')
        warning('Invalid structure in %s', invFile);
        continue;
    end

    perc_inv = S_inv.aha_analysis2.perc_array_mi;   % expected: 20 x M

    if size(perc_inv,1) < nRes.Invivo
        warning('%s: Invivo has only %d rows, expected at least %d.', ...
            subjectName, size(perc_inv,1), nRes.Invivo);
        continue;
    end

    if size(perc_inv,2) ~= M
        warning('%s: segment number mismatch between Avg0016 and Invivo.', subjectName);
        continue;
    end

    for r = 1:nRes.Invivo
        thisScore = perc_inv(r,:).';
        scores.Invivo{r} = [scores.Invivo{r}; thisScore];
    end
end


% Do bootstrap ROC by Patient
result = bootstrapROCByPatient(scores.Invivo{1}, labels_all, patientID_all, 2000);

fprintf('Apparent AUC = %.3f\n', result.AUC);
fprintf('Bootstrap mean AUC = %.3f\n', result.AUC_mean);
fprintf('95%% CI = [%.3f, %.3f]\n', result.AUC_CI(1), result.AUC_CI(2));

% figure;
% plot(result.FPR, result.TPR, 'LineWidth', 2);
% hold on;
% plot([0 1], [0 1], '--');
% xlabel('False Positive Rate');
% ylabel('True Positive Rate');
% title(sprintf('ROC Curve (AUC = %.3f, 95%% CI [%.3f, %.3f])', ...
%     result.AUC, result.AUC_CI(1), result.AUC_CI(2)));
% grid on;
% axis square;

nBoot = 2000;

results_Avg0016 = cell(28,1);
for r = 1:28
    fprintf('Running Avg0016 resolution %d...\n', r);
    results_Avg0016{r} = bootstrapROCByPatient(scores.Avg0016{r}, labels_all, patientID_all, nBoot);
end

results_Invivo = cell(20,1);
for r = 1:20
    fprintf('Running Invivo resolution %d...\n', r);
    results_Invivo{r} = bootstrapROCByPatient(scores.Invivo{r}, labels_all, patientID_all, nBoot);
end
%%  
auc_avg = cellfun(@(x) x.AUC, results_Avg0016);
ci_avg  = cell2mat(cellfun(@(x) x.AUC_CI(:)', results_Avg0016, 'UniformOutput', false));

auc_inv = cellfun(@(x) x.AUC, results_Invivo);
ci_inv  = cell2mat(cellfun(@(x) x.AUC_CI(:)', results_Invivo, 'UniformOutput', false));

%% Figures for Avg0016 and Invivo
lb = 0.5;
figure('Position', [100 0 1600 1600]);

k = ([229, 240, 248] - [0, 113, 188]) / (1 - lb);

CI95_array_avg16 = zeros(28, 2);
auc_mean_array_avg16 = zeros(28, 1);

for i = 1:28
    
    %scores_i = scores.Avg0016{i};
    
    % ===== ROC =====
    %[X,Y,~,AUC] = perfcurve(labels, scores_i, 1);
    
    % ===== Bootstrap =====
    %[AUC_boot, CI95] = bootstrapAUC(scores_i, labels, patientID, 1000);
    
    CI95_array_avg16(i,:) = ci_avg(i,:);
    auc_mean_array_avg16(i) = auc_avg(i);
    
    subplot(4,7,i);
    plot(results_Avg0016{i}.FPR, results_Avg0016{i}.TPR, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    %plot(X,Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    
    % ===== Text =====
    text(0.9,0.5, num2str(round(results_Avg0016{i}.AUC, 2)), ...
        'FontSize',22,'HorizontalAlignment','right');
    
    text(0.95,0.3, ...
        ['[', num2str(round(results_Avg0016{i}.AUC_CI(1),2), '%.2f'), ',', num2str(round(results_Avg0016{i}.AUC_CI(2),2),'%.2f'), ']'], ...
        'FontSize',20,'HorizontalAlignment','right');
    
    % ===== Style =====
    set(gca,'FontSize',16,'Xcolor','w','Ycolor','w');
    set(gca,'XTick',[],'YTick',[]);
    
    pos = get(gca,'Position');
    pos(1) = 0.1 + mod(i-1,7)*0.09;
    pos(2) = 0.8 - fix((i-1)/7)*0.165;
    set(gca,'Position',pos);
    
    % ===== Background color =====
    if results_Avg0016{i}.AUC < lb
        set(gca,'color',[229,240,248]/255)
    else
        set(gca,'color',([229,240,248]-k*(results_Avg0016{i}.AUC-lb))/255)
    end
end

%% Invivo
lb = 0.5;
figure('Position', [100 0 1000 1600]);

k = ([229, 240, 248] - [0, 113, 188]) / (1 - lb);

CI95_array_invivo = zeros(20, 2);
auc_mean_array_invivo = zeros(20, 1);

for i = 1:20
    
    %scores_i = scores.Avg0016{i};
    
    % ===== ROC =====
    %[X,Y,~,AUC] = perfcurve(labels, scores_i, 1);
    
    % ===== Bootstrap =====
    %[AUC_boot, CI95] = bootstrapAUC(scores_i, labels, patientID, 1000);
    
    CI95_array_invivo(i,:) = ci_inv(i,:);
    auc_mean_array_invivo(i) = auc_inv(i);
    
    subplot(4,5,i);
    plot(results_Invivo{i}.FPR, results_Invivo{i}.TPR, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    %plot(X,Y, 'LineWidth', 3, 'color', [246, 101, 72]/255);
    
    % ===== Text =====
    text(0.9,0.5, num2str(round(results_Invivo{i}.AUC, 2)), ...
        'FontSize',22,'HorizontalAlignment','right');
    
    text(0.95,0.3, ...
        ['[', num2str(round(results_Invivo{i}.AUC_CI(1),2), '%.2f'), ',', num2str(round(results_Invivo{i}.AUC_CI(2),2),'%.2f'), ']'], ...
        'FontSize',20,'HorizontalAlignment','right');
    
    % ===== Style =====
    set(gca,'FontSize',16,'Xcolor','w','Ycolor','w');
    set(gca,'XTick',[],'YTick',[]);
   
    pos = get(gca, 'Position');
    pos(1) = 0.1 + mod(i-1, 5) * 0.127;
    pos(2) = 0.8 - fix((i-1)/5) * 0.165;
    set(gca,'Position',pos);

    % ===== Background color =====
    if results_Invivo{i}.AUC < lb
        set(gca,'color',[229,240,248]/255)
    else
        set(gca,'color',([229,240,248]-k*(results_Invivo{i}.AUC-lb))/255)
    end
end

%% Invivo
AUC_avg16 = [1,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan];
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
%% Flatten the voxel size, plot vs AUC (version 2 for SCMR)
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 300 600]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
%mdl = fitnlm(tbl,modelfun,beta0);


% ci = coefCI(mdl);
% b = mdl.Coefficients.Estimate;

%Y_pred = modelfun(b, X)
B_avg16 = B;
B_avg16(1) = 0.3 * 0.3 * 2;
plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
%plot(X, Y_pred); %ylim([0.5 1])
% Y_lb = modelfun(ci(:,1), X);
% Y_ub = modelfun(ci(:,2), X);
%plot(X, Y_lb);
%plot(X, Y_ub);
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_invivo(1,:),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_invivo(2,:),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_invivo(3,:),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_invivo(4,:),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%% Avg0016
AUC_avg16 = [1,nan,nan,nan,nan,nan,nan;nan,nan,nan,nan,nan,nan,nan;nan,nan,nan,nan,nan,nan,nan;nan,nan,nan,nan,nan,nan,nan];
AUC_invivo_allavg = [nan,0.96,0.94,0.93,0.86,0.82,0.70;0.91,0.91,0.87,0.86,0.81,0.78,0.69;0.83,0.81,0.80,0.76,0.75,0.73,0.67;0.74,0.73,0.72,0.72,0.69,0.67,0.66];
%% Flatten the voxel
inplane_res = [0.3, 0.6, 0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 300 600]);
Y = AUC_avg16_array(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
%mdl = fitnlm(tbl,modelfun,beta0);


% ci = coefCI(mdl);
% b = mdl.Coefficients.Estimate;

%Y_pred = modelfun(b, X)
B_avg16 = B;
B_avg16(1) = 0.3 * 0.3 * 2;
plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
%plot(X, Y_pred); %ylim([0.5 1])
% Y_lb = modelfun(ci(:,1), X);
% Y_ub = modelfun(ci(:,2), X);
%plot(X, Y_lb);
%plot(X, Y_ub);
set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_invivo(1,:),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_invivo(2,:),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_invivo(3,:),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_invivo(4,:),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');


% dataStruct = struct();
%
% for k = 1:length(avg_num_cell)
%     key = avg_num_cell{k};
%     dataStruct.(key).scores = [];
%     dataStruct.(key).labels = [];
%     dataStruct.(key).patientID = [];
% end
%
% % scores = [];    % predicted values (Nx1)
% % labels = [];    % 0/1 ground truth
% % patientID = []; % same length, indicates cluster
%
% pid = 0;
%
% for i = 1:length(subject_name_cell)
%
%     subjectName = subject_name_cell{i};
%     subjectPath = fullfile(rootDir, subjectName);
%
%     if ~isfolder(subjectPath)
%         warning('Missing folder: %s', subjectName);
%         continue;
%     end
%
%     pid = pid + 1;
%
%     for k = 1:length(avg_num_cell)
%
%         avg_tag = avg_num_cell{k};
%
%         % find matching file
%         pattern = sprintf('aha_analysis_%s.mat', avg_tag);
%         f = dir(fullfile(subjectPath, pattern));
%
%         if isempty(f)
%             warning('Missing %s for %s', avg_tag, subjectName);
%             continue;
%         end
%
%         % load
%         tmp = load(fullfile(subjectPath, f(1).name));
%
%         if ~isfield(tmp, 'aha_analysis2')
%             warning('aha_analysis2 missing in %s', f(1).name);
%             continue;
%         end
%
%         perc = tmp.aha_analysis2.perc_array_mi; % N x M
%
%
%         % ===== Select resolution =====
%         resIdx = 1;  % <-- CHANGE if needed
%         perc_use = perc(resIdx, :)';
%
%         % ===== Define labels (YOU SHOULD REPLACE THIS) =====
%         thresh = 0.1;
%         labels_use = perc_use > thresh;
%
%         % ===== Append =====
%         dataStruct.(avg_tag).scores = [dataStruct.(avg_tag).scores; perc_use];
%         dataStruct.(avg_tag).labels = [dataStruct.(avg_tag).labels; labels_use];
%         dataStruct.(avg_tag).patientID = [dataStruct.(avg_tag).patientID; ...
%             repmat(pid, length(perc_use), 1)];
%     end
% end