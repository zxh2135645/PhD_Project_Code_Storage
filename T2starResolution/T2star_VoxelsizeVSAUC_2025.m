%% T2* value in different voxel size
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';
idx_mat = repmat(1:5, [4,1]) + 5*repmat([0:3]', [1,5]);

vol_array = vol_mat(:);
idx_array = idx_mat(:);

[B,I] = sort(vol_array)
BI = idx_array(I)

% 6, 11,12,17,18,23

AUC_avg16 = [1,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan];
AUC_avg16 = [nan,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan;nan,nan,nan,nan,nan];

%AUC_avg16 = [1, 0.96, 0.94, 0.92, 0.86, 0.8, 0.71; 0.92, 0.88, 0.88, 0.85, 0.8, 0.76, 0.69; ...
%    0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.7; 0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];
AUC_invivo = [0.57,0.68,0.68,0.77,0.64;0.69,0.76,0.76,0.75,0.66;0.76,0.76,0.71,0.69,0.6;0.74,0.69,0.66,0.66,0.58];
%05/29/2023
AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
AUC_invivo_subjectavg = [0.56,0.67,0.68,0.78,0.64;0.67,0.67,0.70,0.73,0.62;0.65,0.68,0.63,0.62,0.54;0.64,0.58,0.60,0.62,0.52];

AUC_avg16_array = AUC_avg16(:);
AUC_invivo_array = AUC_invivo(:);

figure(); plot(B, AUC_avg16_array(I),'o'); grid on;
figure(); plot(B, AUC_invivo_array(I),'o'); grid on;

%% Flatten the voxel size, plot vs AUC (version 2 for SCMR)
inplane_res = [0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

AUC_invivo = AUC_invivo_allavg;
figure('Position', [100 0 300 300]);
Y = AUC_invivo_array(I);
X = B;
tbl = table(X, Y);
%modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
% modelfun = @(b,x) b(1) + b(2)*x.^b(3);

% gamma function PDF
% k = b(1); lambda = b(2);
modelfun = @(b,x) b(2)^b(1) ./ (gamma(b(1))) .* x.^(b(1)-1) .* exp(-b(2)*x);
beta0 = [0 0];
mdl = fitnlm(tbl,modelfun,beta0);


ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
%B_avg16 = B;
%B_avg16(1) = 0.3 * 0.3 * 2;
%plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
%plot(X, Y_pred); %ylim([0.5 1])
Y_lb = modelfun(ci(:,1), X);
Y_ub = modelfun(ci(:,2), X);
plot(X, Y_lb);
plot(X, Y_ub);

default_color = orderedcolors("gem");
% set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
%     'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,1) = plot(vol_mat(1,:), AUC_invivo(1,:),'-', 'Color', default_color(1,:));
set(plotHandles_auc(:,1), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,2) = plot(vol_mat(2,:), AUC_invivo(2,:),'-', 'Color', default_color(2,:));
set(plotHandles_auc(:,2), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,3) = plot(vol_mat(3,:), AUC_invivo(3,:),'-', 'Color', default_color(3,:));
set(plotHandles_auc(:,3), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,4) = plot(vol_mat(4,:), AUC_invivo(4,:),'-', 'Color', default_color(4,:));
set(plotHandles_auc(:,4), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
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
%%
% x, y are your scatter vectors
M = fit_gamma_like(X, y, true);
disp(M.bestName)
disp(M.params)

% Predict on new points:
xnew = linspace(min(x), max(x), 100).';
ynew = M.predict(xnew);

%%
X = [1.28 2 2.56 3.38 3.84 4 5.12 5.12 6 6.76 8 8.82 10.14 10.24 13.52 15.36 17.64 20.48 26.46 35.28]';
Y = [0.57 0.68 0.69 0.68 0.76 0.76 0.74 0.77 0.76 0.76 0.69 0.64 0.71 0.75 0.66 0.69 0.66 0.66 0.60 0.58]';
% X = [13.52 15.36 17.64 20.48 26.46 35.28]';
% Y = [0.66 0.69 0.66 0.66 0.60 0.58]';

n = numel(X);
degList = 2:2;                    % try linear up to quintic
best = struct('deg',[], 'AICc',inf);

results = struct('deg',{},'p',{},'mu',{},'rss',{},'AICc',{},'rmse',{},'r2',{});
for d = degList
    % ... compute p, mu, rss, AICc, rmse, r2 ...
    results(end+1) = struct('deg',d,'p',p,'mu',mu,'rss',rss,'AICc',AICc,'rmse',rmse,'r2',r2);
end

for d = degList
    [p,S,mu] = polyfit(X,Y,d);    % numerically stable
    Yhat = polyval(p,X,[],mu);
    rss  = sum((Y - Yhat).^2);
    k    = d + 1;
    AIC  = n*log(rss/n) + 2*k;
    AICc = AIC + (2*k*(k+1))/(n-k-1);    % small-sample correction
    rmse = sqrt(mean((Y-Yhat).^2));
    r2   = 1 - rss/sum((Y - mean(Y)).^2);

    results(end+1) = struct('deg',d,'p',p,'mu',mu,'rss',rss,'AICc',AICc,'rmse',rmse,'r2',r2); %#ok<SAGROW>
    if AICc < best.AICc
        best = results(end);
    end
end

fprintf('Best degree by AICc: %d\n', best.deg);
fprintf('RMSE = %.5f, R^2 = %.3f, AICc = %.3f\n', best.rmse, best.r2, best.AICc);

% Plot
xg = linspace(min(X), max(X), 400)';
yg = polyval(best.p, xg, [], best.mu);

figure; hold on
scatter(X,Y,28,'filled'); 
plot(xg,yg,'LineWidth',2);
grid on
xlabel('X'); ylabel('Y'); title(sprintf('Best polynomial fit (degree %d)', best.deg));
legend('Data','Polynomial fit','Location','best');
%%
xg_1 = xg;
xg_2 = yg;

%%
default_color = orderedcolors("gem");
figure('Position', [100 0 300 300]);
plot(xg_1(xg_1<13), xg_2(xg_1<13), 'LineWidth', 2, 'Color',default_color(2,:));

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
%% T2* value in different voxel size
inplane_res = [0.3, 0.6, 0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';
idx_mat = repmat(1:7, [4,1]) + 7*repmat([0:3]', [1,7]);

vol_array = vol_mat(:);
idx_array = idx_mat(:);

[B,I] = sort(vol_array)
BI = idx_array(I)

% 6, 11,12,17,18,23

AUC_avg16 = [1, 0.96, 0.94, 0.92, 0.86, 0.8, 0.71; 0.92, 0.88, 0.88, 0.85, 0.8, 0.76, 0.69; ...
    0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.7; 0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];
AUC_invivo = [0.57,0.68,0.68,0.77,0.64;0.69,0.76,0.76,0.75,0.66;0.76,0.76,0.71,0.69,0.6;0.74,0.69,0.66,0.66,0.58];
%05/29/2023
%AUC_invivo_allavg = [0.57,0.69,0.72,0.82,0.66;0.73,0.74,0.80,0.77,0.65;0.74,0.77,0.72,0.69,0.59;0.71,0.67,0.68,0.70,0.57];
%AUC_invivo_subjectavg = [0.56,0.67,0.68,0.78,0.64;0.67,0.67,0.70,0.73,0.62;0.65,0.68,0.63,0.62,0.54;0.64,0.58,0.60,0.62,0.52];

AUC_avg16_array = AUC_avg16(:);

figure(); plot(B, AUC_avg16_array(I),'o'); grid on;

%% Flatten the voxel size, plot vs AUC (version 2 for SCMR)
inplane_res = [0.3, 0.6, 0.8, 1.0, 1.3, 1.6, 2.1];
thrplane_res = [2, 4, 6, 8];

vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';


figure('Position', [100 0 300 300]);
Y = AUC_avg16(I);
X = B;
tbl = table(X, Y);
modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
beta0 = [0 0 0];
mdl = fitnlm(tbl,modelfun,beta0);


ci = coefCI(mdl);
b = mdl.Coefficients.Estimate;

Y_pred = modelfun(b, X)
%B_avg16 = B;
%B_avg16(1) = 0.3 * 0.3 * 2;
%plotHandles_auc(:,1) = plot(B_avg16, AUC_avg16_array(I),'o'); grid on;
hold on;
% plot(X, Y_pred); %ylim([0.5 1])
Y_lb = modelfun(ci(:,1), X);
Y_ub = modelfun(ci(:,2), X);
%plot(X, Y_lb);
%plot(X, Y_ub);

default_color = orderedcolors("gem");
% set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
%     'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,1) = plot(vol_mat(1,:), AUC_avg16(1,:),'-', 'Color', default_color(1,:));
set(plotHandles_auc(:,1), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,2) = plot(vol_mat(2,:), AUC_avg16(2,:),'-', 'Color', default_color(2,:));
set(plotHandles_auc(:,2), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,3) = plot(vol_mat(3,:), AUC_avg16(3,:),'-', 'Color', default_color(3,:));
set(plotHandles_auc(:,3), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);

plotHandles_auc(:,4) = plot(vol_mat(4,:), AUC_avg16(4,:),'-', 'Color', default_color(4,:));
set(plotHandles_auc(:,4), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12, ...
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

%%
figure('Position', [100 0 300 300]);
plot(X, Y_pred, 'LineWidth', 2);

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