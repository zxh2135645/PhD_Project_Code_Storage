clear all;
close all;

sigmoid = @(x) 1 ./ (1 + exp(-x));
fx = @(x, x0) 1 - 1 ./ (1 + exp(-(x-x0)));

X = -40:0.1:40;
X0 = -20:0.4:25;


figure();
Y = sigmoid(X);
h1 = plot(X,Y, 'LineWidth', 2);
hold on;

n = length(X0);
for i = 1:length(X0);

%for i = n:n

    Y2 = fx(X,X0(i));
    h2 = plot(X,Y2, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
    Y3 = Y.*Y2;
    %Y3 = conv(Y,Y2);
    h3 = plot(X,Y3, '--', 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250]); 
    xlim([-40 40]) ;ylim([0 1]);
    pause(0.2);
    if i ~= n
        set(h2,'visible','off');
        %set(h3,'visible','off');
    end
    axis off;
end
%% Variation of Senario 3
sigmoid = @(x, k1) 1 ./ (1 + exp(-k1*x));
fx = @(x, k2, x0) 1 - 1 ./ (1 + exp(-k2*(x-x0)));
X = -0:0.1:80;
X0 = -20:0.4:100;

k1 = 1;
k2 = 1/16;
figure();
Y = sigmoid(X,k1);
h1 = plot(X,Y, 'LineWidth', 3);
hold on;

n = 300;
%for i = 1:length(X0);

for i = n:n

    Y2 = fx(X,k2,X0(i));
    h2 = plot(X,Y2, 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]); 
    Y3 = Y.*Y2;
    %Y3 = conv(Y,Y2);
    h3 = plot(X,Y3, '--', 'LineWidth', 4, 'Color', [0.9290, 0.6940, 0.1250]); 
    xlim([40 80]) ;ylim([0 1]);
    pause(0.2);
    if i ~= n
        set(h2,'visible','off');
        %set(h3,'visible','off');
    end
    axis off;
end

%% Variation of in scales (noise level)
sigmoid = @(x, x1, k1) 1 ./ (1 + exp(-k1*(x-x1)));
fx = @(x, k2, x0) 1 - 1 ./ (1 + exp(-k2*(x-x0)));
X = -20:0.5:100;
X0 = -20:0.4:100;

k1 = 1/4;
k2 = 1/4;
figure();
X1 = 0;
Y = sigmoid(X,X1,k1);
h1 = plot(X,Y, 'LineWidth', 3);
hold on;

n = 100;
%for i = 1:length(X0);

for i = n:n

    Y2 = fx(X,k2,X0(i));
    h2 = plot(X,Y2, 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]); 
    Y3 = Y.*Y2;
    %Y3 = conv(Y,Y2);
    h3 = plot(X,Y3, '--', 'LineWidth', 4, 'Color', [0.9290, 0.6940, 0.1250]); 
    xlim([-40 40]) ;ylim([0 1]);
    pause(0.2);
    if i ~= n
        set(h2,'visible','off');
        %set(h3,'visible','off');
    end
    axis off;
end

%% Fitting for ex-vivo practical
%% Variation of Senario 3
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-100],...
               'Upper',[1,1,100],...
               'StartPoint',[1 1,0]);
ft = fittype('1 ./ (1 + exp(-k1*x)) * (1 - 1 ./ (1 + exp(-k2*(x-x0))))','independent','x', 'coefficients',{'k1','k2','x0'}, 'options',fo);

auc_mat = [0.57, 0.68, 0.68, 0.77, 0.64;
    0.69, 0.76, 0.76, 0.75, 0.66;
    0.76, 0.76, 0.64, 0.69, 0.60;
    0.74, 0.69, 0.66, 0.66, 0.58];
snr_mat = [3.889101759, 6.796342678, 8.90610959, 11.69427981, 15.14968152;
    7.989617412, 12.37448587, 15.04384232, 18.1386487, 24.18634346;
    11.57426404, 16.37096552, 19.903348, 21.82762437, 27.80253859;
    14.41929372, 18.08539358, 21.71116137, 24.08643505, 28.5705392];
vox_mat = ([0.8; 1.0; 1.3; 1.6; 2.1].^2 * [2, 4, 6, 8]).';
[myfit, gof] = fit(vox_mat(:),auc_mat(:),ft);

X = 1:0.1:40;
figure();
plot(myfit, vox_mat(:),auc_mat(:));

hold on;
y1 = 1 ./ (1 + exp(-myfit.k1*(X)));
y2 = 1 - 1 ./ (1 + exp(-myfit.k2*(X-myfit.x0)));

plot(X,y1);
plot(X,y2);

%% Fitting for ex-vivo ideal
auc_mat = [1, 0.96, 0.94, 0.92, 0.86, 0.80, 0.71;
           0.92, 0.88, 0.88, 0.85, 0.80, 0.76, 0.69;
           0.83, 0.81, 0.82, 0.78, 0.77, 0.74, 0.70;
           0.76, 0.75, 0.74, 0.74, 0.71, 0.67, 0.67];
vox_mat = ([0.3; 0.6; 0.8; 1.0; 1.3; 1.6; 2.1].^2 * [2, 4, 6, 8]).';

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1,0,-100, -1000,-1],...
               'Upper',[1,1,100,1000,1],...
               'StartPoint',[1,1,0,0,0]);
ft = fittype('(1 ./ (1 + exp(-k1*(x-x1)))) * ((1 - 1 ./ (1 + exp(-k2*(x-x0)))) + h2)','independent','x', 'coefficients',{'k1','k2','x0','x1','h2'}, 'options',fo);
[myfit2,gof2] = fit(vox_mat(:),auc_mat(:),ft);
X = 1:0.1:40;
figure();
plot(myfit2, vox_mat(:),auc_mat(:));
hold on;
y1 = 1 ./ (1 + exp(-myfit2.k1*(X-myfit2.x1)));
y2 = 1 - 1 ./ (1 + exp(-myfit2.k2*(X-myfit2.x0))) + myfit2.h2;

plot(X,y1);
plot(X,y2);



