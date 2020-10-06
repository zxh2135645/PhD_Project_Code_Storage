close all;
clear all;

c = [0 0.2 0.4 0.6];
T1 = [1104, 504, 301.1, 179.5];
T2 = [136, 30.3, 21.5, 20.3];
T2star = [128, 25.8, 15.7, 10.8];

% Rescanned numbers:
% 1. From MAPIT
T1 = [236.8, 177.9, 120.4 70.3];
T2 = [152.8 37.8 24.4 12.3];
T2star = [99.7 28.7 19.5 10.9];

% 2. From CMR
T1 = [1181 505.4 309.9 111.9];
T2 = [118.9 35.4 28.6 26.1];
T2star = [123.8 33.2 22.5 11.8];


y = -1./T2 + 1./T2star;
C = [ones(length(c),1), c'];

b = C\y';
y_h = b(2)*c + b(1);

figure();
plot(c, y, 'LineWidth', 2); hold on;
plot(c, y_h, 'LineWidth', 2); grid on;


Csq = [ones(length(c),1), (c.^2)'];
bsq = Csq\y';
y_hsq = bsq(2)*c.^2 + bsq(1);

figure();
plot(c.^2, y, 'LineWidth', 2); hold on;
plot(c.^2, y_hsq, 'LineWidth', 2); grid on;

% The current conclusion is that it is tricy to get 1/T2* - 1/T2

