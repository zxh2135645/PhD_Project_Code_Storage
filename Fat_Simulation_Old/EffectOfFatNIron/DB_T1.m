
clear all; clc; close all;
%% DB 
Mz_0 = 1;
Mz = Mz_0;
TD = 700;
theta = 20;
T1 = 1486;
%Hemo
for i = 1:20
    M0 = Mz * 0.84 + (1 - Mz * 0.84) * (1 - exp(-TD/1486));
    Mz = M0 * cos(theta*pi/180);
    Mxy(i) = M0 * sin(theta*pi/180);
end
%remote 
Mz_remote = 1;
for i = 1:20
    M0 = Mz_remote * 0.84 + (1 - Mz_remote * 0.84) * (1 - exp(-TD/1249));
    Mz_remote = M0 * cos(theta*pi/180);
    Mxy_remote(i) = M0 * sin(theta*pi/180);
end

contrast = Mxy_remote - Mxy;
plot(Mxy);
hold on;
plot(Mxy_remote);

%% BB
Mz_1 = 1;

for i = 1:20
    M0 = Mz_1  + (1 - Mz_1 ) * (1 - exp(-TD/1486));
    Mz_1 = M0 * cos(theta*pi/180);
    Mxy_1(i) = M0 * sin(theta*pi/180);
end
%remote 
Mz_remote_1 = 1;
for i = 1:20
    M0 = Mz_remote_1 + (1 - Mz_remote_1 ) * (1 - exp(-TD/1249));
    Mz_remote_1 = M0 * cos(theta*pi/180);
    Mxy_remote_1(i) = M0 * sin(theta*pi/180);
end

contrast_1 = Mxy_remote_1 - Mxy_1;

plot(Mxy_1);
hold on;
plot(Mxy_remote_1);
