
%%%%%%%%%%%%%%%%%%%%%%%%%% input  file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Analysis of Iron concentration - so far T2 mapping is questionable

addpath('../function/');
labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT', 'IR'};
label = labels{3};


proj_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name = input('Please type subject name here:  ', 's');
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

load(cat(2, subject_data_dir, 'maps.mat'));

%% plot loaded non-fat, 3D surface overview
c = 0:10:50;
ff = [0,2.5,5,10,20,30,40];
delta = 1./(maps.t2star_siemens) - 1./(maps.t2_siemens);
[X,Y] = meshgrid(c, ff);
figure(); surf(X,Y,delta);
colormap(brewermap([],'*YlGnBu'));
%% and linear regression
C = [ones(length(c),1), c'];
b = C\delta';
y_h = b(2)*c + b(1);
figure();
plot(c, delta(1,:), 'LineWidth', 2);
hold on;
plot(c, y_h, 'LineWidth', 2); grid on;

Csq = [ones(length(c),1), (c.^2)'];
bsq = Csq\delta';
y_hsq = bsq(2)*c.^2 + bsq(1);

figure();
plot(c, delta(1,:), 'LineWidth', 2); hold on;
plot(c, y_hsq, 'LineWidth', 2); grid on;