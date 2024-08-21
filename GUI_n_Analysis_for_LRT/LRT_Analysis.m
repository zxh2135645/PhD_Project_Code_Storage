clear all;
close all;

%% Linear regression and Bland-Altman
%% Cardiac Function (Ejection Fraction)
color_cell1 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell2 = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};


addpath('../function/BlandAltman/');
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1', 'WK8'};
nstates = length(states);

% Data preparation
%            Sofia Lisbon Paris  Jesse  George Chili  Nutmeg Ginger Dave  Paprika Cinnamon           
lrt_acute = [35.98, 32.1, 27.35, 36.89, 34.87, 32.05, 24.57, 39.93, 24.62, 27.25, 33.31];
%             Sofia  Lisbon  Paris  Jesse  George Nutmeg Ginger Cinnamon
lrt_chronic = [39.74, 39.74, 28.27, 30.54, 31.87, 35.02, 33.97, 33.97,nan,nan,nan];
cmr_acute = [37, 31.55, 28.51, 38.33, 34.57, 32.31, 23, 39, 24, 27, 34];
cmr_chronic = [38.8, 41, 28.67, 30.88, 34, 35, 35, 33,nan,nan,nan];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2];
% Baseline data with noise
data1 = cat(3, cmr_acute(:), cmr_chronic(:));
data2 = cat(3, lrt_acute(:), lrt_chronic(:));

% BA plot paramters
tit = 'Ejection Fraction (EF) Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'CMR-EF','LRT-EF','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [20 40 20 40];
% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'AxesLimits', axesLimits);


data1_acute = cat(3, cmr_acute(:));
data2_acute = cat(3, lrt_acute(:));

% BA plot paramters
states_acute = {'WK1'};
gnames = {territories, states_acute}; % names of groups in data {dimension 1 and 2}
[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_chronic = cat(3, cmr_chronic(:));
data2_chronic = cat(3, lrt_chronic(:));

% BA plot paramters
states_chronic = {'WK8'};
%colors = 'rb'; 
colors = [color_cell2{5};... % or RGB triplets
		      color_cell1{5}];
cc = [color_cell2{2};...
    color_cell1{3}];
gnames = {territories, states_chronic}; % names of groups in data {dimension 1 and 2}

[cr_chronic, fig_chronic, statsStruct_chronic] = BlandAltman_Customized(data1_chronic, data2_chronic,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', cc,...
    'AxesLimits', axesLimits);


%% Infarct Size

% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1', 'WK8'};
nstates = length(states);

% Data preparation
%            Sofia Lisbon Paris  Jesse  George Chili  Nutmeg Ginger Dave  Paprika Cinnamon Carlos          
lrt_acute = [9.43, 18.27, 19.15, 18.54, 23.55, 22.67, 20.36, 23.31, 24.98, 10.65, 21.06,   47];
%             Sofia  Lisbon  Paris  Jesse  George Nutmeg Ginger Paprika Cinnamon
lrt_chronic = [8.84, 11.18,  8.23,  12.78, 16.99, 10.08, 10.6,  8.46,   12.05,nan,nan,nan];
cmr_acute = [8.27, 19.13, 18.13, 20.69, 22.53, 21.68, 19.37, 20.72, 22.94, 11.34, 23.17,   43.97];
cmr_chronic = [8.72, 11.01, 7.11, 13.1, 16.01, 11.7, 10.32, 9.23, 14.91, nan,nan,nan];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2];
% Baseline data with noise
data1 = cat(3, cmr_acute(:), cmr_chronic(:));
data2 = cat(3, lrt_acute(:), lrt_chronic(:));

% BA plot paramters
tit = 'MI Size Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'CMR-MI size','LRT-MI size','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [0 50 0 50];

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_acute = cat(3, cmr_acute(:));
data2_acute = cat(3, lrt_acute(:));

% BA plot paramters
states_acute = {'WK1'};
gnames = {territories, states_acute}; % names of groups in data {dimension 1 and 2}

[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);

data1_chronic = cat(3, cmr_chronic(:));
data2_chronic = cat(3, lrt_chronic(:));

% BA plot paramters
states_chronic = {'WK8'};
colors = [color_cell2{5};... % or RGB triplets
		      color_cell1{5}];
cc = [color_cell2{2};...
    color_cell1{3}];
gnames = {territories, states_chronic}; % names of groups in data {dimension 1 and 2}

[cr_chronic, fig_chronic, statsStruct_chronic] = BlandAltman_Customized(data1_chronic, data2_chronic,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', cc,...
    'AxesLimits', axesLimits);

%% Infarct Transmurality
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1', 'WK8'};
nstates = length(states);

% Data preparation
%            Sofia Lisbon Paris  Jesse  George Chili  Nutmeg Ginger Dave  Carlos Paprika Cinnamon           
lrt_acute = [43.15, 46.03, 51.03, 53.3, 61.91, 67.94, 68.4,  63.06, 69.84, 81,   55.14,   62.82];
%             Sofia  Lisbon  Paris  Jesse  George Nutmeg Ginger Paprika Cinnamon
lrt_chronic = [43.41, 44.91, 50.19, 52.94, 62.67, 67.53, 65.81,  59.5,   56.08, nan,nan,nan];
cmr_acute = [41.87, 49.56, 42.62, 60, 64.74, 61.26, 64.16, 66.5, 72.46, 76.91, 56.16,  61.35];
cmr_chronic = [43.64, 45.65, 46.86, 49.98, 61.37, 59.09, 58.95, 61.32,  51.84, nan,nan,nan];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2];
% Baseline data with noise
data1 = cat(3, cmr_acute(:), cmr_chronic(:));
data2 = cat(3, lrt_acute(:), lrt_chronic(:));

% BA plot paramters
tit = 'MI Transmurality Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'CMR-MI','LRT-MI','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
axesLimits = [40 90 40 90];

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_acute = cat(3, cmr_acute(:));
data2_acute = cat(3, lrt_acute(:));

% BA plot paramters
states_acute = {'WK1'};
gnames = {territories, states_acute}; % names of groups in data {dimension 1 and 2}

[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_chronic = cat(3, cmr_chronic(:));
data2_chronic = cat(3, lrt_chronic(:));

% BA plot paramters
states_chronic = {'WK8'};
colors = [color_cell2{5};... % or RGB triplets
		  color_cell1{5}];
cc = [color_cell2{2};...
    color_cell1{3}]; 
gnames = {territories, states_chronic}; % names of groups in data {dimension 1 and 2}

[cr_chronic, fig_chronic, statsStruct_chronic] = BlandAltman_Customized(data1_chronic, data2_chronic,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', cc,...
    'AxesLimits', axesLimits);


%% MVO
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1', 'WK8'};
nstates = length(states);

% Data preparation
%           Lisbon Jesse George Chili Nutmeg Ginger Dave  Carlos Paprika            
lrt_acute = [4.35, 5.79, 9.72,  4.81, 3.68,  4.3,   7.48, 20.15, 2.17];
%           Lisbon Jesse George Chili Nutmeg Ginger Dave  Carlos Paprika
cmr_acute = [3.13, 8.08, 7.34,  5.93, 2.53,  3.53,  4.67, 17.96, 1.81];

% Baseline data with noise
data1_acute = cat(3, cmr_acute(:));
data2_acute = cat(3, lrt_acute(:));

% BA plot paramters
tit = 'MVO Size Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'CMR-MVO','LRT-MVO','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end
% axesLimits = [40 90 40 90];

% BA plot paramters
[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48);

%% Hemorrhage
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1', 'WK8'};
nstates = length(states);

% Data preparation
%            Lisbon Jesse  George Chili  Nutmeg Ginger Dave  Carlos Paprika          
lrt_acute = [2.2,   7.65,  9.01,  5.62,  5.82,  6.63,  7.16, 12.26, 2.84];
%             Lisbon  Jesse  George Nutmeg Ginger Paprika
% lrt_chronic = [4.81,     3.96,  5.04,  3.34,  4.22,  1.52, nan, nan, nan];
lrt_chronic = [4.81,     3.96,  5.04,  3.34,  4.22,  1.52, 0, 0, 0];
cmr_acute = [2.62, 7.99, 9.78, 7.79, 5.59, 6.39, 6.34, 11.28, 2.81];
% cmr_chronic = [5.30, 3.02, 5.07, 3.43, 3.2, 1.58, nan, nan, nan];
cmr_chronic = [5.30, 3.02, 5.07, 3.43, 3.2, 1.58, 0, 0, 0];

% Baseline data with noise
data1 = cat(3, cmr_acute(:), cmr_chronic(:));
data2 = cat(3, lrt_acute(:), lrt_chronic(:));

% BA plot paramters
tit = 'Hemorrhage Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'CMR-Hemorrhage','LRT-Hemorrhage','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_cell1{5};... % or RGB triplets
		      color_cell2{5}];
end

axesLimits = [0 14 0 14];

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman_Customized(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_acute = cat(3, cmr_acute(:));
data2_acute = cat(3, lrt_acute(:));

% BA plot paramters
states_acute = {'WK1'};
gnames = {territories, states_acute}; % names of groups in data {dimension 1 and 2}

[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48,...
    'AxesLimits', axesLimits);


data1_chronic = cat(3, cmr_chronic(:));
data2_chronic = cat(3, lrt_chronic(:));

% BA plot paramters
states_chronic = {'WK8'};
colors = [color_cell2{5};... % or RGB triplets
		  color_cell1{5}];
cc = [color_cell2{2};...
    color_cell1{3}]; 
gnames = {territories, states_chronic}; % names of groups in data {dimension 1 and 2}

[cr_chronic, fig_chronic, statsStruct_chronic] = BlandAltman_Customized(data1_chronic, data2_chronic,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', cc,...
    'AxesLimits', axesLimits);

