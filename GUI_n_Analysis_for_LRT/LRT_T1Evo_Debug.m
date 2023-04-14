clear all;
close all;
%%
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');

mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Ny, Nx, Nz);
%if ~exist(mask_f)
    for i = 1:Nz
        dispim = @(x,st) fftshift(x(:,:,i,:),1);
        for j = 1:1
            for k = 1:1
                temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                cw = 0.5*max(vec(abs(temp)));
                figure();
                imagesc(abs(temp(:,:,end,1))/cw); axis image;
                roi = drawpolygon;
                mask(:,:,i) = createMask(roi);
            end
        end
    end
%    save(mask_f, 'mask');
% else
%     load(mask_f);
% end

%% 3D plot
addpath('../function/');
j = 1;
k = 1;
n = 1;
slc = 2;
figure();
hold on;
evos = zeros(sizes(5), sizes(2));
for n = 1:sizes(5)
    dispim = @(x,st) fftshift(x(:,:,slc,:),1);
    temp = Gr\reshape(Phi(:,:,j,k,n), L, []);
    temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, sizes(2));
    cw = 0.5*max(vec(abs(temp)));
    % figure();
    % imshow3D(abs(temp)./cw);
    mask(mask == 0) = nan;
    evo = mean(abs(reshape(mask(:,:,slc) .* temp, [], sizes(2))), 1, 'omitnan');

    plot(evo);
    evos(n,:) = evo;
end


N = 1:sizes(5);
x = 1:sizes(2);
[X,Y] = meshgrid(x,N);
figure();
surf(X,Y,evos); %zlim([0  0.5e-3])

%% Comparing to PostCon MOLLI
addpath('../function/BlandAltman/');
% Lisbon D8 vs Lisbon D6
% 1,  4,  6,  8, 10,  12, 14
t1_molli = [285,433,484,515,538,559,579;...
            1019,746,630,545,516,488,442;...
            138,199,237,265,288,304,320;...
            270,283,290,300,306,315,326];
t1_molli = [285,433,484,515,538,559,579;...
            138,199,237,265,288,304,320;...
            270,283,290,300,306,315,326];
t1_molli = [285,433,484,515,538,559,579;...
            138,199,237,265,288,304,320];

t1_lrt = [247,	324,	351,	374,	394,	399,	416; ...
          665,	577.667,	526.667,	476.333,	451.667,	446.667,	387.333;...
          107.806,	143.581,	166.774,	181.774,	201.226,	220.387,	228.323;...
          144,	157.8,	169.4,	176.4,	188.4,	186.8,	187.4];

t1_lrt = [247,	324,	351,	374,	394,	399,	416; ...
          107.806,	143.581,	166.774,	181.774,	201.226,	220.387,	228.323;...
          144,	157.8,	169.4,	176.4,	188.4,	186.8,	187.4];
t1_lrt = [247,	324,	351,	374,	394,	399,	416; ...
          107.806,	143.581,	166.774,	181.774,	201.226,	220.387,	228.323];% Heart muscle territories per patient

time_point = 'Acute';
territories = {time_point};
nterritories = length(territories);

% Patient states during measurement
states = {'Remote', 'MVO', 'Blood Pool', 'MI'};
states = {'Remote', 'Blood Pool', 'MI'};
states = {'Remote', 'Blood Pool'};

nstates = length(states);

data1 = reshape(t1_molli.', size(t1_molli,2), 1, []);
data2 = reshape(t1_lrt.', size(t1_lrt,2), 1, []);


% BA plot paramters
tit = 'T1 Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'T1 MOLLI','T1 LRT','ms'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [[0, 0.4470, 0.7410];... % or RGB triplets
		      %[0.8500, 0.3250, 0.0980];...
              [0.9290, 0.6940, 0.1250]];
              %[0.4940, 0.1840, 0.5560]];
end

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on');


