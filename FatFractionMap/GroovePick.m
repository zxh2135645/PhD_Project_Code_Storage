%%
% Picking up grooves slice per slice 
clear all;
close all;
clc
current_dir = pwd;
%% 
addpath('..\function\');
addpath('..\AHA16Segment\');

base_dir = GetFullPath(cat(2, current_dir, '\..\..\Data\Diane\ContourData\'));
out_dir = GetFullPath(cat(2, base_dir, '../../Diane/ResultsFolder_180718/'));

sequence_label = {'MAG', 'T1Map', 'MultiEcho'};
output_label = {'LGE', 'T1', 'MultiEcho'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'freeROI'};
anatomy = anatomy_label{1};
overwrite_label = 1;

name_glob = glob(cat(2, base_dir, '/*_*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end


for ll = 1:length(sequence_label)
    for n = 1:length(Names)
        label = sequence_label{ll};
        coords = struct;
        name = Names{n};
        outputFileName = [base_dir, name, '\', label, '_coords.mat'];
        
        switch label
            case {'MAG', 'PSIR'}
                labelo = output_label{1};
            case {'T1Map'}
                labelo = output_label{2};
            case {'MultiEcho'}
                labelo = output_label{3};
            otherwise
                disp('Sequence label is out of range');
        end
        
        mat_glob = glob([base_dir, name, '\', labelo, '\', label, '_vol_img_*D.mat']);
        img = load(mat_glob{1});
        fdnames = fieldnames(img);
        img_3D = getfield(img, fdnames{1});
        
        
        heart_glob = glob(cat(2, base_dir, name, '/', labelo, '/Heart/mask*.mat'));
        heart = load(heart_glob{1});
        fdnames = fieldnames(heart);
        heart_3D = getfield(heart, fdnames{1});
        
        if numel(size(img_3D)) > 3
            % Then it is 4D (multiecho), Rows x Columns x TEs x Slices
            img_3D = squeeze(img_3D(:,:,1,:));
            heart_3D = squeeze(heart_3D(:,:,1,:));
        end
        
        [x, y, x_centroid, y_centroid] = GroovePickNCheckFunc3(img_3D, heart_3D, outputFileName, overwrite_label);
        
        coords.x = x; % Rows
        coords.y = y; % Column
        coords.x_centroid = x_centroid;
        coords.y_centroid = y_centroid;
        save(outputFileName, 'coords');
    end
end