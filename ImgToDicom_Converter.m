clear all;
close all;
% Please refer to GUI2020_Update10112020 for writing Dicom

%% convert mat file to Dicom for Khalid
% load the data
data_dir = uigetdir;
load(cat(2, data_dir, '/VE_RES_DATA.mat'));

%% Try to plot it
addpath('function/');
figure(); imagesc(squeeze(D1{1}(1,:,:)));


save_dir = GetFullPath(cat(2, data_dir, '/../Dicom/'));
for d = 1:3
    for i = 1:length(D1)
        sub_dir = cat(2, save_dir, 'Subject_', num2str(i), '/');
        if ~exist(sub_dir, 'dir')
            mkdir(sub_dir);
        end
        D1_reorder = permute(D1{i}, [2,3,1]);
        D2_reorder = permute(D2{i}, [2,3,1]);
        D3_reorder = permute(D3{i}, [2,3,1]);
        
        for slc = 1:size(D1{i}, 1)
            D1_min = min(min(squeeze(D1_reorder(:,:,slc))));
            D1_max = max(max(squeeze(D1_reorder(:,:,slc))));
            D2_min = min(min(squeeze(D2_reorder(:,:,slc))));
            D2_max = max(max(squeeze(D2_reorder(:,:,slc))));
            D3_min = min(min(squeeze(D3_reorder(:,:,slc))));
            D3_max = max(max(squeeze(D3_reorder(:,:,slc))));
            
            D1_rescale = (D1_reorder(:,:,slc) - D1_min) ./ (D1_max - D1_min);
            D2_rescale = (D2_reorder(:,:,slc) - D2_min) ./ (D2_max - D2_min);
            D3_rescale = (D3_reorder(:,:,slc) - D3_min) ./ (D3_max - D3_min);
            
            fname = cat(2, 'D1_', num2str(i), '_slice', num2str(slc), '.dcm');
            dicomwrite(squeeze(D1_rescale), cat(2, sub_dir, fname));
            fname = cat(2, 'D2_', num2str(i), '_slice', num2str(slc), '.dcm');
            dicomwrite(squeeze(D2_rescale), cat(2, sub_dir, fname));
            fname = cat(2, 'D3_', num2str(i), '_slice', num2str(slc), '.dcm');
            dicomwrite(squeeze(D3_rescale), cat(2, sub_dir, fname));
        end
    end
end
