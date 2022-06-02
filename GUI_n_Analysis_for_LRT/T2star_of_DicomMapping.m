clear all;
close all;
%% generate T2* maps of single-slice mGRE sequences of D8-P1 or D6-P2
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

labels = {'T2_MULTIECHO'};
label = labels{1};
qMRinfo('mono_t2'); 

idx_array = contains(folder_glob, label);
[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    'TriggerTime',...
    'RepetitionTime',...
    'EchoTime', 
    };

whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end
%%
t2star_w_4D = zeros(size(whatsinit{1},1), size(whatsinit{1},2), size(whatsinit{1},3), length(whatsinit));
for i = 1:size(t2star_w_4D,3)
    t2star_w_4D(:,:,:,i) = whatsinit{i};
end
t2star_w_4D = permute(t2star_w_4D, [1,2,4,3]);

TE_array = zeros(size(whatsinit{1},3),1);
for i = 1:length(TE_array)
   TE_array(i) = slice_data(i).EchoTime; 
end

Model = mono_t2;  % Create class from model
Model.Prot.SEdata.Mat = TE_array; %
Model.st = [100 2000];
Model.lb = [1 2000];
Model.fx = [0 0];
Model.voxelwise = 1;
Model.options.FitType = 'Linear';
data = struct;  % Create data structure
data.SEdata = t2star_w_4D;
%data.Mask = ones(size(whatsinit{1},1), size(whatsinit{1},2), length(whatsinit));
FitResults = FitData(data, Model, 0); %fit data


fh = figure('Position', [100 100 600 800]);
axis tight manual
for i = 1:length(whatsinit)
    subplot(3,4,i); imagesc(FitResults.T2(:,:,i)); axis image; caxis([0 100]);
end
%%
figure();
imagesc(FitResults.T2(:,:,1)); axis image; caxis([0 100]);
roi = drawpolygon;
remote = createMask(roi);

figure();
imagesc(FitResults.T2(:,:,1)); axis image; caxis([0 100]);
roi = drawpolygon;
hemo = createMask(roi);

mean_t2_remote_array = zeros(length(whatsinit),1);
mean_t2_hemo_array = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    mean_t2_remote_array(i) = mean(nonzeros(FitResults.T2(:,:,i) .* remote));
    mean_t2_hemo_array(i) = mean(nonzeros(FitResults.T2(:,:,i) .* hemo));
end
%%
figure();
x = [1,4,6,8,10,12,14,16,NaN,NaN,NaN,NaN];
plot(x, mean_t2_remote_array, 'LineWidth', 1.5); hold on;
plot(x, mean_t2_hemo_array, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('T2star (ms)');
grid on;
legend({'Remote', 'Hemorrhage'});
