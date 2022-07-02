clear all; 
close all;
%% 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');
% dispim = @(x)fftshift(x(:,:,SliceSlider.Value,:),1);
qMRinfo('mono_t2'); % set it up first in qMRLab
% T2star
% sizes(2) -> T1
% sizes(3) -> Cardiac
% sizes(4) -> resp
TE_array = [1.41, 3.38, 5.39, 7.40, 9.41, 11.42]; % 13.43, 15.44 ms
IR_array = repmat(TE_array, [sizes(2), 1]) + repmat((0:1:(sizes(2)-1)) * params.lEchoSpacing, [length(TE_array), 1]).'*1000;

mask_f = cat(2, fid_path, 'mask_rect.mat');
mask = zeros(Ny, Nx, Nz);
if ~exist(mask_f)
    for i = 1:Nz
        dispim = @(x,st)fftshift(x(:,:,i,:),1);
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
    save(mask_f, 'mask');
else
    load(mask_f);
end

%% fitting T2* maps and T1 maps (Skip to next section for SingleEcho Recon)
t1_map = zeros(Ny, Nx, Nz, 1, 1, 1);
t2star_map = zeros(Ny, Nx, Nz, 1, 1, sizes(2));
for i = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    temp = Gr\reshape(Phi(:,:,1,1,:), L, []);
    temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
    for j = 1:1
        for k = 1:1
            temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
            temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
            for l = 1:sizes(2)
                % Reshape matrix as [Width x Height x #Slice x #TE]
                ipt = abs(temp(:,:,l,:));
                Model = mono_t2;  % Create class from model
                %Model = Custom_OptionsGUI(Model);
                Model.Prot.SEdata.Mat = TE_array.'; %
                Model.st = [100 2000];
                Model.lb = [1 2000];
                Model.fx = [0 0];
                Model.voxelwise = 1;
                Model.options.FitType = 'Linear';
                data = struct;  % Create data structure
                data.SEdata = ipt;
                data.Mask = mask(:,:,i);
                FitResults = FitData(data, Model); %fit data
                t2star_map(:,:,i,j,k,l) = FitResults.T2;
            end
            
            for m = 1:1
                % a - create object
                Model = inversion_recovery;
                % Reshape matrix as [Width x Height x #Slice x #TE]
                ipt = abs(reshape(temp(:,:,:,m), Ny, Nx, 1, []));
                
                data = struct;
                data.IRData= double(ipt);
                data.Mask= double(mask(:,:,i));

                Model.Prot.IRData.Mat = IR_array(:,m);
                Model.voxelwise = 1;

                % b- fit dataset
                FitResults = FitData(data,Model,0);
                t1_map(:,:,i,j,k,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
            end
        end
    end
end

map_to_save = struct;
map_to_save.mask = mask;
map_to_save.t2star_map = t2star_map;
map_to_save.t1_map = t1_map;
save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings.mat');
save(save_f, 'map_to_save');

%% T2star and T1 altogether
% for i = 1:Nz
% for i = 3:3
%     dispim = @(x,st)fftshift(x(:,:,i,:),1);
%     for j = 1:1
%         for k = sizes:sizes(4)
%             for l = 1:sizes(2)
%                 % Reshape matrix as [Width x Height x #Slice x #TE]
%                 ipt = abs(temp(:,:,l,:));
%                 Model = mono_t2;  % Create class from model
%                 %Model = Custom_OptionsGUI(Model);
%                 Model.Prot.SEdata.Mat = TE_array.'; %
%                 Model.st = [100 2000];
%                 Model.lb = [1 2000];
%                 Model.fx = [0 0];
%                 Model.voxelwise = 1;
%                 Model.options.FitType = 'Linear';
%                 data = struct;  % Create data structure
%                 data.SEdata = ipt;
%                 data.Mask = mask(:,:,i);
%                 FitResults = FitData(data, Model); %fit data
%                 t2star_map(:,:,i,1,1,l) = FitResults.T2;
%             end
%             
%             for m = 1:sizes(5)
%                 % a - create object
%                 Model = inversion_recovery;
%                 % Reshape matrix as [Width x Height x #Slice x #TE]
%                 ipt = abs(reshape(temp(:,:,:,m), Ny, Nx, 1, []));
%                 
%                 data = struct;
%                 data.IRData= double(ipt);
%                 data.Mask= double(mask(:,:,i));
% 
%                 Model.Prot.IRData.Mat = IR_array(:,m);
%                 Model.voxelwise = 1;
% 
%                 % b- fit dataset
%                 FitResults = FitData(data,Model,0);
%                 t1_map(:,:,i,1,1,l) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
%             end
%         end
%     end
% end
% 
% map_to_save = struct;
% map_to_save.mask = mask;
% map_to_save.t2star_map = t2star_map;
% map_to_save.t1_map = t1_map;
% save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings.mat');
% save(save_f, 'map_to_save');
%% visualize slices
slc = 4;
dispim = @(x)fftshift(x(:,:,slc,:),1);

num_seg = 21;
figure();
for i = 1:size(Phi,3)
    temp = Gr\reshape(Phi(:,num_seg,i,1,:), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    cw = max(vec(abs(temp)));
    
    % figure();
    subplot(4,6,i);
    temp = imrotate(temp, 90);
    temp = abs(flip(temp(:,:,1),2)/cw); % First Section
    temp = uint16(temp*4095);
    ax2 = imagesc(temp); axis image; colormap gray;axis off;
    title(cat(2, 'Cardiac Phase ', num2str(i)));
end
%% For single-echo mapping
% T1 map only
t1_map = zeros(Ny, Nx, Nz, sizes(5));
%i = input(sprintf('Select Slice of Interest [%d]: ', 3));
%for i = 1:Nz
for i = 4:4
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    %for j = 1:1 % Cardiac
    for j = 11:11
        for k = 1:1 % Respiratory
            for m = 1:sizes(5)
                % a - create object
                Model = inversion_recovery;
                % Reshape matrix as [Width x Height x #Slice x #TE]
                
                temp = Gr\reshape(Phi(:,:,j,k,m), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                ipt = abs(reshape(temp(:,:,:), Ny, Nx, 1, []));
                
                data = struct;
                data.IRData= double(ipt);
                data.Mask= double(mask(:,:,i));
                
                Model.Prot.IRData.Mat = IR_array(:,1);
                Model.voxelwise = 1;
                
                % b- fit dataset
                FitResults = FitData(data,Model,0);
                
                mask_temp = double(mask(:,:,i));
                if any(mask_temp(:))
                    t1_map(:,:,i,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
                end
            end
        end
    end
end

map_to_save = struct;
map_to_save.mask = mask;
map_to_save.t1_map = t1_map;

% Plot
figure();
%for i = 1:sizes(5)
for i = 1:Nz
   subplot(3,5,i);
   imagesc(squeeze(t1_map(:,:,i,15))); axis image; axis off; caxis([200 900]);
   % title(cat(2, 'Section ', num2str(i)));
   title(cat(2, 'Slice ', num2str(i)));
end

% Load mask
%load(cat(2, fid_path, 'FID25768_21D35_D3_6mm_USR14%_L64_results_2021_07_28_23_41_Echo1_T1mapping_Analysis.mat'));
BW = zeros(size(t1_map,1), size(t1_map,2), size(t1_map,3));
for nt = 1:size(t1_map,3)
    dispim = @(x,st)fftshift(x(:,:,nt,:),1);
    temp = Gr\reshape(Phi(:,:,1,1,end), L, []);
    temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
    cw = 0.5*max(vec(abs(temp)));
    figure();
    imagesc(abs(temp(:,:,end))/cw); axis image;
    
    roi = drawpolygon;
    BW(:,:,nt) = createMask(roi);
end

t1_masked = BW.* t1_map;
t1_masked(isnan(t1_masked)) = 0;

N_nt = sizes(5);
% T1 values of ROI
t1_mean = zeros(N_nt,size(BW,3));
t1_sd = zeros(N_nt, size(BW,3));
for slc = 1:Nz
    BW_temp = BW(:,:,slc);
    if any(BW_temp(:))
        for i = 1:N_nt
            t1_mean(i,slc) = mean(nonzeros(t1_masked(:,:,slc,i)));
            t1_sd(i,slc) = std(nonzeros(t1_masked(:,:,slc,i)));
        end
    end
end

metrics.t1_mean = t1_mean;
metrics.t1_sd = t1_sd;
metrics.BW = BW;
save_fname = cat(2, fid_path, fid_file(1:end-4), '_T1mapping_MOLLI_Analysis.mat');
save(save_fname, '-struct', 'metrics');

save_f = cat(2, fid_path, fid_file(1:19), 'LRT_T1MOLLI_Mappings_Seg15.mat');
save(save_f, 'map_to_save', '-v7.3');

%% T2star mapping
NEco_old = params.NEco_old; % 6
% for i = 1:Nz
addpath('../function/');
echo_f_glob = glob(cat(2, fid_path, '*Echo???????.mat'));
N_seg = 15;
t2star_map = zeros(Ny, Nx, Nz, N_seg, sizes(5));

% i = input(sprintf('Select Slice of Interest [%d]: ', 3));
%
for slc = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,slc,:),1);
    mask_temp = mask(:,:,slc);
    if any(mask_temp(:))
        for nt = 1:sizes(5)
            for j = 1:1
                for k = 1:1
                    
                    for neco = 1:NEco_old
                        if neco == 1
                            temp_4D = zeros([size(temp), NEco_old]);
                        end
                        load(echo_f_glob{neco}, 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'params', 'sizes');
                        temp = Gr\reshape(Phi(:,:,j,k,nt), L, []);
                        temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                        %cw = 0.5*max(vec(abs(temp)));
                        %figure();
                        %implay(abs(temp)/cw); axis image; colormap gray;
                        temp_4D(:,:,:,neco) = temp;
                    end
                    
                    for l = 1:N_seg
                        % Reshape matrix as [Width x Height x #Slice x #TE]
                        ipt = abs(temp_4D(:,:,(sizes(2)-N_seg+l),:));
                        
                        Model = mono_t2;  % Create class from model
                        %Model = Custom_OptionsGUI(Model);
                        Model.Prot.SEdata.Mat = TE_array.'; %
                        Model.st = [100 2000];
                        Model.lb = [1 2000];
                        Model.fx = [0 0];
                        Model.voxelwise = 1;
                        Model.options.FitType = 'Linear';
                        data = struct;  % Create data structure
                        data.SEdata = ipt;
                        data.Mask = squeeze(mask(:,:,slc));
                        FitResults = FitData(data, Model); %fit data
                        t2star_map(:,:,slc,l,nt) = FitResults.T2;
                    end
                end
            end
        end
    end
end
% Plot
figure();
for i = 1:Nz
   subplot(3,5,i);
   %figure();
   imagesc(squeeze(t2star_map(:,:,i,1,end))); axis image; axis off; caxis([0 100]);
   title(cat(2, 'Section ', num2str(i)));
end



% Load mask
%load(cat(2, fid_path, 'FID25768_21D35_D3_6mm_USR14%_L64_results_2021_07_28_23_41_Echo1_T1mapping_Analysis.mat'));
t2star_masked = BW.*squeeze(t2star_map(:,:,:,1,:));
t2star_masked(isnan(t2star_masked)) = 0;
t2star_masked(t2star_masked < 0) = 0;
t2star_masked(t2star_masked > 100) = 100;
% T1 values of ROI
t2star_mean = zeros(size(t2star_map, 5),Nz);
t2star_sd = zeros(size(t2star_map, 5), Nz);
for slc = 1:Nz
    for i = 1:size(t2star_map, 5)
        t2star_mean(i,slc) = mean(nonzeros(t2star_masked(:,:,slc,i)));
        t2star_sd(i,slc) = std(nonzeros(t2star_masked(:,:,slc,i)));
    end
end

metrics.t2star_mean = t2star_mean;
metrics.t2star_sd = t2star_sd;
metrics.BW = BW;
save_fname = cat(2, fid_path, fid_file(1:end-4), '_T2starmapping_MOLLI_Analysis_Seg15.mat');
save(save_fname, '-struct', 'metrics');


% map_to_save.t2star_map = t2star_map;
% save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings_Seg15.mat');
% save(save_f, 'map_to_save', '-v7.3');



% Display image
t2star = squeeze(t2star_map(:,:,:,:,end));
t2star(t2star<0) = 0;
t2star(t2star>100) = 100;
t2star = mean(t2star, 4);
figure(); 
for i = 1:size(t2star,3)
    subplot(4,4,i);
    imagesc(t2star(:,:,i)); caxis([0 50]);
    axis image; colormap gray;
end

map_to_save.t2star = t2star;
save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings_Seg15.mat');
save(save_f, 'map_to_save', '-v7.3');


save_f = cat(2, fid_path, fid_file(1:15), 'LRT_T2starMappings_Seg15.mat');
save(save_f, 't2star', '-v7.3');


%% Single-echo mapping based on PSIR
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'psir_mat');
save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings_Seg15.mat');

qMRinfo('mono_t2');
Ny = size(psir_mat,1);
Nx = size(psir_mat,2);
Nseg = size(psir_mat,3);
Nsect = size(psir_mat,4);

t1_map = zeros(size(psir_mat,1),size(psir_mat,2),size(psir_mat,4));

%i = input(sprintf('Select Slice of Interest [%d]: ', 3));
% for m = 1:sizes(5)
% data = struct;
for m = sizes(5):sizes(5)
    % a - create object
    Model = inversion_recovery;
    % Reshape matrix as [Width x Height x #Slice x #TE]
    
    %temp = Gr\reshape(Phi(:,:,j,k,m), L, []);
    %temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
    ipt = abs(reshape(psir_mat(:,:,:,m), Ny, Nx, 1, []));
    
    data.IRData= double(ipt);
    
    if m == 1
        figure();
        imagesc(ipt(:,:,1,end));
        roi = drawpolygon;
        mask = createMask(roi);
        data.Mask= double(mask);
    end
    
    Model.Prot.IRData.Mat = IR_array(:,1);
    Model.voxelwise = 1;
    
    % b- fit dataset
    FitResults = FitData(data,Model,0);
    
    mask_temp = double(mask);
    if any(mask_temp(:))
        t1_map(:,:,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
    end
end

map_to_save = struct;
map_to_save.mask = mask;
map_to_save.t1_map = t1_map;

%% ImageJ
% Convert mat to dicom (T1 map)
load(strcat(fid_path, fid_file(1:15), 'LRT_Mappings_Seg15_CardiacPhase11.mat'), 'map_to_save');
t1_map = map_to_save.t1_map;

figure();
for i = 1:size(t1_map, 3)
    subplot(4,4,i)
    imagesc(t1_map(:,:,i,1).*mask(:,:,i)); caxis([100 1000]); axis image;
end

save_dir = cat(2, fid_path, 'DICOM_T1/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));
label = 'LRT';
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

slc = 3;
for m = 1:sizes(5)
    %metadata = dicominfo(slice_data{1}(m).Filename);
    metadata = dicominfo(slice_data(m).Filename);
    
    t1 = uint16(t1_map(:,:,slc,m));
    metadata.WindowCenter = 500;
    metadata.WindowWidth =  800;
    metadata.SmallestImagePixelValue = min(t1(:));
    metadata.LargestImagePixelValue = max(t1(:));
    
    fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_Section', num2str(m), '_ForImageJ', '.dcm');
    dicomwrite(t1, fname, metadata, 'CreateMode', 'copy');
end

