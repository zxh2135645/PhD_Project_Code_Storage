clear all;
close all;
%%
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');

TE_array = [1.41, 3.38, 5.39, 7.40, 9.41, 11.42]; % 13.43, 15.44 ms
TR = params.TR / 1e3 / 192;

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

%%
TR = params.lEchoSpacing * 1e3;
% TR = params.TR / 1e3 / 192;
alpha = params.adFlipAngleDegree * pi/180;
t1_map = zeros(Ny, Nx, Nz, sizes(5));
%%
t1_map_array = zeros(Ny*Nx,1);
for i = 3:3
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    mask_2d = mask(:,:,i);
    for j = 1:1
        for k = 1:1
            tic;
            %for m = sizes(5):sizes(5)
            for m = 1:sizes(5)
                % Reshape matrix as [Width x Height x #Slice x #TE]
                temp = Gr\reshape(Phi(:,:,j,k,m), L, []);
                temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                ipt = abs(reshape(temp(:,:,:), [], 192));
                mask_1d = mask_2d(:);
                
                n = (21:192)-20;

                
                lb = [100, -1, 0];
                ub = [2000, -0.5, 1];
                x0 = [100, -1, 1];

                for kk = 1:length(mask_1d)
                    if mask_1d(kk) == 1
                        ipt_kk = ipt(kk, :);
                        [~,b] = min(ipt_kk);
                        ipt_kk(1:b-1) = -ipt_kk(1:b-1);

                        y = ipt_kk(21:192);
                        fun = @(x) x(3) .* (1 - exp(-TR/x(1)))./(1 - cos(alpha).*exp(-TR/x(1))) .* ...
                            (1 + (x(2)-1).*(exp(-TR/x(1)).*cos(alpha)).^(n-1)) - y;
                        x = lsqnonlin(fun,x0,lb,ub);
                        t1_map_array(kk) = x(1);
                    end
                end
                t1_map(:,:,i,m) = reshape(t1_map_array, Ny, Nx);
            end
            toc;
            
        end
    end
end

t1_map_test = reshape(t1_map_array, Ny, Nx);
%%
addpath('../function/');
figure();
for i = 1:size(t1_map, 3)
    subplot(4,4,i)
    imagesc(t1_map(:,:,i,1).*mask(:,:,i)); caxis([100 1000]); axis image;
end

save_dir = cat(2, fid_path, 'DICOM_T1_Analytical/');
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
    'EchoTime',...
    'EchoTrainLength',
    };

whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end

NumEcho = slice_data(1).EchoTrainLength;
NumSlc = length(slice_data) / NumEcho;
slc = 3;

for m = 1:sizes(5)
    %metadata = dicominfo(slice_data{1}(m).Filename);
    
    i = ceil(m/NumSlc);
    slc_virtual = m - (i-1) * NumSlc;
    
    idx = (slc_virtual - 1)*NumEcho + i;
    metadata = dicominfo(slice_data(idx).Filename);
    
    t1 = uint16(t1_map(:,:,slc,m));
    metadata.WindowCenter = 500;
    metadata.WindowWidth =  800;
    metadata.SmallestImagePixelValue = min(t1(:));
    metadata.LargestImagePixelValue = max(t1(:));
    
    fname = cat(2, save_dir, fid_file(1:19), 'Echo', num2str(1), '_Section', num2str(m), '_Slice', num2str(slc), '_ForImageJ', '.dcm');
    dicomwrite(t1, fname, metadata, 'CreateMode', 'copy');
end

%%
% X(1) -> T1
% x(2) -> B 
% x(3) -> C