clear all;
close all;

addpath('../function/')
fpath = uigetdir();

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
    'EchoTime',...
    'MagneticFieldStrength'
    };

f_glob = glob(cat(2, fpath, '/*/'));

for i = 1:length(f_glob)
    strings = strsplit(f_glob{i}, '/');
    fname = strings{end-1};

    [whatsinit, slice_data] = dicom23D(f_glob{i}, dicom_fields);

    output_struct = struct;
    output_struct.MagneticFieldStrength = slice_data(1).MagneticFieldStrength;

    eco = slice_data(1).EchoTime;
    eco_array = [];
    eco_array(1) = eco;
    n = 2;

    while eco ~= slice_data(n).EchoTime
        eco_array(n) = slice_data(n).EchoTime;
        n = n + 1;
    end

    output_struct.EchoTime = eco_array * 1e-3;
    temp_data = reshape(whatsinit, size(whatsinit,1), size(whatsinit, 2), ...
        length(eco_array), []);
    output_struct.image_data = permute(temp_data, [1,2,4,3]);
    
    save_f = cat(2, fpath, '/', fname, '.mat');
    save(save_f, '-struct', 'output_struct');
end

%% VIDA XA30
fpath = uigetdir();
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
    'EchoTime',...
    'MagneticFieldStrength'
    };

f_glob = glob(cat(2, fpath, '/*/'));

for i = 1:length(f_glob)
    
    strings = strsplit(f_glob{i}, '/');
    fname = strings{end-1};
    dir_glob = f_glob{i};
    dicom_glob = glob(cat(2, dir_glob, '*'));
    num_eco = length(dicom_glob);
    eco_array = zeros(1, num_eco);
    for j = 1:num_eco
        whatsinit = squeeze(dicomread(dicom_glob{j}));
        if j == 1
            whatsinit_4d = zeros([size(squeeze(whatsinit)), num_eco]);
        end
        whatsinit_4d(:,:,:,j) = whatsinit;
        info = dicominfo(dicom_glob{j});
        eco_array(j) = info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime;
    end

    output_struct = struct;
    output_struct.MagneticFieldStrength = info.MagneticFieldStrength;
    output_struct.EchoTime = eco_array * 1e-3;
    output_struct.image_data = whatsinit_4d;
    save_f = cat(2, fpath, '/', fname, '.mat');

    save(save_f, '-struct', 'output_struct');
end



