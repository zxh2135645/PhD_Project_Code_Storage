function Func_ConvertMat_Dicom_ImageJ(t1_map, fid_file, dicom_dir, save_dir, slc, sizes)
%% ImageJ
% Convert mat to dicom (T1 map) for Image J
% t1_map should be 4D: Ny * Nx * Nz * NSeg
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
end