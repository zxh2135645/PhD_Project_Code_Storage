function [volume_image, slice_data] = convert_vida_header(dicom_f, info)
slice_data = info;
if ~isfield(info, 'SliceLocation')
    slice_data.SliceLocation = info.PerFrameFunctionalGroupsSequence.Item_1.Private_0021_11fe.Item_1.Private_0021_1188;
end
volume_image = double(dicomread(dicom_f));

end