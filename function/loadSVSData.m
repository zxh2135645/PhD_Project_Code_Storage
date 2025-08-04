function ss = loadSVSData(dst_dir, SeriesNum, InstanceNum)
    

    dt = Spectro.dicomTree('dir',dst_dir);

    matched = dt.searchForSeriesInstanceNumber(SeriesNum, InstanceNum);

    ss = Spectro.Spec(matched);

    % ss = Spectro.dicom(matched);
end