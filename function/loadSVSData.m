function ss = loadSVSData(dst_dir, SeriesNum, InstanceNum)
    oxsa_dir = 'D:\OXSA\';
    addpath(oxsa_dir);

    
    dt = Spectro.dicomTree('dir',dst_dir);

    matched = dt.searchForSeriesInstanceNumber(SeriesNum, InstanceNum);

    ss = Spectro.Spec(matched);
end