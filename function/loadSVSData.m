function ss = loadSVSData(dst_dir, SeriesNum, InstanceNum)
    addpath('../../OXSA/');

    dt = Spectro.dicomTree('dir',dst_dir);

    matched = dt.searchForSeriesInstanceNumber(SeriesNum, InstanceNum);

    ss = Spectro.Spec(matched);
end