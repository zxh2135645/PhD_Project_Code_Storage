clear all;
close all;

base_dir = uigetdir;

bl_excel = load(cat(2, base_dir, '/Results/Metrics_FatNIron_Analysis_BL_excel.mat'));
bl_good = load(cat(2, base_dir, '/Results/Metrics_FatNIron_Analysis_BL_good.mat'));
fu_excel = load(cat(2, base_dir, '/Results/Metrics_FatNIron_Analysis_excel.mat'));
fu_good = load(cat(2, base_dir, '/Results/Metrics_FatNIron_Analysis_good.mat'));

ff_array_bl_excel_p = bl_excel.Chord50_BL.mean_ff_array_50chord_hemo_p;
r2star_array_bl_excel_p = bl_excel.Chord50_BL.mean_r2star_array_50chord_hemo_p;
ff_array_bl_excel_n = bl_excel.Chord50_BL.mean_ff_array_50chord_hemo_n;
r2star_array_bl_excel_n = bl_excel.Chord50_BL.mean_r2star_array_50chord_hemo_n;

ff_array_bl_good_p = bl_good.Chord50_BL.mean_ff_array_50chord_hemo_p;
r2star_array_bl_good_p = bl_good.Chord50_BL.mean_r2star_array_50chord_hemo_p;
ff_array_bl_good_n = bl_good.Chord50_BL.mean_ff_array_50chord_hemo_n;
r2star_array_bl_good_n = bl_good.Chord50_BL.mean_r2star_array_50chord_hemo_n;

ff_array_fu_excel_p = fu_excel.Chord50.mean_ff_array_50chord_hemo_p;
r2star_array_fu_excel_p = fu_excel.Chord50.mean_r2star_array_50chord_hemo_p;
ff_array_fu_excel_n = fu_excel.Chord50.mean_ff_array_50chord_hemo_n;
r2star_array_fu_excel_n = fu_excel.Chord50.mean_r2star_array_50chord_hemo_n;

ff_array_fu_good_p = fu_good.Chord50.mean_ff_array_50chord_hemo_p;
r2star_array_fu_good_p = fu_good.Chord50.mean_r2star_array_50chord_hemo_p;
ff_array_fu_good_n = fu_good.Chord50.mean_ff_array_50chord_hemo_n;
r2star_array_fu_good_n = fu_good.Chord50.mean_r2star_array_50chord_hemo_n;


r2star_array_bl_excel_p(isnan(ff_array_bl_excel_p)) = [];
r2star_array_bl_good_p(isnan(ff_array_bl_good_p)) = [];
r2star_array_bl_excel_n(isnan(ff_array_bl_excel_n)) = [];
r2star_array_bl_good_n(isnan(ff_array_bl_good_n)) = [];
ff_array_bl_excel_p(isnan(ff_array_bl_excel_p)) = [];
ff_array_bl_good_p(isnan(ff_array_bl_good_p)) = [];
ff_array_bl_excel_n(isnan(ff_array_bl_excel_n)) = [];
ff_array_bl_good_n(isnan(ff_array_bl_good_n)) = [];

r2star_array_fu_excel_p(isnan(ff_array_fu_excel_p)) = [];
r2star_array_fu_good_p(isnan(ff_array_fu_good_p)) = [];
r2star_array_fu_excel_n(isnan(ff_array_fu_excel_n)) = [];
r2star_array_fu_good_n(isnan(ff_array_fu_good_n)) = [];
ff_array_fu_excel_p(isnan(ff_array_fu_excel_p)) = [];
ff_array_fu_good_p(isnan(ff_array_fu_good_p)) = [];
ff_array_fu_excel_n(isnan(ff_array_fu_excel_n)) = [];
ff_array_fu_good_n(isnan(ff_array_fu_good_n)) = [];

ff_array_bl_p = [ff_array_bl_excel_p, ff_array_bl_good_p].';
ff_array_bl_n = [ff_array_bl_excel_n, ff_array_bl_good_n].';
r2star_array_bl_p = [r2star_array_bl_excel_p, r2star_array_bl_good_p].';
r2star_array_bl_n = [r2star_array_bl_excel_n, r2star_array_bl_good_n].';

ff_array_fu_p = [ff_array_fu_excel_p, ff_array_fu_good_p].';
ff_array_fu_n = [ff_array_fu_excel_n, ff_array_fu_good_n].';
r2star_array_fu_p = [r2star_array_fu_excel_p, r2star_array_fu_good_p].';
r2star_array_fu_n = [r2star_array_fu_excel_n, r2star_array_fu_good_n].';

ff_array_bl_excel_p = ff_array_bl_excel_p.';
ff_array_bl_excel_n = ff_array_bl_excel_n.';
r2star_array_bl_excel_p = r2star_array_bl_excel_p.';
r2star_array_bl_excel_n = r2star_array_bl_excel_n.';
ff_array_fu_excel_p = ff_array_fu_excel_p.';
ff_array_fu_excel_n = ff_array_fu_excel_n.';
r2star_array_fu_excel_p = r2star_array_fu_excel_p.';
r2star_array_fu_excel_n = r2star_array_fu_excel_n.';

%% ROI

ff_array_bl_excel_p = bl_excel.ROI_BL.mean_ff_roi_array_hemo_p;
r2star_array_bl_excel_p = bl_excel.ROI_BL.mean_r2star_roi_array_hemo_p;
ff_array_bl_excel_n = bl_excel.ROI_BL.mean_ff_roi_array_hemo_n;
r2star_array_bl_excel_n = bl_excel.ROI_BL.mean_r2star_roi_array_hemo_n;

ff_array_bl_good_p = bl_good.ROI_BL.mean_ff_roi_array_hemo_p;
r2star_array_bl_good_p = bl_good.ROI_BL.mean_r2star_roi_array_hemo_p;
ff_array_bl_good_n = bl_good.ROI_BL.mean_ff_roi_array_hemo_n;
r2star_array_bl_good_n = bl_good.ROI_BL.mean_r2star_roi_array_hemo_n;

ff_array_fu_excel_p = fu_excel.ROI.mean_ff_roi_array_hemo_p;
r2star_array_fu_excel_p = fu_excel.ROI.mean_r2star_roi_array_hemo_p;
ff_array_fu_excel_n = fu_excel.ROI.mean_ff_roi_array_hemo_n;
r2star_array_fu_excel_n = fu_excel.ROI.mean_r2star_roi_array_hemo_n;

ff_array_fu_good_p = fu_good.ROI.mean_ff_roi_array_hemo_p;
r2star_array_fu_good_p = fu_good.ROI.mean_r2star_roi_array_hemo_p;
ff_array_fu_good_n = fu_good.ROI.mean_ff_roi_array_hemo_n;
r2star_array_fu_good_n = fu_good.ROI.mean_r2star_roi_array_hemo_n;

ff_array_bl_excel_p(isnan(ff_array_bl_excel_p)) = [];
ff_array_bl_good_p(isnan(ff_array_bl_good_p)) = [];
ff_array_bl_excel_n(isnan(ff_array_bl_excel_n)) = [];
ff_array_bl_good_n(isnan(ff_array_bl_good_n)) = [];
r2star_array_bl_excel_p(isnan(r2star_array_bl_excel_p)) = [];
r2star_array_bl_good_p(isnan(r2star_array_bl_good_p)) = [];
r2star_array_bl_excel_n(isnan(r2star_array_bl_excel_n)) = [];
r2star_array_bl_good_n(isnan(r2star_array_bl_good_n)) = [];

ff_array_fu_excel_p(isnan(ff_array_fu_excel_p)) = [];
ff_array_fu_good_p(isnan(ff_array_fu_good_p)) = [];
ff_array_fu_excel_n(isnan(ff_array_fu_excel_n)) = [];
ff_array_fu_good_n(isnan(ff_array_fu_good_n)) = [];
r2star_array_fu_excel_p(isnan(r2star_array_fu_excel_p)) = [];
r2star_array_fu_good_p(isnan(r2star_array_fu_good_p)) = [];
r2star_array_fu_excel_n(isnan(r2star_array_fu_excel_n)) = [];
r2star_array_fu_good_n(isnan(r2star_array_fu_good_n)) = [];

ff_array_bl_p = [ff_array_bl_excel_p, ff_array_bl_good_p].';
ff_array_bl_n = [ff_array_bl_excel_n, ff_array_bl_good_n].';
r2star_array_bl_p = [r2star_array_bl_excel_p, r2star_array_bl_good_p].';
r2star_array_bl_n = [r2star_array_bl_excel_n, r2star_array_bl_good_n].';

ff_array_fu_p = [ff_array_fu_excel_p, ff_array_fu_good_p].';
ff_array_fu_n = [ff_array_fu_excel_n, ff_array_fu_good_n].';
r2star_array_fu_p = [r2star_array_fu_excel_p, r2star_array_fu_good_p].';
r2star_array_fu_n = [r2star_array_fu_excel_n, r2star_array_fu_good_n].';

ff_array_bl_excel_p = ff_array_bl_excel_p.';
ff_array_bl_excel_n = ff_array_bl_excel_n.';
r2star_array_bl_excel_p = r2star_array_bl_excel_p.';
r2star_array_bl_excel_n = r2star_array_bl_excel_n.';
ff_array_fu_excel_p = ff_array_fu_excel_p.';
ff_array_fu_excel_n = ff_array_fu_excel_n.';
r2star_array_fu_excel_p = r2star_array_fu_excel_p.';
r2star_array_fu_excel_n = r2star_array_fu_excel_n.';
