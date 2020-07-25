clear all;
close all;


base_dir = 'D:\Data\Cedars_2020_pig\DHARMAKUMAR_20P11_20P11_EXVIVO\RawData\';
fname = cat(2, base_dir, 'meas_MID00560_FID11379_3D_mGRE_sa_hiRES_0_9x0_9x1_Patient.dat');

twix_obj_in = mapVBVD(fname);