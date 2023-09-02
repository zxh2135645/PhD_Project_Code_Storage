clear all;
close all;
% addpath('/Users/jameszhang/Documents/Rohan/RawData/mapVBVD/');
base_dir = '/Users/jameszhang/Documents/Rohan/RawData/RawData_HeartPhantom_01312019/';
fname = cat(2, base_dir, 'meas_MID00132_FID44424_GRE_TI400_FA30.dat');

base_dir = 'D:\Data\Cedars_2020_pig\DHARMAKUMAR_20P11_20P11_EXVIVO\RawData\';
fname = cat(2, base_dir, 'meas_MID00560_FID11379_3D_mGRE_sa_hiRES_0_9x0_9x1_Patient.dat');
%[headers,protocol]=read_dat_headers_PK(fname)
%fid=fopen(fname,'r')

twix_obj_in = mapVBVD(fname); 
% what read in is k-space raw data, do I really want to do it?
