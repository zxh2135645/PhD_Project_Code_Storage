clear all;
close all;
%%
% addpath('/Users/jameszhang/Documents/Rohan/RawData/mapVBVD/');
base_dir = '/Users/jameszhang/Documents/Rohan/RawData/RawData_HeartPhantom_01312019/';
fname = cat(2, base_dir, 'meas_MID00132_FID44424_GRE_TI400_FA30.dat');

base_dir = 'D:\Data\Cedars_2020_pig\DHARMAKUMAR_20P11_20P11_EXVIVO\RawData\';
fname = cat(2, base_dir, 'meas_MID00560_FID11379_3D_mGRE_sa_hiRES_0_9x0_9x1_Patient.dat');

base_dir = '/Users/jameszhang/Downloads/';
%base_dir = '/Users/jameszhang/Documents/RYLab/Results/QSM_ucla/';
fname = cat(2, base_dir, 'meas_MID2518_CV_QSM_XZ_FID655969.dat');

fname = cat(2, base_dir, 'meas_MID177_MUSIC_freerun_2015_v3_FID658099.dat');
%fname = cat(2, base_dir, 'meas_MID172_CV_QSM_XZ2_FID658094.dat');
%[headers,protocol]=read_dat_headers_PK(fname)
%fid=fopen(fname,'r')


twix_obj_in = mapVBVD(fname); 
% what read in is k-space raw data, do I really want to do it?

%% See Sequence File Name
twix_obj_in.hdr.Config.SequenceFileName
%%
%figure;plot(twix_obj_in.image.Lin(9:8:end),'.')
figure;plot(twix_obj_in.RTfeedback.Lin(:),'.')

title(twix_obj_in.hdr.Config.SequenceFileName)
%% This code is for check the histogram for PE lines
%Ny = twix_obj_in.hdr.Meas.NLin;
%Nz = twix_obj_in.hdr.Meas.NImagePar;

Ny = twix_obj_in.hdr.Meas.iNoOfFourierLines
Nz = twix_obj_in.hdr.Meas.iNoOfFourierPartitions

figure;hist(twix_obj_in.image.Lin(9:16:end),Ny)

figure;hist(twix_obj_in.image.Par(9:16:end),Nz)


%%
Ny = twix_obj_in.hdr.Meas.iNoOfFourierLines

Nz = twix_obj_in.hdr.Meas.iNoOfFourierPartitions