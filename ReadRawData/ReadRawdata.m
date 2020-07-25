%%Read Siemens .dat Rawdata- Randy Yang 2017 Jan
%%Fnction to read Rawdata with 2-4Dimensions
%%Run Rawdata=ReadRawdata(cd,'meas_MID29_mGRE_3D_TR_TE_40_1echoes_dizon_noFS_TE3_3_FID27207.dat');
%%Rawdata=>[NumProj,BaseResol,NumPartitions,NumChannels] 
%%Confidential; Do not distribute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set path and flags
% path_datfile=['/Users/XGuan/Documents/MATLAB/ReadRawdata/Twix/insaline'];
% 
% filename_datfile=['meas_MID49_mGRE_3D_TR_TE_40_1echoes_dizon_noFS_TE3_3_FID27227.dat'];
% cd(path_datafile)

function RawData_All_3D_Org= ReadRawdata(path_datfile,filename_datfile)

%path_datfile=['/Users/XGuan/Documents/MATLAB/ReadRawdata/Twix/insaline'];

%filename_datfile=['meas_MID49_mGRE_3D_TR_TE_40_1echoes_dizon_noFS_TE3_3_FID27227.dat'];

flag_saveMatFile = 1;
flag_loadsavedMat = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Reading Raw Data --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~flag_loadsavedMat)
    [headers,protocol] = read_dat_headers_PK([filename_datfile]);
    disp(['Read Phantom Data File: ' filename_datfile]); disp(' ' );
    disp(['BaseResolution = ' num2str(protocol.sKSpace.lBaseResolution)]);
    disp([    ', #ReadoutSamples = ' num2str(protocol.sKSpace.lBaseResolution*2)]);
    disp([    ', #RadialViews = ' num2str(protocol.sKSpace.lRadialViews)]);
    disp([    ', #coils = ' num2str(length(protocol.asCoilSelectMeas{1}.asList))]);
    disp([    ', #segments = ' num2str(protocol.sFastImaging.lSegments)]);
    disp([    ', #shots = ' num2str(protocol.sFastImaging.lShots)]);
    disp([    ', #slices = ' num2str(length(protocol.sSliceArray.asSlice))]);
    TotalScanTime = protocol.lTotalScanTimeSec;
    disp([    ', Total Scan Time = ' num2str(TotalScanTime)]);
    disp([    ', Receiver dwelll time = ' num2str(protocol.sRXSPEC.alDwellTime{1}/1000)]);
    disp([    ', Readout FOV = ' num2str(protocol.sSliceArray.asSlice{1}.dReadoutFOV)]);
    disp([    ', Slice thickness = ' num2str(protocol.sSliceArray.asSlice{1}.dThickness)]);
    disp([    ', Flip Angle = ' num2str(protocol.adFlipAngleDegree{1})]);
    TRperRO = protocol.alTR{1}/1000 / protocol.sKSpace.lRadialViews;
    disp([    ', TR per r/o = ' num2str(TRperRO)]);
    disp([    ', TE = ' num2str(protocol.alTE{1}/1000)]);
    dwell_time = protocol.sRXSPEC.alDwellTime{1}/1000;
    rBW = round(1/(dwell_time*1e-6*2)/protocol.sKSpace.lBaseResolution);
    disp([    ', rBW = ' num2str(rBW)]);  disp(' ');
    %%%---MY WIP PARAMETER--%%%
%    disp('UI Gradient Delay Correction Parameters: '); disp([protocol.sWiPMemBlock.adFree ]);
%    disp('UI drop-down menu Parameter: '); disp([protocol.sWiPMemBlock.alFree]); disp(' ');
%    R_inshot = protocol.sWiPMemBlock.alFree{2};
%    disp([' *** R_inshot = ' num2str(R_inshot) ' ***']);
    BaseResol = protocol.sKSpace.lBaseResolution*2;
    NumChannels = length(protocol.asCoilSelectMeas{1}.asList);
    interleavetrigger=1;
    NumProj = protocol.sKSpace.lPhaseEncodingLines*protocol.sPhysioImaging.lPhases*interleavetrigger;
    NumShots = protocol.sFastImaging.lShots; % # of shots in the scan prescription (# heart-beats needed to acq. the full # of projections) -- determines size of sliding window
    NumSegments = protocol.sFastImaging.lSegments; % # of segments -- projections per time frame
    NumSlices = length(protocol.sSliceArray.asSlice);
    NumPartitions = protocol.sKSpace.lPartitions;
    NumEcho=protocol.lContrasts;
    NumRepetitions = 1;
%     if (protocol.lContrasts>1)
%     NumContrasts=protocol.lContrasts.*(size(list,1)-2);
%     end
if (protocol.lContrasts>1)
   IS_mGRE=1;
    end
    disp(['Reading Phantom Data File: ' filename_datfile]); disp(' ' );

    NumDim=str2num(protocol.sKSpace.ucDimension(4))/2+1;%2 is 2D  4 is 3D
    RawReadFlags.Dimension = NumDim; 
   
    RawReadFlags.Version = 21;
    Tarj=str2num(protocol.sKSpace.ucTrajectory(4));%2 is radial 1 is cartisian
    RawReadFlags.Trajectory = Tarj; 
    RawReadFlags.NumProj = NumProj; RawReadFlags.NumNoiseAcqs = 3; % ### needed to work with ReadRaw_GratioBS ###
    RawReadFlags.NumSegments=NumSegments;
    RawReadFlags.BaseResol=BaseResol;
    RawReadFlags.NumEcho=NumEcho;
    RawReadFlags.Presto=0;
    if (RawReadFlags.Dimension==2)
    RawReadFlags.NumPartitions=1;
    else
            RawReadFlags.NumPartitions=NumPartitions;
    end
    RawReadFlags.NumChannels=NumChannels;
%    tmpRaw = ReadRaw_GratioBS([path_datfile filename_datfile '.dat'],RawReadFlags,0); %@@@
    if (~exist([path_datfile,'RawData.mat']))
 %if (1)
   % [tmpRaw, NoiseData, RefData, RawReadFlags,RawData_index,SG_Im_data,SG_Im_data_nonft] = Read_RawData_BS_presto([filename_datfile],RawReadFlags,0); %@@@
   % [tmpRaw, NoiseData, RefData, RawReadFlags,RawData_index,SG_Im_data,SG_Im_data_nonft] = ReadData_presto_mEchoPhase([filename_datfile],RawReadFlags,0); %@@@
%     [tmpRaw, NoiseData, RefData, RawReadFlags,RawData_index_tmp,SG_Im_data_tmp,SG_Im_data_nonft] = ReadData_presto_mEchoPhaseNAV([filename_datfile],RawReadFlags,1); %@@@
     [tmpRaw, NoiseData, RefData, RawReadFlags,RawData_index_tmp] = ReadData_presto_mEchoPhaseNAV([filename_datfile],RawReadFlags,1); %@@@
    %flag_doROreflect=0;
   %FileName= [filename_datfile];
%    Read_RawData_BS_presto
%    if (RawReadFlags.Dimension==3 )
filenum=1;
        tmpRaw2D=tmpRaw(:,:,3,:,1);
        tmpRaw=tmpRaw( (RawData_index_tmp>0),:,:,:,:);%preserve lines with data
        RawData_index_tmp=RawData_index_tmp((RawData_index_tmp>0)); %preserve lines with data
        RawData_All_3D_Org(:,:,:,:,:,filenum)=permute(tmpRaw,[2,1,4,3,5,6]);
         if (RawReadFlags.Dimension==3 )
         RawData_All_3D(:,:,:,:,:,filenum)=permute(squeeze(fftshift(ifft(ifftshift(tmpRaw,3),[],3),3)), [2 1 4 3 5,6]);
         else
           RawData_All_3D(:,:,:,:,:,filenum)= tmpRaw;
         end
%         SG_Im_data(:,:,:,(filenum-1)*size(SG_Im_data_tmp,4)+1:filenum*size(SG_Im_data_tmp,4))=SG_Im_data_tmp;
%         RawData_index(:,:,(filenum-1)*size(SG_Im_data_tmp,4)+1:filenum*size(SG_Im_data_tmp,4))=repmat(RawData_index_tmp,[1,1,size(SG_Im_data_tmp,4)]);
        RawData_index=repmat(RawData_index_tmp,[1,1,1]);
        %save([list(filenum).name(1:end-4),'RawData.mat'],'RawData_All_3D_Org','RawData_All_3D_Org','RawData_index','SG_Im_data');
            RawData_All_3D_Org=squeeze(RawData_All_3D_Org);
            RawData_All_3D=squeeze(RawData_All_3D);
            RawData_index=squeeze(RawData_index);
%            SG_Im_data=squeeze(SG_Im_data);
        save([path_datfile,'\RawData.mat'],'RawData_All_3D_Org','RawData_All_3D','RawData_index');%disable
        %save data
%    else
       % tmpRaw2D=tmpRaw;
%    end

    RawData_All = permute(squeeze(tmpRaw2D), [2 1 3]); % of size: (BaseResol, NumProj, NumChannels)
    clearvars tmpRaw NoiseData RefData RawReadFlags tmpRaw2D;
        else
        load([path_datfile,'\RawData.mat'])
    end
else
    if (~flag_doAltMATfileName)
        savedMat_filename_list = dir([SaveMatPATH filesep 'GoldenRatio_' patient_id_tag_file '*' 'Nprj2Recon' num2str(Nprj2Recon) '_SW' num2str(SWwidth) '*.mat' ]);
    else
        savedMat_filename_list = dir([SaveMatPATH filesep 'GoldenRatio_' dataset_index '*' 'Nprj2Recon' num2str(Nprj2Recon) '_SW' num2str(SWwidth) '*.mat' ]);
    end
    savedMat_filename = savedMat_filename_list(1).name
    load([SaveMatPATH filesep savedMat_filename]);
end


