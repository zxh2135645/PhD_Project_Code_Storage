if (size(RawData_All_3D_allcoil,3)>12 )
    for necho=1
       RawData_All_Recon = squeeze(RawData_All_3D_allcoil(:,:,:,:,necho));
   NumProj_comb=size(RawData_All_Recon,2);
    %Calib_proj_indices=RawData_index(1:NumProj);
    tic
    patitiontocount=1:NumPartitions;
parfor patcount=patitiontocount
%if (NumDim==3)
   patcount
  % RawData_All(:,:,:)=RawData_All_Recon(:,:,:,patcount); 
%end
%RawData_All_Calib = RawData_All;
img_grid_temp(:,:,:,patcount)=GratioRecon_only (RawData_All_Recon(:,:,:,patcount),NumProj_comb,BaseResol,[1:NumProj_comb],RawData_index,NumChannels,gridOS_ndc,kernelWidth_ndc,beta,AzimuthalAngle,flag_calibDoGaussApod,Calib_apod_FWHM_ratio2BaseResol);
    %GratioRecon_inline_CoilEstimation_Jan2013;
   %CalibRegridRadial_thisSlice_3D_nonbin(:,:,patcount)=SOS(CalibRegridRadial_thisSlice);
end
img_grid=permute(img_grid_temp,[1,2,4,3]);
clear img_grid_temp;
toc
N=BaseResol;
for patcount=patitiontocount
    img_grid(:,:,patcount,:) = img_grid(:,:,patcount,:)./repmat(gridKernel_ndc,[1 1 1 NumChannels]);
end
    img_grid_combine_si = sqrt(sum(abs(img_grid).^2,4));
    slicer(shiftdim(abs(img_grid_combine_si),2),shiftdim(abs(img_grid_combine_si),3))
    % ---Channel compression
    sroi = img_grid(N/4+1:N/4+N/2,N/4+1:N/4+N/2,1:NumPartitions,:);
    sroi = reshape(sroi,[N*N*NumPartitions/4 NumChannels]);
    clear img_grid
    P = zeros(NumChannels,NumChannels);
    NN = length(sroi);
    norms = sum(abs(sroi).^2,2);
    for ii = 1:NumChannels
        for jj = 1:NumChannels
            P(ii,jj) = sum(sroi(:,ii).*conj(sroi(:,jj))./norms,1); % w/ normalization
        end
    end
    [UU,SS,VV] = svd(P);
    channels_red = 12;
    VV = VV';
    A = VV(1:channels_red,:);
    % ---Channel compression
    Raw_data_red = zeros(NumProj_comb,channels_red,N);
    Raw_data=permute(RawData_All_3D_allcoil,[2 3 1 4 5]);
    for ee=1:NumContrasts% apply same matrix to all echo
    for kk=1:NumPartitions
    for ii = 1:NumProj_comb
        for jj = 1:BaseResol
            Raw_data_red(ii,:,jj,kk,ee) = A*squeeze(Raw_data(ii,:,jj,kk,ee).');
        end
    end
    end
    end
    RawData_All_3D(:,:,:,:,:) = permute(Raw_data_red,[3 1 2 4 5]);
    %
    Raw_data_red = zeros(NumProj_comb,channels_red,N);
    Raw_data=permute(RawData_All_3D_Org,[2 3 1 4 5]);
    for ee=1:NumContrasts% apply same matrix to all echo
    for kk=1:NumPartitions
    for ii = 1:NumProj_comb
        for jj = 1:BaseResol
            Raw_data_red(ii,:,jj,kk,ee) = A*squeeze(Raw_data(ii,:,jj,kk,ee).');
        end
    end
    end
    end
    RawData_All_3D_Org = permute(Raw_data_red,[3 1 2 4 5]);
    %
    %NumChannels = channels_red;
  %  clear sroi Raw_data Raw_data_red
  %AR(:,:,necho)=A;
    end
   NumChannels=size(RawData_All_3D,3);
else
    RawData_All_3D=RawData_All_3D_allcoil;
end   
    clear RawData_All_3D_allcoil Raw_data_red;