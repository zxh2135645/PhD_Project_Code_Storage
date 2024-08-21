clear all;
close all;
% XZ 11/23/2023
addpath(genpath(pwd));
rmpath(genpath('./Example/'));
addpath(genpath('../reconstructionPipeline/'));
addpath('../function/');
FlagMask=1;
FlagSaveB0=1;
FlagLoadFolder=1;
IsManual_Mask=1;
Scan_Type='Cardiac';

dicom_dir = uigetdir();
resdir=GetFullPath([dicom_dir, '/Result/']);
if ~exist(resdir)
    mkdir(resdir);
end

Names = {'DMR0466920', 'DMR0475552', 'DMR0481292', 'DMR0489758', 'DMR0493560', 'DMR0506800', 'DMR0540385'};
time_points = {'FU'};
label = 'T2star';
dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    'EchoTime',...
    'MagneticFieldStrength',...
    'ImagingFrequency',...
    'ImageType','RescaleSlope','RescaleIntercept',...
    'EchoNumber',...
    };

%% Main body
% dicom_folder = '3D_IDEAL_Avg4_Bipolar';
% dicom_folder = 'T2STARMAP_ANATOMICAL_0014';              % 21P18
%dicom_folder = '3D_MGRE_SA_HIRES_0_9X0_9X1_PATIENT_0128'; % 21P18
%dicom_folder = '3D_MGRE_SA_HIRES_0_9X0_9X1_PATIENT_0006'; % 21P17
%dicom_folder = 'T2STARMAP_ANATOMICAL_0015';               % 21P17
%dicom_folder = 'T2STARMAP_ANATOMICAL_0010';                % RYN

for nn = 1:length(Names)
    name = Names{nn};
    namedir=GetFullPath([resdir, name, '/']);
    if ~exist(namedir)
        mkdir(namedir);
    end
    for tp = 1:length(time_points)
        time_point = time_points{tp};
        tp_dir = GetFullPath([namedir, time_point, '/']);
        if ~exist(tp_dir)
            mkdir(tp_dir);
        end
        dicom_folder_glob = glob(cat(2, dicom_dir, '/KY_Patient_Data/', name, '/', time_point, '/', label, '/*/'));
        
        for nnn = 1:length(dicom_folder_glob)
            fname = dicom_folder_glob{nnn};
        [img_mag slice_data_mag] = dicom23D(fname, dicom_fields);

        matrix_size(1) = single(slice_data_mag(1).Width);
        matrix_size(2) = single(slice_data_mag(1).Height);
        matrix_size(3) = 1;

        eco = slice_data_mag(1).EchoTime;
        eco_array = [];
        eco_array(1) = eco;

        n = 2;

        while eco ~= slice_data_mag(n).EchoTime
            eco_array(n) = slice_data_mag(n).EchoTime;
            n = n + 1;
            if n > length(slice_data_mag)
                break;
            end
        end

        num_eco = length(eco_array);

        for i = 1:length(slice_data_mag)
            if mod(i, num_eco) == 0
                slice_data_mag(i).EchoNumber = num_eco;
            else
                slice_data_mag(i).EchoNumber = mod(i, num_eco);
            end
        end

        voxel_size(1,1) = single(slice_data_mag(1).PixelSpacing(1));
        voxel_size(2,1) = single(slice_data_mag(1).PixelSpacing(2));
        voxel_size(3,1) = single(slice_data_mag(1).SliceThickness);
        CF = slice_data_mag(1).ImagingFrequency *1e6;

        minSlice = 1e10;
        maxSlice = -1e10;
        rctr = 0; ictr=0;
        progress = '';
        TE = eco_array;
        NumEcho = num_eco;


        p_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);
        m_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);

        clear iMag;
        for i = 1:length(slice_data_mag)
            if slice_data_mag(i).SliceLocation<minSlice
                minSlice = slice_data_mag(i).SliceLocation;
                minLoc = slice_data_mag(i).ImagePositionPatient;
            end
            if slice_data_mag(i).SliceLocation>maxSlice
                maxSlice = slice_data_mag(i).SliceLocation;
                maxLoc = slice_data_mag(i).ImagePositionPatient;
            end

            if (slice_data_mag(i).ImageType(18)=='M')||(slice_data_mag(i).ImageType(18)=='m')
                ictr = ictr + 1;
                for ii=1:length(progress); fprintf('\b'); end
                progress=sprintf('Reading file %d', rctr+ictr);
                fprintf(progress);
                m_ImagePositionPatient(:,ictr) = slice_data_mag(i).ImagePositionPatient;
                i_EchoNumber(ictr) = slice_data_mag(i).EchoNumber;
                iMag(:,:,ictr)  = transpose(single(img_mag(:,:,ictr))); %magnitude
            end
        end

        fprintf('\n');

        matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3))+1;

        Affine2D = reshape(slice_data_mag(1).ImageOrientationPatient,[3 2]);
        Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
        B0_dir = Affine3D\[0 0 1]';

        m_slice = int32(round(sqrt(sum((m_ImagePositionPatient-minLoc).^2,1))/voxel_size(3)) +1);
        m_ind = sub2ind([matrix_size(3) NumEcho], m_slice(:), int32(i_EchoNumber(:)));
        iMag(:,:,m_ind) = iMag;


        iField = reshape(iMag, ...
            [matrix_size(1) matrix_size(2) matrix_size(3) NumEcho]);
        % iField = permute(iField,[2 1 3 4 5]); %This is because the first dimension is row in DICOM but COLUMN in MATLAB
        % iField(:,:,1:2:end,:) = -iField(:,:,1:2:end,:);
        if length(TE)==1
            delta_TE = TE;
        else
            delta_TE = TE(2) - TE(1);
        end


        % iField_trunc = iField(:,:,13:16,:);
        iField_trunc = iField(:,:,:,:);
        % figure(); imagesc(angle(iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))); caxis([-1 1]);
        % figure(); imagesc(angle(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2))); caxis([-1 1]);
        % figure(); imagesc(angle(iField_trunc(:,:,1,4)./iField_trunc(:,:,1,3))); caxis([-1 1]);
        % figure(); imagesc(angle(iField_trunc(:,:,1,5)./iField_trunc(:,:,1,4))); caxis([-1 1]);
        % figure(); imagesc(angle((iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))./(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2)))); caxis([-1 1]);

        f_central = CF;
        Mask_pre = ones((size(iField_trunc,1)-2), (size(iField_trunc, 2)-2), size(iField_trunc, 3));
        Mask = zeros((size(iField_trunc,1)), (size(iField_trunc, 2)), size(iField_trunc, 3));
        Mask(2:end-1, 2:end-1, :) = Mask_pre;

        % correct for eddy current
        mag_iField = abs(iField_trunc);
        mag_iField_norm = mag_iField ./ max(mag_iField(:));

        Fieldmap_eddy  = sqrt(sqrt((iField_trunc(:,:,:,2)./iField_trunc(:,:,:,1))./(iField_trunc(:,:,:,3)./iField_trunc(:,:,:,2))));
        iField_uneddy = zeros(size(iField_trunc));
        iField_uneddy(:,:,:,1:2:end) = iField_trunc(:,:,:,1:2:end).*Fieldmap_eddy;
        iField_uneddy(:,:,:,2:2:end) = iField_trunc(:,:,:,2:2:end)./Fieldmap_eddy;
        iField_trunc = mag_iField_norm.*exp(1i*angle(iField_uneddy));

        figure(); imagesc(angle(iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))); caxis([-1 1]);
        figure(); imagesc(angle(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2))); caxis([-1 1]);
        figure(); imagesc(angle(iField_trunc(:,:,1,4)./iField_trunc(:,:,1,3))); caxis([-1 1]);
        figure(); imagesc(angle(iField_trunc(:,:,1,5)./iField_trunc(:,:,1,4))); caxis([-1 1]);
        figure(); imagesc(angle((iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))./(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2)))); caxis([-1 1]);

        %
        % mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_roi.mat'));
        % mask_roi = zeros(size(iField_trunc,1), size(iField_trunc,2), size(iField_trunc,3));
        % if ~exist(mask_f)
        %     for slc = 1:size(iField_trunc,3)
        %         figure();
        %         imagesc(abs(iField_trunc(:,:,slc,3)));axis image;
        %         roi = drawpolygon;
        %         mask_roi(:,:,slc) = createMask(roi);
        %     end
        %     save(mask_f, 'mask_roi');
        % else
        %     load(mask_f);
        % end

        %
        SUBSAMPLE = 1;
        %dfat = [-244.3, -221.7, -175.4, -119.3, -32.1, 34] / (42.58*1.5) * 10^(-6) * f_central;
        dfat = -221.7 / (42.58*1.5) * 10^(-6) * f_central;
        LABEL = 1;

        [water,fat,iFreq,unwph_uf,unwph,N_std, R2s, fitting_error] = ...
            spurs_gc(iField_trunc(:,:,:,:),TE*1e-3,f_central,voxel_size, Mask>0, SUBSAMPLE, dfat, LABEL);
        % Save data
        slc = 1;
        figure;
        subplot(2,2,1);
        imagesc(abs(fat(:,:,slc)));colorbar;colormap jet; caxis([0,0.2]);
        subplot(2,2,2);
        imagesc(abs(water(:,:,slc)));colorbar;colormap jet;caxis([0,0.5]);
        subplot(2,2,3);
        %imagesc(abs(fib(:,:,slc)));colorbar;colormap jet;caxis([0 1])
        imagesc(abs(iFreq(:,:,slc)));colorbar;colormap jet;caxis([-pi pi])
        subplot(2,2,4);
        imagesc(abs(fitting_error(:,:,slc)));colorbar;colormap jet;

        ff = zeros(size(fat));
        fat_flag = fat > water;
        ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
        ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));

        figure(); imagesc(ff(:,:,1));caxis([0 0.5]);

        figure(); imagesc(fitting_error(:,:,1));caxis([0 0.2]);


        AllOtherMaps = struct;
        %AllOtherMaps.inputphase = inputphase;
        dTE_2echo=diff(TE)/1e-3;
        fmap_2echo = unwph/(2*pi*dTE_2echo(1)); % Hertz %using trimmed weighted mean
        B0map = fmap_2echo/f_central; %ppm Central frequency

        % AllOtherMaps.inputphase = inputphase;
        AllOtherMaps.B0map = B0map;
        AllOtherMaps.fmap_2echo = fmap_2echo;
        AllOtherMaps.Voxelsize = voxel_size;
        AllOtherMaps.water = water;
        AllOtherMaps.fat = fat;
        %AllOtherMaps.fib = fib;
        AllOtherMaps.unwph = unwph;
        AllOtherMaps.fitting_error = fitting_error;
        AllOtherMaps.R2s = R2s;
        AllOtherMaps.ff = ff;

        %%
        if FlagSaveB0
            strings = strsplit(fname, '/');
            filename = strings{end-1};
            save([tp_dir,  filename, '_SinglePeak.mat'], '-struct', 'AllOtherMaps','-v7.3');
        end
        
        end
    end
    close all;
end

