clear all 
close all

% This is initially used for generating Phase map from AllPhaseMap, for the
% uses of Greg, can be potentially applied to other senarios

% Initialization & and image recon
addpath(genpath(pwd));
rmpath(genpath('./Example/'));
addpath(genpath('../reconstructionPipeline/'));
addpath('../function/');


SOS = @(x) sqrt(sum(abs(x).^2,3));

base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

time_label = {'8WK'};
names_to_rule_out = {'FF_Data_FromIDEAL', 'FreshHeart_Exvivo', 'Dharmakumar_Cedars_18D16_10month_Invivo'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
%% Load image
clear phase_unwrapped phase_unwrapped_mn_3d phase_unwrapped_wt_3d B0mapstd;
subject_name_cell = {'18D15_12month_Invivo', 'Gobi_12month_Invivo', 'Hope_12month_Invivo', 'Merry_12month_Invivo', 'Mojave_12month_Invivo',...
    'Ryn_12month_Invivo', 'Sahara_12month_Invivo', 'Sunny_12month_Invivo'};

for n = 1:length(subject_name_cell)
    name = subject_name_cell{n};
    raw_glob = glob(cat(2, name_glob{n}, 'RawData/meas_*'));
    save_dir = GetFullPath(cat(2, name_glob{n}, 'Result/'));
    clear AllPhasemap_ForGreg;
    for i = 1:length(raw_glob)

        filename = raw_glob{i};
        twix_obj_in = mapVBVD(filename, 'removeOS'); % return all image-data:
        if (length(twix_obj_in)>1)% R.Y. avoid adj coil sensitivity
            for  k=1:length(twix_obj_in)
                if (~strcmp(twix_obj_in{k}.hdr.MeasYaps.tSequenceFileName,'%AdjustSeq%/AdjCoilSensSeq') )
                    twix_obj=twix_obj_in{k};
                end
            end
        else
            twix_obj=twix_obj_in;
        end

        rawdata = squeeze(twix_obj.image(''));

        if twix_obj.image.NSeg == 1 && twix_obj.image.NAve == 1
            [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);
        elseif twix_obj.image.NSeg == 1 && twix_obj.image.NAve > 1
            [NumRO, NumCh, NumPE, NumSlices, NumAve, NumEchos] = size(rawdata);
        elseif twix_obj.image.NAve == 1 && (twix_obj.image.NPar && twix_obj.image.NPar == 1)
            [NumRO, NumCh, NumPE, NumEchos, NSeg] = size(rawdata);
            NumSlices = 1;
        end
        clear twix_obj_in rawdata;
        nameSession         = ['-' num2str(twix_obj.hdr.Config.MeasUID)];
        temp = strsplit(filename, '/');
        strings = strsplit(temp{end-2},'_');
        PatientName_new = cat(2, strings{end-2}, '_', strings{end-1}, '_', strings{end});
        new_path = GetFullPath(cat(2, name_glob{n}, '../FF_Data_FromIDEAL/', PatientName_new, '/'));

        if exist(cat(2, new_path, 'xspace', nameSession), 'dir')
            xspace_glob = glob(cat(2, new_path, 'xspace', nameSession, '/*'));
        end

        xspace_cell = cell(length(xspace_glob), 1);
        for eco = 1:length(xspace_glob)
            xspace_cell{eco} = load(xspace_glob{eco});
        end

        RawdataFT = zeros([size(xspace_cell{1}.data), 1, NumEchos]);

        for eco = 1:NumEchos
            RawdataFT(:,:,:,:,eco) = xspace_cell{eco}.data;
        end

        RawdataFT = permute(RawdataFT, [2,3,4,1,5]);% NumRO, NumPE, NumSlices, NumCh, NumEcho
        RawdataFT = flip(flip((RawdataFT),1),2);

        AllPhasemap(i).Name = temp{end}(1:end-4);
        AllPhasemap(i).compleximg = RawdataFT;% In object domain
        for necho=1:NumEchos
            AllPhasemap(i).TE(necho)=twix_obj.hdr.MeasYaps.alTE{necho} / 1000;%usec
        end
        AllPhasemap(i).RoFOV = twix_obj.hdr.Config.RoFOV;%usec
        AllPhasemap(i).PeFOV = twix_obj.hdr.Config.PeFOV;%usec
        AllPhasemap(i).Voxelsize=[AllPhasemap(i).RoFOV/twix_obj.hdr.Meas.NImageCols, AllPhasemap(i).PeFOV/twix_obj.hdr.Meas.NImageLins ];

        iField = permute(RawdataFT,[1,2,3,5,4])*100;

        if size(iField,5)>1
            % combine multiple coils together, assuming the coil is the fifth dimension
            iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
            mag_iField = abs(iField);
            mag_iField_norm = mag_iField ./ max(mag_iField(:));

            Fieldmap_eddy  = sqrt(sqrt((iField(:,:,:,2)./iField(:,:,:,1))./(iField(:,:,:,3)./iField(:,:,:,2))));
            iField_uneddy = zeros(size(iField));
            iField_uneddy(:,:,:,1:2:end) = iField(:,:,:,1:2:end).*Fieldmap_eddy;
            iField_uneddy(:,:,:,2:2:end) = iField(:,:,:,2:2:end)./Fieldmap_eddy;
            iField = mag_iField_norm.*exp(1i*angle(iField_uneddy));
        end

        AllPhasemap(i).MAG = squeeze(abs(iField));
        AllPhasemap(i).PHASE = squeeze(angle(iField));

    end
    save([save_dir,'AllPhasemap_ForGreg.mat'], 'AllPhasemap','-v7.3');
end


