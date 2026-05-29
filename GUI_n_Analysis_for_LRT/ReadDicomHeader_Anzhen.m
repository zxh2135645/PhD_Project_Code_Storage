% ReadDicomHeader_Anzhen.m
% Read all DICOM files in /Users/jameszhang/Documents/Data/20251127/*/*/*/I*
% Print full filename when SeriesDescription == 'NOMI-DCE10'
%% Convert MAT file to DICOM (CINE and LGE)

% pattern = '/Users/jameszhang/Documents/Data/20251127_DICOM/*/*/*/I*';
% dirs = dir(pattern);
% 
% if isempty(dirs)
%     fprintf('No directories matched pattern: %s\n', pattern);
% end
% 
% matchCount = 0;
% % First glob the first-level folders under 20251127, then collect I* entries two levels down
% topPattern = '/Users/jameszhang/Documents/Data/20251127_DICOM/*/';
% topDirs = dir(topPattern);
% 
% 
% % Collect subject names as the parent directory of each matched I* entry
% subjectNames = {};
% for k = 1:numel(topDirs)
%     if isempty(topDirs(k).folder)
%         continue;
%     end
%     parts = strsplit(topDirs(k).folder, filesep);
%     parts = parts(~cellfun('isempty', parts)); % remove empty entries
%     if ~isempty(parts)
%         subjectNames{end+1} = parts{6}; %#ok<AGROW>
%     end
% end
% subjectNames = unique(subjectNames);
% fprintf('Extracted %d subject(s) (one directory above I* entries).\n', numel(subjectNames));
%% Convert MAT file to DICOM (CINE and LGE) New Version

pattern = '/Volumes/Extreme SSD/Anzhen/DICOM/nomidicom/*/*/*/I*';
dirs = dir(pattern);

if isempty(dirs)
    fprintf('No directories matched pattern: %s\n', pattern);
end

matchCount = 0;
% First glob the first-level folders under 20251127, then collect I* entries two levels down
topPattern = '/Volumes/Extreme SSD/Anzhen/DICOM/nomidicom/*/';
topDirs = dir(topPattern);


% Collect subject names as the parent directory of each matched I* entry
subjectNames = {};
for k = 1:numel(topDirs)
    if isempty(topDirs(k).folder)
        continue;
    end
    parts = strsplit(topDirs(k).folder, filesep);
    parts = parts(~cellfun('isempty', parts)); % remove empty entries
    if ~isempty(parts)
        subjectNames{end+1} = parts{6}; %#ok<AGROW>
    end
end
subjectNames = unique(subjectNames);
fprintf('Extracted %d subject(s) (one directory above I* entries).\n', numel(subjectNames));
%% Find matching subjects
% Read mapping table and match to subjectNames
% collect FID IDs from processed .mat files
procDir = '/Volumes/Extreme SSD/Anzhen/ProcessedData/20251130/20251130_ProcessedData/';
%procDir = '/Volumes/Extreme SSD/Anzhen/ProcessedData/20251201/20251201_ProcessedData/';

matFiles = dir(fullfile(procDir, '*.mat'));

processedFIDs = {};
for f = 1:numel(matFiles)
    name = matFiles(f).name;
    tok = regexp(name, 'FID(\d{5})', 'tokens', 'once');
    if ~isempty(tok)
        processedFIDs{end+1} = tok{1}; %#ok<AGROW>
    end
end
processedFIDs = unique(processedFIDs);
if isempty(processedFIDs)
    fprintf('No FIDxxxxx IDs found in %s\n', procDir);
else
    fprintf('Found %d unique FID IDs.\n', numel(processedFIDs));
end

% numeric form if needed
processedFIDs_num = str2double(processedFIDs);
%%
excelFile = '/Volumes/Extreme SSD/Anzhen/CENO_updated.xlsx';
T = readtable(excelFile,'ReadVariableNames',true);

col1 = T{:,1};   % first column (keys)
col2 = T{:,2};   % second column (values)

% canonical string and numeric representations of keys
keysStr = cellstr(string(col2));
keysNum = str2double(keysStr);

matchedValues = cell(numel(processedFIDs_num),1);

for i = 1:numel(processedFIDs_num)
    subj = processedFIDs_num(:,i);
    subjStr = string(subj);
    subjNum = str2double(subjStr);

    % try exact string match first
    idx = find(strcmp(keysStr, subjStr), 1);

    % if no string match and subject is numeric-like, try numeric match
    if isempty(idx) && ~isnan(subjNum)
        idx = find(keysNum == subjNum, 1);
    end

    if ~isempty(idx)
        matchedValues{i} = col1(idx); % preserve original type (numeric or cell)
    else
        matchedValues{i} = []; % no match
    end
end

% Export results to Excel (one row per subject)
outT = table(processedFIDs_num(:), matchedValues(:), 'VariableNames', {'Subject', 'MatchedValue'});
% outFile = '/Volumes/Extreme SSD/Anzhen/CENO_matched.xlsx';
% writetable(outT, outFile);
%fprintf('Wrote %d subject matches to %s\n', numel(subjectNames), outFile);
%% CINE
TriggerTime = 'TriggerTime';
NominalInterval = 'NominalInterval';
NumberOfTemporalPositions = 'NumberOfTemporalPositions';
InstanceNumber = 'InstanceNumber';
AcquisitionNumber = 'AcquisitionNumber';
SpacingBetweenSlices = 'SpacingBetweenSlices';
InversionTime = 'InversionTime';

fields_to_copy = {'IntervalsAcquired', 'IntervalsRejected', 'LowRRValue', 'HighRRValue', 'ImageComments'};

%for i = 1:numel(subjectNames)
for i = 1:numel(matchedValues)

    % find all entries in 'dirs' that belong to this subject
    % subj = string(matchedValues{i});
    % Normalize matchedValues{i} into a 10-digit, zero-padded string
    val = matchedValues{i};
    if iscell(val)
        val = val{1};
    end

    if isempty(val)
        subj = string(repmat('0',1,10));
    else
        if isnumeric(val)
            n = round(val);
        else
            sraw = char(string(val));
            md = regexp(sraw, '\d+', 'match');
            if ~isempty(md)
                n = str2double(md{1});
            else
                n = NaN;
            end
        end

        if ~isnan(n)
            subj = string(sprintf('%010d', n));
        else
            % fallback: extract digits from original text, pad/truncate to 10
            digits_only = regexp(char(string(val)), '\d', 'match');
            if isempty(digits_only)
                s = '0';
            else
                s = [digits_only{:}];
            end
            if numel(s) > 10
                s = s(end-9:end);
            else
                s = [repmat('0',1,10-numel(s)), s];
            end
            subj = string(s);
        end
    end

    %%
    folders = {dirs.folder};
    % match folder components that contain the subject name as a separate path part
    % pat = [filesep regexptranslate('escape', subj) filesep];
    idx = find(~cellfun('isempty', regexp(folders, subj)));
    if isempty(idx)
        fprintf('No I* entries found for subject %s\n', subj);
        continue;
    end
    % save all matched dirs for this subject
    subjectDirs{i} = dirs(idx);

    cine_check = 0;
    lge_check = 0;
    dicom_headers = {};
    lrt_match = 0;
    for ii = 1:numel(dirs(idx))
        tmp = dirs(idx(ii));
        filePath = fullfile(tmp.folder, tmp.name);

        try
            info = dicominfo(filePath);
        catch
            % not a DICOM file or cannot read header -> skip
            continue;
        end

        if cine_check == 0
            if isfield(info, 'SeriesDescription') && contains(strtrim(info.SeriesDescription), 'cine_')
                fprintf('%s\n', filePath);
                %matchCount = matchCount + 1;
                dicom_headers_cine = info;
                cine_check = 1;
            end
        end

        if lge_check == 0
            if isfield(info, 'SeriesDescription') && contains(strtrim(info.SeriesDescription), 'sax_MAG')
                fprintf('%s\n', filePath);
                %matchCount = matchCount + 1;
                dicom_headers_lge = info;
                lge_check = 1;
            end
        end


        if isfield(info, 'SeriesDescription') && strcmp(strtrim(info.SeriesDescription), 'NOMI-DCE10')
            fprintf('%s\n', filePath);
            matchCount = matchCount + 1;
            lrt_match = lrt_match + 1;
            dicom_headers{lrt_match} = info;

        end

        % sort dicom_headers by SliceLocation (fall back to ImagePositionPatient(3))
        nHdr = numel(dicom_headers);
        locs = nan(1,nHdr);
        for jj = 1:nHdr
            hdr = dicom_headers{jj};
            if isfield(hdr, 'SliceLocation') && ~isempty(hdr.SliceLocation)
                locs(jj) = hdr.SliceLocation;
            elseif isfield(hdr, 'ImagePositionPatient') && numel(hdr.ImagePositionPatient) >= 3
                locs(jj) = hdr.ImagePositionPatient(3);
            end
        end
        [~, order] = sort(locs, 'ascend'); % use 'descend' if reverse order desired
        dicom_headers = dicom_headers(order);
    end

    %% Load LRT
    splitted_str = strsplit(procDir,'/');
    fid_path = cat(2, GetFullPath(cat(2, topDirs(1).folder, '/../../..')), '/ProcessedData/', splitted_str{6}, '/', splitted_str{7}, '/');
    fid = num2str(processedFIDs_num(:,i));
    fid_file = cat(2, 'FID', fid, '*.mat');
    f_glob = glob(strcat(fid_path, fid_file));
    f_name = strsplit(f_glob{1},'/');
    f_name = f_name{end};
    load(f_glob{1}, 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int', 'Phi_rt_small_init');

    if ~exist('RR_int')
        SGblock = params.NavInterval;
        fs=1/(params.lEchoSpacing*SGblock);
        df=fs/size(Phi_rt_small_init,2);
        total_time = params.lEchoSpacing*SGblock * size(Phi_rt_small_init,2);
        % RR_int = 1000/((find(Hidx(56:end) == max(Hidx(56:end/2)),1)-1)*df)
        [idx_min, b_min] = find(Hidx(1:end) == min(Hidx(1:end/2)), 1);
        RR_int = 60/((find(Hidx(b_min:end) == max(Hidx(b_min:end/2)),1)-1)*df) 
    end
    %%
    % Copy selected fields from dicom_headers_cine{1} to dicom_headers{1}
    dicom_headers_to_copy = cell(1, length(dicom_headers)*size(Phi,3));

    for slc = 1:length(dicom_headers)
        for frame = 1:size(Phi,3)
            idxx = (slc-1)*size(Phi,3) + frame;
            dicom_headers_to_copy{idxx} = dicom_headers{slc};
            dicom_headers_to_copy{idxx}.(TriggerTime) = (frame-1) * RR_int / size(Phi,3);
            dicom_headers_to_copy{idxx}.(NominalInterval) = RR_int;
            dicom_headers_to_copy{idxx}.(NumberOfTemporalPositions) = size(Phi,3);
            dicom_headers_to_copy{idxx}.(InstanceNumber) = idxx;
            dicom_headers_to_copy{idxx}.(AcquisitionNumber) = idxx;
            dicom_headers_to_copy{idxx}.(SpacingBetweenSlices) = 6;

            for ii = 1:length(fields_to_copy)
                field = fields_to_copy{ii};
                if isfield(dicom_headers_cine, field)
                    dicom_headers_to_copy{idxx}.(field) = dicom_headers_cine.(field);
                end
            end
        end
    end

    %% Display image (Write DICOM)
    slc_array = fftshift(1:Nz);
    % slc_array = 1:Nz;
    % num_seg_array = [61, 81, 101, 121, 201];
    num_seg_array = [101];
    % num_seg_array = [141, 151, 161, 171, 181, 191];
    temtemp_4D = zeros(Ny, Nx, size(Phi,3), length(slc_array));

    % cardiac phase and resp phase needs to be encoded
    resp_phase = 1;
    card_phase = 24;

    for iii = 1:length(slc_array)
        slc = slc_array(iii);
        dispim = @(x)fftshift(x(:,:,slc,:),1);
        for num = 1:length(num_seg_array)
            temp = Gr\reshape(Phi(:,num_seg_array(num),:,resp_phase,2), L, []);
            temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);

            temtemp_4D(:,:,:,iii) = temp;
        end
    end

    % Crop first dimension of temtemp_4D to dicom_headers_to_copy{1}.Rows, centered
    target_rows = dicom_headers_to_copy{1}.Rows;
    current_rows = size(temtemp_4D, 1);
    row_start = floor((current_rows - target_rows)/2) + 1;
    row_end = row_start + target_rows - 1;
    temtemp_4D = temtemp_4D(row_start:row_end, :, :, :);

    %% Write CINE into DICOM
    X = dicom_headers{1}.Filename;
    NEco_old = params.NEco_old;
    %len = length(slice_data{1})/NEco_old;

    save_dir = GetFullPath(cat(2, fid_path, 'DICOM_CINE_FID', fid, '/'));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    for te = 1:length(num_seg_array)
        num_seg = num_seg_array(te);
        for slc = 1:length(dicom_headers)
            %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
            %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
            

            f_dir = GetFullPath(cat(2, save_dir, 'Slice_', num2str(slc,'%02.0f'), '/'));

            if ~exist(f_dir, 'dir')
                mkdir(f_dir);
            end

            for j = 1:size(Phi,3)
                fname = GetFullPath(cat(2, f_dir, f_name(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(slc,'%02.0f'), '_CINE', num2str(j,'%02.0f'), '.dcm'));
                %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));

                metadata = dicom_headers_to_copy{(slc-1)*size(Phi,3) + j};

                %slc_array_flip = flip(slc_array);
                %slc = slc_array_flip(i);
                if length(dicom_headers) < Nz
                    slc_temp = slc + (Nz - length(dicom_headers));
                else
                    slc_temp = slc;
                end

                temp = imrotate(temtemp_4D(:,:,j,slc_temp), 90);
                cw = 0.8*max(vec(abs(temtemp_4D(:,:,j,slc_temp))));

                % temp = temp_4D(:,:,slc,te);
                temp = abs(flip(temp,2)./cw);
                temp = uint16(temp*4095);
                metadata.WindowCenter = 2048;
                metadata.WindowWidth = 4095;

                dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
                metadata.SmallestImagePixelValue = min(temp(:));
                metadata.LargestImagePixelValue = max(temp(:));
                dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
            end
        end
    end


    %% Convert to LGE
    num_seg_array_lge = [11,16,21,26,31,41,51,61,71,81,101];
    dicom_headers_to_copy_lge = cell(1, length(dicom_headers)*length(num_seg_array_lge));

    for slc = 1:length(dicom_headers)
        for frame = 1:length(num_seg_array_lge)
            idxx = (slc-1)*length(num_seg_array_lge) + frame;
            inv_idx = num_seg_array_lge(frame);
            dicom_headers_to_copy_lge{idxx} = dicom_headers{slc};
            dicom_headers_to_copy_lge{idxx}.(InversionTime) = inv_idx * params.lEchoSpacing*SGblock*1000; % in ms

            for ii = 1:length(fields_to_copy)
                field = fields_to_copy{ii};
                if isfield(dicom_headers_lge, field)
                    dicom_headers_to_copy_lge{idxx}.(field) = dicom_headers_lge.(field);
                end
            end
        end
    end

    %% Display image (Write DICOM)
    % slc_array = fftshift(1:Nz);
    slc_array = 1:Nz;
    % num_seg_array = [61, 81, 101, 121, 201];
    % num_seg_array = [101];
    % num_seg_array = [141, 151, 161, 171, 181, 191];
    enhancement_array = [1:size(Phi,5)];
    temtemp_5D = zeros(Ny, Nx, length(num_seg_array_lge), length(slc_array), size(Phi,5));

    % cardiac phase and resp phase needs to be encoded
    resp_phase = 1;
    card_phase = 24;

    for iii = 1:length(slc_array)
        slc = slc_array(iii);
        dispim = @(x)fftshift(x(:,:,slc,:),1);
        for num = 1:length(num_seg_array_lge)
            temp = Gr\reshape(Phi(:,num_seg_array_lge(num),card_phase,resp_phase,:), L, []);
            temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
            temtemp_5D(:,:,num,iii,:) = temp;
        end
    end

    % Crop first dimension of temtemp_4D to dicom_headers_to_copy{1}.Rows, centered
    target_rows = dicom_headers_to_copy_lge{1}.Rows;
    current_rows = size(temtemp_5D, 1);
    row_start = floor((current_rows - target_rows)/2) + 1;
    row_end = row_start + target_rows - 1;
    temtemp_5D = temtemp_5D(row_start:row_end, :, :, :, :);

    %% Write CINE into DICOM
    X = dicom_headers{1}.Filename;
    NEco_old = params.NEco_old;
    %len = length(slice_data{1})/NEco_old;

    save_dir = GetFullPath(cat(2, fid_path, 'DICOM_LGE_FID', fid, '/'));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    for enh = 1:length(enhancement_array)
        enh_dir = GetFullPath(cat(2, save_dir, 'Enhancement_', num2str(enh,'%02.0f'), '/'));

        if ~exist(enh_dir, 'dir')
            mkdir(enh_dir);
        end

        for slc = 1:length(dicom_headers)
            %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
            %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);

            f_dir = GetFullPath(cat(2, enh_dir, 'Slice_', num2str(slc,'%02.0f'), '/'));

            if ~exist(f_dir, 'dir')
                mkdir(f_dir);
            end

            for j = 1:length(num_seg_array_lge)
                num_seg = num_seg_array_lge(j);
                % for j = 1:size(Phi,3)
                fname = GetFullPath(cat(2, f_dir, f_name(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(slc,'%02.0f'), '_CINE', num2str(j,'%02.0f'), '.dcm'));
                %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));

                metadata = dicom_headers_to_copy_lge{(slc-1)*length(num_seg_array_lge) + j};

                %slc_array_flip = flip(slc_array);
                %slc = slc_array_flip(i);
                if length(dicom_headers) < Nz
                    slc_temp = slc + (Nz - length(dicom_headers));
                else
                    slc_temp = slc;
                end

                temp = imrotate(temtemp_5D(:,:,j,slc_temp,enh), 90);
                cw = 0.8*max(vec(abs(temtemp_5D(:,:,j,slc_temp,enh))));

                % temp = temp_4D(:,:,slc,te);
                temp = abs(flip(temp,2)./cw);
                temp = uint16(temp*4095);
                metadata.WindowCenter = 2048;
                metadata.WindowWidth = 4095;

                dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
                metadata.SmallestImagePixelValue = min(temp(:));
                metadata.LargestImagePixelValue = max(temp(:));
                dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
            end
        end
    end
end
fprintf('Found %d file(s) with SeriesDescription ''NOMI-DCE10''.\n', matchCount);
% %% Convert to LGE
% %%
% matchCount = 0;
% InversionTime = 'InversionTime';
%
% fields_to_copy = {'IntervalsAcquired', 'IntervalsRejected', 'LowRRValue', 'HighRRValue', 'ImageComments'};
%
% %for i = 1:numel(subjectNames)
% for i = 1:numel(subjectNames)
%
%     % find all entries in 'dirs' that belong to this subject
%     subj = subjectNames{i};
%     folders = {dirs.folder};
%     % match folder components that contain the subject name as a separate path part
%     pat = [filesep regexptranslate('escape', subj) filesep];
%     idx = find(~cellfun('isempty', regexp(folders, pat)));
%     if isempty(idx)
%         fprintf('No I* entries found for subject %s\n', subj);
%         %continue;
%     end
%     % save all matched dirs for this subject
%     subjectDirs{i} = dirs(idx);
%
%     lge_check = 0;
%     dicom_headers = {};
%     lrt_match = 0;
%     for ii = 1:numel(dirs(idx))
%         tmp = dirs(idx(ii));
%         filePath = fullfile(tmp.folder, tmp.name);
%
%         try
%             info = dicominfo(filePath);
%         catch
%             % not a DICOM file or cannot read header -> skip
%             continue;
%         end
%
%         if lge_check == 0
%             if isfield(info, 'SeriesDescription') && contains(strtrim(info.SeriesDescription), 'sax_MAG')
%                 fprintf('%s\n', filePath);
%                 %matchCount = matchCount + 1;
%                 dicom_headers_lge = info;
%                 lge_check = 1;
%             end
%         end
%
%
%         if isfield(info, 'SeriesDescription') && strcmp(strtrim(info.SeriesDescription), 'NOMI-DCE10')
%             fprintf('%s\n', filePath);
%             matchCount = matchCount + 1;
%             lrt_match = lrt_match + 1;
%             dicom_headers{lrt_match} = info;
%         end
%
%     end
%
%     %% Load LRT
%     fid_path = cat(2, GetFullPath(cat(2, topDirs(1).folder, '/..')), '_ProcessedData/');
%     fid_file = cat(2, 'FID', subj, '*.mat');
%     f_glob = glob(strcat(fid_path, fid_file));
%     f_name = strsplit(f_glob{1},'/');
%     f_name = f_name{end};
%     load(f_glob{1}, 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int', 'Phi_rt_small_init');
%
%     if ~exist('RR_int')
%         SGblock = params.NavInterval;
%         fs=1/(params.lEchoSpacing*SGblock);
%         df=fs/size(Phi_rt_small_init,2);
%         total_time = params.lEchoSpacing*SGblock * size(Phi_rt_small_init,2);
%         % RR_int = 1000/((find(Hidx(56:end) == max(Hidx(56:end/2)),1)-1)*df)
%         [idx_min, b_min] = find(Hidx(1:end) == min(Hidx(1:end/2)), 1);
%         RR_int = 60/((find(Hidx(b_min:end) == max(Hidx(b_min:end/2)),1)-1)*df) 
%     end
%     %%
%     % Copy selected fields from dicom_headers_lge{1} to dicom_headers{1}
%     % dicom_headers_to_copy = cell(1, length(dicom_headers)*size(Phi,3));
% 
%     num_seg_array_lge = [11,16,21,26,31,41,51,61,71,81,101];
%     dicom_headers_to_copy_lge = cell(1, length(dicom_headers)*length(num_seg_array_lge));
% 
%     for slc = 1:length(dicom_headers)
%         for frame = 1:length(num_seg_array_lge)
%             idxx = (slc-1)*length(num_seg_array_lge) + frame;
%             inv_idx = num_seg_array_lge(frame);
%             dicom_headers_to_copy_lge{idxx} = dicom_headers{slc};
%             dicom_headers_to_copy_lge{idxx}.(InversionTime) = inv_idx * params.lEchoSpacing*SGblock*1000; % in ms
% 
%             for ii = 1:length(fields_to_copy)
%                 field = fields_to_copy{ii};
%                 if isfield(dicom_headers_lge, field)
%                     dicom_headers_to_copy_lge{idxx}.(field) = dicom_headers_lge.(field);
%                 end
%             end
%         end
%     end
% 
%     %% Display image (Write DICOM)
%     % slc_array = fftshift(1:Nz);
%     slc_array = 1:Nz;
%     % num_seg_array = [61, 81, 101, 121, 201];
%     % num_seg_array = [101];
%     % num_seg_array = [141, 151, 161, 171, 181, 191];
%     enhancement_array = [1:size(Phi,5)];
%     temtemp_5D = zeros(Ny, Nx, length(num_seg_array_lge), length(slc_array), size(Phi,5));
% 
%     % cardiac phase and resp phase needs to be encoded
%     resp_phase = 1;
%     card_phase = 24;
% 
%     for iii = 1:length(slc_array)
%         slc = slc_array(iii);
%         dispim = @(x)fftshift(x(:,:,slc,:),1);
%         for num = 1:length(num_seg_array_lge)
%             temp = Gr\reshape(Phi(:,num_seg_array_lge(num),card_phase,resp_phase,:), L, []);
%             temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
%             temtemp_5D(:,:,num,iii,:) = temp;
%         end
%     end
% 
%     % Crop first dimension of temtemp_4D to dicom_headers_to_copy{1}.Rows, centered
%     target_rows = dicom_headers_to_copy_lge{1}.Rows;
%     current_rows = size(temtemp_5D, 1);
%     row_start = floor((current_rows - target_rows)/2) + 1;
%     row_end = row_start + target_rows - 1;
%     temtemp_5D = temtemp_5D(row_start:row_end, :, :, :, :);
% 
%     %% Write CINE into DICOM
%     X = dicom_headers{1}.Filename;
%     NEco_old = params.NEco_old;
%     %len = length(slice_data{1})/NEco_old;
% 
%     save_dir = GetFullPath(cat(2, fid_path, 'DICOM_LGE_FID', subj, '/'));
%     if ~exist(save_dir, 'dir')
%         mkdir(save_dir);
%     end
% 
% 
%     for slc = 1:Nz
%         %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
%         %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
% 
%         f_dir = GetFullPath(cat(2, save_dir, 'Slice_', num2str(slc,'%02.0f'), '/'));
% 
%         if ~exist(f_dir, 'dir')
%             mkdir(f_dir);
%         end
%         for enh = 1:length(enhancement_array)
%             enh_dir = GetFullPath(cat(2, f_dir, 'Enhancement_', num2str(enh,'%02.0f'), '/'));
% 
%             if ~exist(enh_dir, 'dir')
%                 mkdir(enh_dir);
%             end
%             for j = 1:length(num_seg_array_lge)
%                 num_seg = num_seg_array_lge(j);
%                 % for j = 1:size(Phi,3)
%                 fname = GetFullPath(cat(2, enh_dir, f_name(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(slc,'%02.0f'), '_CINE', num2str(j,'%02.0f'), '.dcm'));
%                 %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));
% 
%                 metadata = dicom_headers_to_copy_lge{(slc-1)*length(num_seg_array_lge) + j};
% 
%                 %slc_array_flip = flip(slc_array);
%                 %slc = slc_array_flip(i);
%                 temp = imrotate(temtemp_5D(:,:,j,slc,enh), 90);
%                 cw = 0.8*max(vec(abs(temtemp_5D(:,:,j,slc,enh))));
% 
%                 % temp = temp_4D(:,:,slc,te);
%                 temp = abs(flip(temp,2)./cw);
%                 temp = uint16(temp*4095);
%                 metadata.WindowCenter = 2048;
%                 metadata.WindowWidth = 4095;
% 
%                 dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
%                 metadata.SmallestImagePixelValue = min(temp(:));
%                 metadata.LargestImagePixelValue = max(temp(:));
%                 dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
%             end
%         end
%     end
% end
% 
%fprintf('Found %d file(s) with SeriesDescription ''NOMI-DCE10''.\n', matchCount);