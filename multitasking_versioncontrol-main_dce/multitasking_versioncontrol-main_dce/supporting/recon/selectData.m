function [params, status, msg] = selectData(flagLogOn)

if nargin < 1
    flagLogOn = true;
end

status = 0; msg = ' File selected';
params = [];

com.mathworks.mwswing.MJFileChooserPerPlatform.setUseSwingDialog(1);
[fid_file, fid_path] = uigetfile('*.dat;*.mat');
drawnow; pause(0.1);
fidString  = [fid_path fid_file];
datetimeStr = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
try
    if strcmp(fid_file((end-3):end),'.dat')
        params.filePath   = [fid_path fid_file(1:(end-4)) '_' datetimeStr];
        params.fileString = fullfile(params.filePath,'multitasking');
        mkdir(params.filePath);    
        if flagLogOn
            diary([params.fileString '.log']);
        end
    elseif strcmp(fid_file((end-3):end),'.mat')
        load(fidString, 'Params');
        params = Params;
        params.filePath   = [fid_path fid_file(1:(end-4)) '_' datetimeStr];
        params.fileString = fullfile(params.filePath,'multitasking');
        mkdir(params.filePath);    
        if flagLogOn
            diary([params.fileString '.log']);
        end
    else
        status = 1; msg = ' Wrong file path';
    end
    params.isSortPart   = false;
    params.isDownsample = false;
    params.fidString    = fidString;
catch
    status = 1; msg = ' Wrong file path';
end

    

