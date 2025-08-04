function [licenseID, flag, IDtag] = loadLicenseFile(reconOptions,flagDebug)

flag = 0;
IDtag = '';
fileName = [];

if nargin < 2
    flagDebug = false;
end

try
    licenseIDtag = reconOptions.licenseID;
catch
    licenseIDtag = 'None';
end

[licensePath,~]   = fileparts(mfilename('fullpath'));
licenseFolderInfo = dir(licensePath);
%licenseFolderInfo = dir([reconOptions.mainpath '/supporting/utils/license']);

for n = numel(licenseFolderInfo):-1:1   
    currentName = licenseFolderInfo(n).name;
    if contains(currentName,'Cedars') && contains(currentName,'.bin') && ~strncmp(currentName,'.',1)
        fileName = currentName;
        break;
    elseif ~contains(currentName,'.bin') || strncmp(currentName,'.',1)
        licenseFolderInfo(n) = [];
        if flagDebug
            disp(['Skip ' currentName]);
        end
    end
end
if isempty(fileName)
    for n = numel(licenseFolderInfo):-1:1  
        currentName = licenseFolderInfo(n).name;
        tempIDtag = [licenseIDtag '_Offline'];
        if contains(currentName,tempIDtag)
            fileName = currentName;
            break;
        end
        if contains(currentName,licenseIDtag)
            fileName = currentName;
        end
    end
end
if isempty(fileName) && ~isempty(licenseFolderInfo)
    for n = numel(licenseFolderInfo):-1:1  
        fileName = licenseFolderInfo(n).name;
        if contains(fileName,'Offline')
            break;
        end
    end
end
if ~isempty(fileName)
    try
        if (flagDebug)
            fprintf('License file detected. Loading %s...\n',fileName);
        end
        fileID = fopen([licensePath '/' fileName]);
        temp = fread(fileID,'double');
        fclose(fileID);
        IDtag = fileName(11:end-4);
        licenseID = reshape(temp,1,[]);
        if (flagDebug)
            fprintf('done.\n');
        end
    catch errormsg
        fprintf(2,'Error loading license file!\n');
        fprintf(2,'%s\n', errormsg.message);
        IDtag = 'TEMP';
    end
else
    fprintf(2,'License file not found!\n');
    IDtag = 'TEMP';
end

if strncmp(char(licenseID/sum(double('multitasking'))),'http',4)
    flag = 0;
else
    flag = 1;
end

