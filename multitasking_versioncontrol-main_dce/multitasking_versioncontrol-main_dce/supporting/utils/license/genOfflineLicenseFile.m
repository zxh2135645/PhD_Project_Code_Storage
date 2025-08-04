function fileName = genOfflineLicenseFile(licenseID,info,yyyy,mm,dd)
% Usage:
%         Generate license file for offline check using machine info
%         and expiration date yyyy-mm-dd, 
%         with file name ['licenseID_' licenseID '_Offline.bin'].
%
%         If only less than 5 input arguments are given, 
%         the expiration date will be yyyy days from the current date.
%
%         If a license file with the same licenseID already exists and 
%         hasn't expired, new machine info will be appended to the file
%

if isempty(info)
    if ispc
        [~,info] = system('getmac');
    elseif ismac
        [~,info] = system('ifconfig en0 | grep ether');
    elseif isunix
        [~,info] = system('ip addr | grep ether');
    else
        error('OS not recognized.');
    end
end

if nargin < 5
    datestr = datenum(date) + yyyy;
else
    datestr = datenum(yyyy,mm,dd);
end

[licensePath,~]   = fileparts(mfilename('fullpath'));
licenseFolderInfo = dir(licensePath);
try
    for n = 1:numel(licenseFolderInfo)   
        if contains(licenseFolderInfo(n).name,licenseID) && contains(licenseFolderInfo(n).name,'Offline.bin')
            fileID = fopen(licenseFolderInfo(n).name);
            temp = fread(fileID,'double');
            fclose(fileID);
            break;
        end
    end
    oldhash = reshape(temp,1,[]) / sum(double('multitasking')); 
    if oldhash(end) - datenum(date) > 0
        datestr = oldhash(end);
        oldhash(end) = [];
    else
        oldhash = [];
    end
catch
    oldhash = [];
end

if datestr < datenum(date)
    error('Expiration date must be later than the current date!');
else
    ID   = double(info) * datestr;
    md   = java.security.MessageDigest.getInstance('MD5');
    hash = dec2hex(uint8(double(md.digest(ID))+128));
    hash = hash(:).';
            
    if ~contains(char(oldhash),hash)
        hash = [char(oldhash) hash];
    else
        hash = char(oldhash);
    end
    temp = [double(hash) datestr] * sum(double('multitasking'));

    fileName =  ['licenseID_' licenseID '_Offline.bin'];
    fileID = fopen(erase(join(fileName)," "),'w');
    fwrite(fileID, temp, 'double');
    fclose(fileID);
end


% info used to generate Cedars_BIRI license
% datestr = 748783
% recon2 hash: '295892E4E070B5C2EA976BB656DA66C3'
% recon3 hash: '5B7FA010DD2B67EAFE818D9097B95F2C'
