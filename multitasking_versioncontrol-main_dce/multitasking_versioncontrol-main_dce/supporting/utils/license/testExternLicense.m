function isValid = testExternLicense(reconOptions)

if nargin < 1
    reconOptions = [];
end

flagDebug = true;

isValid = false;

% get machine ID
[info,strAdd] = genID;
info = char(info/sum(double('multitasking')));

md = java.security.MessageDigest.getInstance('MD5');    

try
    % load license file
    [licenseID,flagOffline] = loadLicenseFile(reconOptions,1);
    fprintf('Parsing license ID... \n');
    [hashes,dateNr] = getHash(licenseID,flagDebug);
    if ~isempty(hashes)
        fprintf('checking hash ');
        for n = 1:length(strAdd)
            ID   = double(strAdd{n})*sum(dateNr);
            hash = dec2hex(uint8(double(md.digest(ID))+128));
            hash = hash(:).';
            fprintf('#%d.. ',n);
            isValid = contains(hashes,hash);
            if isValid; break; end
        end
        fprintf('\n');
    end
    if ~isValid && flagOffline == 1
        fprintf(2,'Offline license check failed. System info does not match.\n');
    elseif ~isValid
        fprintf(2,'Online license check failed. System info does not match.\n');
    else    % isValid
        disp('License check successful.');
        if flagOffline == 1
            datestr = char(datetime(dateNr,'ConvertFrom','datenum','Format','dd-MMM-yyyy'));
            disp(['License expiration date: ' datestr]);
        end
    end
catch
    isValid = false;
    fprintf(2,'License check failed.\n');
end

disp(['Current date: ' char(datetime("today"))]);
disp('Current system info:');
disp(info);

