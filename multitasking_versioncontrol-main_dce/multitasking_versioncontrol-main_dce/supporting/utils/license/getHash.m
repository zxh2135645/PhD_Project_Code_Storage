function [hashes,dateNr] = getHash(licenseID,flagDebug)

if nargin < 2
    flagDebug = false;
end

licenseID = licenseID/sum(double('multitasking'));

if strncmp(char(licenseID),'http',4)    % online license check
    % get hashes from server
    % hashes  = webread(['https://agchristodoulou.github.io/MTcheck/' htmlID '.html']);
    URL = char(licenseID);
    dateNr = double(date);
    try
        hashes = webread(URL);
        if flagDebug
            fprintf('Online license check: hashes retrieved from server: \n');
            fprintf('%s', hashes);
        end
    catch errormsg
        hashes = [];
        fprintf(2,'Online license check error: %s\n', errormsg.message);
    end 
else                                    % offline license check
    dateNr = licenseID(end);
    licenseID(end) = [];
    licenseDays = dateNr - datenum(date);
    if licenseDays > 0 
        hashes = char(licenseID);
        if flagDebug
            fprintf('Offline license check: hashes retrieved from local: \n');
            fprintf('%s\n', hashes);
        end
    else
        hashes = [];
        fprintf(2,'License expired.\n');
    end
end