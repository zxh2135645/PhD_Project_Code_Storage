
directory = '/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/';      % Full path of the directory to be searched in
stringToBeFound = input("Specify the string you would like to be searched for: ", 's');

dir_glob = glob(cat(2, directory, '*'));

% TODO make it recursive
for dd = 1:(length(dir_glob)+1)
    if dd ~= (length(dir_glob)+1)
        if isdir(dir_glob{dd})
            % Returns all the files and folders in the directory
            filesAndFolders = dir(dir_glob{dd}); 
        else
            filesAndFolders = [];
        end
    else
        filesAndFolders = dir(directory);
    end
    if ~isempty(filesAndFolders)
        % foldersinDir = filesAndFolders(([filesAndFolders.isdir]));
        filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
        numOfFiles = length(filesInDir);
        i=1;
        strings = strsplit(filesAndFolders(1).folder, '/');
        dirString = strcat('In Folder  ', strings{end});
        disp(dirString);
        while(i<=numOfFiles)
            filename = filesInDir(i).name;                              % Store the name of the file
            fid = fopen(cat(2, filesAndFolders(1).folder, '/', filename));
            while(~feof(fid))                                           % Execute till EOF has been reached
                contentOfFile = fgetl(fid);                             % Read the file line-by-line and store the content
                found = strfind(contentOfFile,stringToBeFound);         % Search for the stringToBeFound in contentOfFile
                if ~isempty(found)
                    foundString = strcat('Found in file ------ ', filename);
                    disp(foundString);
                    break;
                end
            end
            fclose(fid);                                                % Close the file
            i = i+1;
        end
    end
end