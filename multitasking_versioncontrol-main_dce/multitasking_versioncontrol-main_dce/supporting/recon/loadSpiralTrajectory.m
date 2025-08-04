function kTraj = loadSpiralTrajectory(params,showTraj)

if nargin < 2
    showTraj = false;
end

fileName = [];

if isstruct(params)
    [filePath,~] = fileparts(params.fidString);
    filePathInfo = dir(filePath);  
    strMID = ['MID' num2str(params.MID)];
    for n = 1:numel(filePathInfo)
        currentName = filePathInfo(n).name;
        if contains(currentName,'spiralTraj') && contains(currentName,strMID) && contains(currentName,'.traj')
            fileName = currentName;
            break;
        end
    end

    currentPath = pwd;
    if isempty(fileName)
        cd(filePath);
        [fileName, filePath] = uigetfile('*.traj;*.mat', 'Select spiral trajectory file');
    end
    cd(currentPath);
elseif ischar(params)
    [filePath,fileName,ext] = fileparts(params);
    fileName = [fileName ext];
    if isempty(filePath)
        filePath = '.';
    end
end

try
    underscore = strfind(fileName,'_');
    dext       = strfind(fileName,'.');
    spiralLength = str2double(fileName(underscore(2)+1:underscore(3)-1));
    spiralNum    = str2double(fileName(underscore(3)+1:dext(1)-1));
    
    fileID = fopen([filePath '/' fileName],'rb');
    kTraj  = fread(fileID,[2 spiralLength*spiralNum],'float');
    fclose(fileID);
    
    kTraj = reshape(kTraj,2,spiralLength,spiralNum);
    kTraj = permute(kTraj,[3 2 1]);
    temp = kTraj;
    temp(:,:,1) = kTraj(:,:,2);
    temp(:,:,2) = kTraj(:,:,1);
    kTraj = temp;
    if showTraj
        figure;
        for n = 1:10:spiralNum
            plot(kTraj(n,:,1),kTraj(n,:,2));hold on;
        end
        hold off; axis equal tight; title(['Spiral Trajectory ' num2str(spiralLength) 'x' num2str(spiralNum) ', showing ' num2str(ceil(spiralNum/10)) ' interleaves']);
    end
catch
    fileID = fopen([filePath '/' fileName],'rb');
    kTraj  = fread(fileID,'float');
    fclose(fileID);
    kTraj = reshape(kTraj,2,[]); 
end