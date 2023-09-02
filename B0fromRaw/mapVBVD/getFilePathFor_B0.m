%% Get the file path for B0 assuming the file hierarchy is
%  Invivo stuides --->  Code
%                  |->  Data ---> B0
%                             |-> Bz

% Created by Charles on 29 Jan 2019
% Updated by Leon on 20 Aug 2019

function filePath = getFilePathFor_B0()
    % find where is the last slash
    lastSlashPosition = strfind(pwd, '\');
    lastSlashPosition = lastSlashPosition ( length(lastSlashPosition) );

    % Extract the parent path of pwd, "*:\\***\Invivo studies"
    filePath = extractBefore(pwd, lastSlashPosition);
    
    % Add the sub-path
    %filePath = char ( join([filePath, '\Data\B0'], '') ); 
    filePath = char ( join([filePath, '\Data\kspace'], '') );  %update the invivo data sub-path
   
end

