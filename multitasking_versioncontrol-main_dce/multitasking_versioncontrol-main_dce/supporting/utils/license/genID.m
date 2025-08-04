function [myID, strAdd] = genID()

if ispc
    [~,info] = system('getmac');
elseif ismac
    [~,info] = system('ifconfig en0 | grep ether');
elseif isunix
    [~,info] = system('ip addr | grep ether');
else
    error('OS not recognized.');
end

myID = double(info) * sum(double('multitasking'));

expression = '\w\w(-|:)\w\w(-|:)\w\w(-|:)\w\w(-|:)\w\w(-|:)\w\w';
strAdd = regexp(info,expression,'match');
for n = length(strAdd):-1:1
    if strcmp(strAdd{n},'ff:ff:ff:ff:ff:ff') || strcmp(strAdd{n},'ff-ff-ff-ff-ff-ff')
        strAdd(n) = [];
    end
end
