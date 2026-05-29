function [Names] = ExtractNames(name_glob)
    Names = cell(length(name_glob), 1);
    
    for i = 1:length(Names)
        if ispc % Check if operation system is windows
            strings = strsplit(name_glob{i}, '\');
        else
            strings = strsplit(name_glob{i}, '/');
        end
        Names{i} = strings{end-1};
    end
end