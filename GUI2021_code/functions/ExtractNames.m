function [Names] = ExtractNames(name_glob)
Names = cell(1, length(name_glob));
if ispc
    for i = 1:length(Names)
        strings = strsplit(name_glob{i}, '\');
        Names{i} = strings{end-1};
    end
else
    for i = 1:length(Names)
        strings = strsplit(name_glob{i}, '/');
        Names{i} = strings{end-1};
    end
end
end