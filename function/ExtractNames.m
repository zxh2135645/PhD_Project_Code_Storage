function [Names] = ExtractNames(name_glob)
    Names = cell(length(name_glob), 1);
    
    for i = 1:length(Names)
        strings = strsplit(name_glob{i}, '\');
        Names{i} = strings{end-1};
    end
end