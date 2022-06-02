function [sloc_array] = ExtractSloc(corres_glob)
    sloc_array = zeros(numel(corres_glob), 1);
    for i = 1:numel(corres_glob)
        strings = strsplit(corres_glob{i}, '_');
        last_string = strings{end};
        
        neg_strings = strsplit(last_string, '-');
        if numel(neg_strings) == 2
            last_neg_string = neg_strings{end};
            strings = strsplit(last_neg_string, '.');
            sloc = -(str2double(strings{1}) + str2double(strings{2}) / 100);
        elseif numel(neg_strings) == 1
            pos_string = neg_strings{1};
            pos_strings = strsplit(pos_string, '+');
            last_pos_string = pos_strings{end};
            strings = strsplit(last_pos_string, '.');
            sloc = (str2double(strings{1}) + str2double(strings{2}) / 100);
        end
        
        sloc_array(i) = sloc;
    end
end