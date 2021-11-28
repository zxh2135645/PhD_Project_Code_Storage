function [list_to_read, order_to_read] = NamePicker(f_glob, label, sel_array)
if nargin == 1
    label = 0;
    dst_name = ExtractNames(f_glob);
elseif label == 1
    dst_name = f_glob;
elseif label == 0
    dst_name = ExtractNames(f_glob);
end

disp(dst_name);

if ~exist('sel_array', 'var')
    sel_array = input('Please add an array here:  ');
end

char_array = num2str(sel_array', '%04.f');


ind_array2 = zeros(size(dst_name, 1), 1);
for i = 1:size(char_array, 1)
    cha = char_array(i, :);
    ind_array2 = ind_array2 + contains(dst_name, cha);
end
% Checking without 0s ahead
if sum(ind_array2) == 0
    char_array = num2str(sel_array');
    for i = 1:size(char_array, 1)
        cha = strtrim(char_array(i, :));
        ind_array2 = ind_array2 + contains(dst_name, cha);
    end
end

ind_array3 = find(ind_array2 == 1);
list_to_read = f_glob(ind_array3);

if label == 0
    name_to_compare = ExtractNames(list_to_read);
elseif label == 1
    name_to_compare = list_to_read;
end

order_to_read = zeros(length(list_to_read), 1);
for i = 1:length(list_to_read)
    order_to_read(i) = find(contains(name_to_compare, strtrim(char_array(i, :))) == 1);
end

end