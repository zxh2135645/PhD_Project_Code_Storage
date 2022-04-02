function dicom_idx = ExtractNum(f_glob)

% To include glob function
% addpath('C:/Users/ZhangX1/Documents/MATLAB/cviParser/');

% Reordering
idx_array = zeros(length(f_glob), 1);
for i = 1:length(f_glob)
    B = regexp(f_glob(i),'\d*','Match');
    
    for ii= 1:length(B)
        if ~isempty(B{ii})
            Num(ii,1)=str2double(B{ii}(end));
        else
            Num(ii,1)=NaN;
        end
    end
    % fprintf('%d\n', Num)
    idx_array(i) = Num;
end

% dicom_idx = sort(idx_array);
dicom_idx = idx_array;

end