function coords = GroovePickNCheckFunc2(Names, label, anatomy, base_dir, out_dir, overwrite_label, doublecheck_label)

if nargin == 5
    overwrite_label = 1;
    doublecheck_label = 0;
elseif nargin == 6
    doublecheck_label = 0;
end

% addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
p = 1; % number of landmarks to pick
x = cell(length(Names), 1); % Rows
y = cell(length(Names), 1); % Columns
x_centroid = cell(length(Names), 1);
y_centroid = cell(length(Names), 1);
coords = struct;
outputFileName = [out_dir, label, 'coords2.mat'];

if ~(exist(outputFileName, 'file') && overwrite_label == 0)
    for n = 1:length(Names)
        name = Names{n};
        
        %%
        mat_glob = glob([base_dir, name, '\', label, '\*\*']);
        mat_file = GetFullPath(cat(2, mat_glob{1}, 'VOLUME_IMAGE.mat'));
        load(mat_file, 'volume_image');
        myocardium_glob = glob(cat(2, out_dir, name, '/', label, '/', anatomy, '/masked_*.mat'));

        % Reordering
        idx_array = zeros(length(myocardium_glob), 1);
        for i = 1:length(myocardium_glob)
            B = regexp(myocardium_glob(i),'\d*','Match');
            
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
        %%
        x_array = zeros(length(idx_array), p);
        y_array = zeros(length(idx_array), p);
        x_centroid_array = zeros(length(idx_array), 1);
        y_centroid_array = zeros(length(idx_array), 1);
        for slice_num = 1:length(idx_array)
            % Load Myocardium mask
            mask_path = cat(2, out_dir, name, '\', label, '\', anatomy, '\', 'masked_heart', num2str(idx_array(slice_num)), '.mat');
            load(mask_path, 'mask_heart');
            
            [x_heart, y_heart] = find(mask_heart ~= 0);
            x_centroid_array(slice_num) = round(mean(x_heart),1);
            y_centroid_array(slice_num) = round(mean(y_heart),1);
            
            figure();
            imagesc(volume_image(:,:,idx_array(slice_num)))
            truesize([3*size(volume_image,1), 3*size(volume_image,2)]);
            axis off;
            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % =============================== code for picking points =======
         % ===============================================================
            for i = 1:p
                [y_array(slice_num,i), x_array(slice_num,i)] = ginput(1);
                y_array(slice_num,i) = round(y_array(slice_num,i), 1);
                x_array(slice_num,i) = round(x_array(slice_num,i), 1);
                h1 = text(y_array(slice_num,i), x_array(slice_num,i), '*', 'HorizontalAlignment', 'center', 'Color', [0 0 0], 'FontSize', 16);
                h2 = text(y_array(slice_num,i), x_array(slice_num,i), num2str(i), 'HorizontalAlignment', 'center', 'Color', [1 0 0], 'FontSize', 14);
            end
        end
        close all;
        x{n} = x_array;
        y{n} = y_array;
        x_centroid{n} = x_centroid_array;
        y_centroid{n} = y_centroid_array;
    end
    coords.x = x;
    coords.y = y;
    coords.x_centroid = x_centroid;
    coords.y_centroid = y_centroid;
    save(outputFileName, 'coords');
else
    fprintf('The file already exists for %s\n', label);
end

if doublecheck_label == 1
    % double check saved coordinates
    prompt = 'Please type in a patient name: ';
    name = input(prompt, 's');
    load(outputFileName, 'coords');
    idx = find(strcmp(name, Names), 1);
    if ~isempty(idx)
        
        mat_glob = glob([base_dir, name, '\', label, '\*\*']);
        mat_file = GetFullPath(cat(2, mat_glob{1}, 'VOLUME_IMAGE.mat'));
        load(mat_file, 'volume_image');
        myocardium_glob = glob(cat(2, out_dir, name, '/', label, '/', anatomy, '/masked_*.mat'));
        
        % Reordering
        idx_array = zeros(length(myocardium_glob), 1);
        for j = 1:length(myocardium_glob)
            B = regexp(myocardium_glob(j),'\d*','Match');
            
            for ii= 1:length(B)
                if ~isempty(B{ii})
                    Num(ii,1)=str2double(B{ii}(end));
                else
                    Num(ii,1)=NaN;
                end
            end
            % fprintf('%d\n', Num)
            idx_array(j) = Num;
        end
        
        for slice_num = 1:length(idx_array)
            figure(), imagesc(volume_image(:,:,idx_array(slice_num))), truesize([3*size(volume_image,1), 3*size(volume_image,2)]);
            axis off;
            hold on; plot(coords.y{idx}(slice_num), coords.x{idx}(slice_num),'r*', 'MarkerSize', 12)
            hold off;
        end
    else 
        fprintf("The name is invalid!");
    end
end

end