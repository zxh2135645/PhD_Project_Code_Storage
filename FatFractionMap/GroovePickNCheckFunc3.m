function [x_array, y_array, x_centroid_array, y_centroid_array] = GroovePickNCheckFunc3(img_3D, heart_3D, outputFileName, overwrite_label, doublecheck_label)

if nargin == 2
    outputFileName = '';
    overwrite_label = 1;
    doublecheck_label = 0;
elseif nargin == 3
    overwrite_label = 1;
    doublecheck_label = 0;
elseif nargin == 4
    doublecheck_label = 0;
end

% addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
p = 1; % number of landmarks to pick


slc_num = size(img_3D, 3);
if ~(exist(outputFileName, 'file') && overwrite_label == 0)
        x_array = zeros(slc_num, p);
        y_array = zeros(slc_num, p);
        x_centroid_array = zeros(slc_num, 1);
        y_centroid_array = zeros(slc_num, 1);
        for slc = 1:slc_num
            % Load Myocardium mask
            heart_2D = heart_3D(:,:,slc);
            
            [x_heart, y_heart] = find(heart_2D ~= 0);
            x_centroid_array(slc) = round(mean(x_heart),1);
            y_centroid_array(slc) = round(mean(y_heart),1);
            
            figure();
            imagesc(img_3D(:,:,slc))
            truesize([3*size(img_3D,1), 3*size(img_3D,2)]);
            caxis([0 100]);
            axis off;
            
            for i = 1:p
                [y_array(slc,i), x_array(slc,i)] = ginput(1);
                y_array(slc,i) = round(y_array(slc,i), 1);
                x_array(slc,i) = round(x_array(slc,i), 1);
                h1 = text(y_array(slc,i), x_array(slc,i), '*', 'HorizontalAlignment', 'center', 'Color', [0 0 0], 'FontSize', 16);
                h2 = text(y_array(slc,i), x_array(slc,i), num2str(i), 'HorizontalAlignment', 'center', 'Color', [1 0 0], 'FontSize', 14);
            end
        end
        
        close all;
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
        

        
        for slc = 1:slc_num
            figure(), imagesc(img_3D(:,:,slc)), truesize([3*size(img_3D,1), 3*size(img_3D,2)]);
            axis off;
            hold on; plot(coords.y{idx}(slc), coords.x{idx}(slc),'r*', 'MarkerSize', 12)
            hold off;
        end
    else 
        fprintf("The name is invalid!");
    end
end

end