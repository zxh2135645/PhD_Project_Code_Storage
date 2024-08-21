clear all;
close all;
%%
% Data is copied from excel spreadsheet 02/26/2023
% TMRP_Repitch

% data = [5	6;	2	3;	2	3;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         4	5;	2	3;	1	2;	1	2;	1	2;	1	2;
%         8	9;	3	4;	3	4;	2	3;	2	3;	2	3;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         4	5;	2	3;	2	3;	1	2;	1	2;	1	2;
%         3	4;	1	2;	1	2;	1	2;	1	2;	1	2;]
% Before the lower bound

hemo_width = [1.2746, 0.5316, 0.9214, 2.3979, 0.5514, 0.5931, 0.5642, 0.5893, 1.1047, 0.6901];
[hemo_width_sorted, I] = sort(hemo_width, 'descend');
data = [4	6;	2	3;	1	3;	1	2;	1	2;	1	2;
        2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
        3	5;	1	3;	1	2;	1	2;	1	2;	1	2;
        6	9;	3	4;	2	4;	2	3;	2	3;	1	3;
        2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
        2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
        2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
        2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
        3	5;	1	3;	1	3;	1	2;	1	2;	1	2;
        2	4;	1	2;	1	2;	1	2;	1	2;	1	2;];

%  pixel size
res_array = [0.3, 0.8, 1.0, 1.3, 1.6, 2.1];

lb = zeros(length(hemo_width), length(res_array));
ub = zeros(length(hemo_width), length(res_array));

for i = 1:length(res_array)
    lb(:,i) = ceil(hemo_width ./ res_array(i));
    ub(:,i) = ceil(sqrt(2)*hemo_width ./ res_array(i)) + 2;
end

y = permute(cat(3, permute(lb, [2 1]), permute(ub, [2 1])), [3 1 2]);

% data = [4	6;	2	3;	1	3;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         3	5;	1	3;	1	2;	1	2;	1	2;	1	2;
%         6	9;	3	4;	2	4;	2	3;	2	3;	1	3;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         2	3;	1	2;	1	2;	1	2;	1	2;	1	2;
%         3	5;	1	3;	1	3;	1	2;	1	2;	1	2;
%         2	4;	1	2;	1	2;	1	2;	1	2;	1	2;];
% 
% y = permute(permute(reshape(data, [6 10 2]),[1 3 2]), [2 1 3]);

y_sorted = y(:,:,I);


x = repmat(1:6, [10 1]) + repmat([0:9]*10, [6 1]).';
x3 = permute(repmat(x, [1 1 2]), [3 2 1]);
%%
figure();
hold on;
for i = 1:size(x3,3)
line(x3(:,1,i), y_sorted(:,1,i),'Color', 'k' , 'LineWidth', 2);
line(x3(:,2:6,i), y_sorted(:,2:6,i), 'Color', [192, 192, 192]/255, 'LineWidth', 2);
plot(x3(:,:,i), y_sorted(2,:,i), 's', 'Color', [231, 126, 144]/255, 'MarkerSize', 6, 'MarkerFaceColor', [231, 126, 144]/255);
plot(x3(:,:,i), y_sorted(1,:,i), 's', 'Color', [85, 188, 194]/255, 'MarkerSize', 6, 'MarkerFaceColor', [85, 188, 194]/255);

end
ylim([0 15]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
