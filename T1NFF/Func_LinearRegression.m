function [x1, b1, X1, yCalc1, yCalc1_for_plot, Rsq1, str1] = ...
    Func_LinearRegression(Mipix_mean_x_nnz, Mipix_mean_y_nnz, x1_lb, x1_ub)
X1 = [ones(length(Mipix_mean_x_nnz), 1), Mipix_mean_x_nnz'];
b1 = X1\Mipix_mean_y_nnz';


if ~isempty(b1)
    x1 = x1_lb:x1_ub;
    X1_for_plot = [ones(length(x1), 1), x1'];
    yCalc1 = X1*b1;
    yCalc1_for_plot = X1_for_plot*b1;
    Rsq1 = 1 - sum((Mipix_mean_y_nnz' - yCalc1).^2)/sum((Mipix_mean_y_nnz' - mean(Mipix_mean_y_nnz)).^2);
    str1 = cat(2, 'Y = ', num2str(round(b1(2),1)), 'X + ', num2str(round(b1(1),1)), ', R^2 = ', num2str(round(Rsq1,2)));
else
    x1 = [];
    yCalc1 = [];
    yCalc1_for_plot = [];
    Rsq1 = 'NA';
    str1 = 'NA';
end

end