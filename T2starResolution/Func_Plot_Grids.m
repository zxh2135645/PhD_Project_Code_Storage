function Func_Plot_Grids(res, half_height)
lim = floor(half_height / res);
for i = 1:lim
    xline(res*(2*i-1)/2, 'Color', [1,1,1]);
    yline(res*(2*i-1)/2, 'Color', [1,1,1]);
    xline(-res*(2*i-1)/2, 'Color', [1,1,1]);
    yline(-res*(2*i-1)/2, 'Color', [1,1,1]);
end
end