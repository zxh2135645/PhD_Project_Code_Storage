function center_coords = DivideMyoInHalf(epi_coords, endo_coords)
num_pt = length(endo_coords);
x1 = endo_coords(:,1);
y1 = endo_coords(:,2);
x2 = epi_coords(:,1);
y2 = epi_coords(:,2);

center_coords = zeros(num_pt, 2);
for i = 1:num_pt
    distances = sqrt((x1(i) - x2) .^ 2 + (y1(i) - y2) .^ 2);
    [minDistance, indexOfMin] = min(distances);
    center_coords(i, 1) = mean([x1(i), x2(indexOfMin)]);
    center_coords(i, 2) = mean([y1(i), y2(indexOfMin)]);
end

end