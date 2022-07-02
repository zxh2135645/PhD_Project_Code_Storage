function LocPixCount = AHA_16Seg(ff, myo_ff, BaseGroove, idx_array)
% 06/22/2022 Xinheng(James) Zhang

n = size(ff, 3);
mode = mod(n,3);
integ = fix(n/3);
if n >= 3
    switch mode
        case {0}
            aha_slice = cat(2, repmat([1], [1, integ]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {1}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {2}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ+1]), repmat([3], [1, integ]));
    end
else
    error("Available slice numbers are smaller than 3.");
end

% Basal

%SegTotalPixCount1 = zeros(6, 1);

Groove = BaseGroove + 60;
basal_idx = idx_array(aha_slice == 1);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(ff(:,:,basal_idx), myo_ff(:,:,basal_idx), 6, Groove);

LocPixCount1 = zeros(6, 1);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount1(i) = LocPixCount1(i) + mean(Segmentpix{i,j});

        % SegTotalPixCount1(i) = SegTotalPixCount1(i) + length(Segmentpix{i,j});
    end
end

LocPixCount1 = LocPixCount1 ./ size(Segmentpix, 2);

% Mid-ventricular
LocPixCount2 = zeros(6, 1);
% SegTotalPixCount2 = zeros(6, 1);
Groove = BaseGroove + 60;
mid_idx = idx_array(aha_slice == 2);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(ff(:,:,mid_idx), myo_ff(:,:,mid_idx), 6, Groove);
for i = 1:6
    for j = 1:size(Segmentpix, 2)
        LocPixCount2(i) = LocPixCount2(i) + mean(Segmentpix{i,j});
        % SegTotalPixCount2(i) = SegTotalPixCount2(i) + length(Segmentpix{i,j});
    end
end

LocPixCount2 = LocPixCount2 ./ size(Segmentpix, 2);

% Apical
LocPixCount3 = zeros(4, 1);
SegTotalPixCount3 = zeros(4, 1);
Groove = BaseGroove + 75;
apical_idx = idx_array(aha_slice == 3);

[Segmentpix, stats, Mask_Segn] = AHASegmentation(ff(:,:,apical_idx), myo_ff(:,:,apical_idx), 4, Groove);

for i = 1:4
    for j = 1:size(Segmentpix, 2)
        LocPixCount3(i) = LocPixCount3(i) + mean(Segmentpix{i,j});
        % SegTotalPixCount3(i) = SegTotalPixCount3(i) + length(Segmentpix{i,j});
    end
end

LocPixCount3 = LocPixCount3 ./ size(Segmentpix, 2);
LocPixCount = [LocPixCount1; LocPixCount2; LocPixCount3];