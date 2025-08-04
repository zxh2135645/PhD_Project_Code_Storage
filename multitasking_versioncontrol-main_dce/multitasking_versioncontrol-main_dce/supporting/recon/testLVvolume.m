function diastoleIdx = testLVvolume(img,maskseed)

cbins = size(img,3);
img = abs(img)/max(abs(img(:)));
if nargin < 2 || size(maskseed,1)~=size(img,1) || size(maskseed,2)~=size(img,2)
    maskseed = zeros(size(img(:,:,1)));
    maskseed(floor(size(img,1)/2),floor(size(img,2)/2)) = 1;
end
useMaskseed = true;

b = zeros(size(img));
mask = zeros(size(img));
masksum = zeros(1,cbins);
clusters = 3;
for n = 1:cbins
    temp = imfill(img(:,:,n));
    [a,~] = kmeans(temp(:),clusters);
    a = reshape(a,size(img,1),size(img,2));
    b(:,:,n) = a;
    for nn = 1:clusters
        width(nn) = find(sum(a==nn)>0,1,'last') - find(sum(a==nn)>0,1) + 1;
    end
    [~,idx] = min(width);
    temp = zeros(size(img,1),size(img,2));
    temp(a==idx) = 0.6;
    temp(maskseed==1) = 1;
    b(:,:,n) = temp;
    if useMaskseed
%         clustern = zeros(1,3);
%         for nn = 1:clusters
%             clustern(nn) = sum(maskseed(a==nn));
%         end
%         [~,idx] = max(clustern);
        idxSeed = find(a.*maskseed==idx,1);
        maskseednew = maskseed;
        while isempty(idxSeed)
            se = strel('disk',5);
            maskseednew = imdilate(maskseednew,se);
            idxSeed = find(a.*maskseednew==idx,1);
        end

        roiPosition(1) = mod(idxSeed,size(img,1));
        roiPosition(2) = ceil(idxSeed/size(img,1));
        mask(:,:,n) = grayconnected(a,roiPosition(1),roiPosition(2));
        
        se = strel('disk',15);
        mask(:,:,n) = imclose(mask(:,:,n),se);
        masksum(n)  = sum(sum(mask(:,:,n)>0));
%             sumn(n) = sum(a(:)==idx);
%             mask(:,:,n) = (a==idx);
%             se = strel('disk',7);
%             mask(:,:,n) = imopen(mask(:,:,n),se);
%             mask(:,:,n) = imclose(mask(:,:,n),se);
        %se = strel('disk',7);
        %mask(:,:,n) = imopen(mask(:,:,n),se);
        %masksum(n) = sum(sum(mask(:,:,n)>0));
    else
        sumn(n) = sum(a(:)==idx);
        mask(:,:,n) = (a==idx);
        se = strel('disk',7);
        mask(:,:,n) = imclose(mask(:,:,n),se);
        mask(:,:,n) = imopen(mask(:,:,n),se);
    end
end
if useMaskseed
    tempmean = mean(masksum);
    for n = 2:cbins
        if masksum(n)<tempmean/4
            masksum(n) = masksum(n-1);
        end
    end
    
    masksum = repmat(masksum,1,3);
    masksum = conv(conv(masksum,ones(1,5),'same'),ones(1,3),'same');
    masksum = masksum(cbins+1:cbins*2);
    [~,sumMax] = max(masksum);
    if abs(sumMax - cbins/2 - 1) <= cbins/4
        diastoleIdx = floor(cbins/2) + 1;
    else
        diastoleIdx = 1;
    end
%     implayZoom(b);
%     implayZoom(mask);
else
    diff = abs(mask(:,:,1) - mask(:,:,end));
    diffsum(1) = sum(diff(:));
    for n = 2:cbins
        diff = abs(mask(:,:,n) - mask(:,:,n-1));
        diffsum(n) = sum(sum(diff(:)));
    end
    diffsum = repmat(diffsum,1,3);
    diffsum = conv(conv(diffsum,ones(1,5),'same'),ones(1,5),'same');
    diffsum = diffsum(cbins+1:cbins*2);
    %figure;plot(diffsum);
    %implayZoom(img/max(abs(img(:))));
    [~,diffMin] = min(diffsum);
    if abs(diffMin - cbins/2 - 1) <= cbins/4
        diastoleIdx = floor(cbins/2) + 1;
    else
        diastoleIdx = 1;
    end
end

% temp1 = img(:,:,1);
% [a,~] = kmeans(temp1(:),3);
% a = reshape(a,size(img,1),size(img,2));
% for n = 1:3
%     width(n) = find(sum(a==n)>0,1,'last') - find(sum(a==n)>0,1) + 1;
% end
% [~,i1] = min(width);
% sum1 = sum(a(:)==i1);
% 
% temp2 = img(:,:,floor(size(img,3)/2)+1);
% [a,~] = kmeans(temp2(:),3);
% a = reshape(a,size(img,1),size(img,2));
% for n = 1:3
%     width(n) = find(sum(a==n)>0,1,'last') - find(sum(a==n)>0,1) + 1;
% end
% [~,i2] = min(width);
% sum2 = sum(a(:)==i2);
% 
% if sum2 > sum1
%     diastoleIdx = floor(size(img,3)/2)+1;
% else
%     diastoleIdx = 1;
% end
