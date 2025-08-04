function res = im2row(im, winSize)
%res = im2row(im, winSize)

if numel(winSize) == 2 %2D
    [sx,sy,sz] = size(im);

    res = single(zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz));
    count=0;
    for y=1:winSize(2)
        for x=1:winSize(1)
            count = count+1;
            res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
                (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
        end
    end
elseif numel(winSize) == 3 %3D
    [sy,sz,sx,sc] = size(im);
    
    res = single(zeros((sy-winSize(1)+1)*(sz-winSize(2)+1)*(sx-winSize(3)+1),prod(winSize),sc));
    count = 0;
    for x = 1:winSize(3)
        for z = 1:winSize(2)
            for y = 1:winSize(1)
                count = count + 1;
                res(:,count,:) = reshape(im(y:sy-winSize(1)+y,z:sz-winSize(2)+z,x:sx-winSize(3)+x,:),...
                                 (sy-winSize(1)+1)*(sz-winSize(2)+1)*(sx-winSize(3)+1),1,sc);
            end
        end
    end
end
