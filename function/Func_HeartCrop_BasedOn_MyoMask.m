function [cropped_img] = Func_HeartCrop_BasedOn_MyoMask(myomask, img)
% myomask and img are 2D
    sz1 = mean(myomask);
    idx1 = sz1 > 0;
    sz2 = mean(myomask);
    idx2 = sz2 > 0;
    
    
end