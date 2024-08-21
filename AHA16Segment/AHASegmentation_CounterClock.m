function [Segmentpix, stats, Mask_Segn] =AHASegmentation_CounterClock(Imgin,Maskin,Segn,Groove,Endoin)
% Created in 04/19/2023 for counterclockwise AHA Segmentation
% Modified by Xinheng Zhang on 07/30/2020
%  Mistakes in assigning Mask_Segn, resolved. 06/22/2022
% Added an alternative way to find center, esp useful for imbalanced
% myocardium shape
Mask_Segn = zeros(size(Imgin));

for m=1:size(Imgin,3)
    %find center
    Mask=Maskin(:,:,m);
    if nargin == 5 % 5 input arguments
        Endo = Endoin(:,:,m);
        C = regionprops(Endo);
        center = C.Centroid;
    else
        [X,Y] = meshgrid(1:size(Mask,2),1:size(Mask,1));
        center=[sum(reshape(Mask.*X,1,[]))/sum(Mask(:)),sum(reshape(Mask.*Y,1,[]))/sum(Mask(:))];
    end
    
    %derive angle
    [X,Y] = meshgrid(1:size(Mask,2),1:size(Mask,1));
    X=X-center(1);
    Y=Y-center(2);
    % AngleMask=angle(X+1i*Y);
    AngleMask=angle(X-1i*Y); % Change to counterclockwise AngleMask
    nintv= 2*pi/Segn;
    Img=Imgin(:,:,m);
    
    
    for n=1:Segn
        nmax=n*nintv+Groove/180*pi;
        nmin=(n-1)*nintv+Groove/180*pi;
        %project to pi and -pi
        nmaxproj=nmax-2*pi*floor((nmax+pi)/2/pi);
        nminproj=nmin-2*pi*floor((nmin+pi)/2/pi);
        %
        if (nminproj>nmaxproj)
            nMask=Mask.*((AngleMask<=nmaxproj)+(AngleMask>=nminproj));
        else
            nMask=Mask.*(AngleMask<=nmaxproj).*(AngleMask>=nminproj);
        end
        Segmentpix{n,m}=Img(nMask==1);
        stats(:,n,m)=[mean(Img(nMask==1)),std(Img(nMask==1)),median(Img(nMask==1)),length(Img(nMask==1))];
        
        Mask_Segn(:,:,m) = Mask_Segn(:,:,m) + nMask * n;
        
    end
end
end



