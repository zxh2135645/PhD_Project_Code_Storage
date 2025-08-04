function [Segmentpix, stats, Mask_index] =AHASegmentation_qsm_invivo(Imgin,Maskin,Seg_number,Groove)
% CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
Mask_index = zeros(size(Maskin));
for m=1:size(Imgin,3)
    %find center
    Mask=Maskin(:,:,m);
    [X,Y] = meshgrid(1:size(Mask,2),1:size(Mask,1));
    center=[sum(reshape(Mask.*X,1,[]))/sum(Mask(:)),sum(reshape(Mask.*Y,1,[]))/sum(Mask(:))];
    %dirive angle
    [X,Y] = meshgrid(1:size(Mask,2),1:size(Mask,1));
    X=X-center(1);
    Y=Y-center(2);
    AngleMask=angle(X+1i*Y);
    if m > round(size(Imgin,3)/3*2)
        Segn = 4;
    else
        Segn = Seg_number;
    end
    nintv= 2*pi/Segn;
    %measure amplitude
    Img= (Imgin(:,:,m));
    for n=1:Segn
        nmax=n*nintv+Groove(m)/180*pi;
        nmin=(n-1)*nintv+Groove(m)/180*pi;
        %project to pi and -pi
        nmaxproj=nmax-2*pi*floor((nmax+pi)/2/pi);
        nminproj=nmin-2*pi*floor((nmin+pi)/2/pi);
    %
        if (nminproj>nmaxproj)
            nMask=Mask.*((AngleMask<=nmaxproj)+(AngleMask>=nminproj));
        else
        nMask=Mask.*(AngleMask<=nmaxproj).*(AngleMask>=nminproj);
        end
        Mask_index(:,:,m) = Mask_index(:,:,m) + nMask*n;
        Segmentpix{n,m}=Img(nMask==1);
    
    end
end

% N = round(size(Imgin,3)/3);
N = size(Imgin,3);
for k=1:N
% for k=1:size(Imgin,3)
    
    if (k) > round(N/3*2)
        Segn = 4;
        for n=1:Segn
%             data_temp = [cell2mat(Segmentpix(n,3*k-2));cell2mat(Segmentpix(n,3*k-1)); cell2mat(Segmentpix(n,3*k))];
%             data_temp = [cell2mat(Segmentpix(n,3*k-2));cell2mat(Segmentpix(n,3*k-1))];
%             data_temp = cell2mat(Segmentpix(n,3*k-1));
            data_temp = cell2mat(Segmentpix(n,k));
            if isempty(data_temp)
                stats(:,n,k)=[0,0,0,0,0];
            else
                stats(:,n,k)=[mean(abs(double(data_temp))),mean(double(data_temp)),std(double(data_temp)),median(double(data_temp)),length(double(data_temp))];
            end
            %stats(:,5:6,k) = zeros([2,5])';
        end
    else
        Segn = Seg_number;
        for n=1:Segn
%              data_temp = [cell2mat(Segmentpix(n,3*k-2));cell2mat(Segmentpix(n,3*k-1)); cell2mat(Segmentpix(n,3*k))];
%              data_temp = cell2mat(Segmentpix(n,3*k - 1));
            data_temp = cell2mat(Segmentpix(n,k));
%             data_temp = [cell2mat(Segmentpix(n,1)); cell2mat(Segmentpix(n,2));cell2mat(Segmentpix(n,3));cell2mat(Segmentpix(n,4));cell2mat(Segmentpix(n,5));cell2mat(Segmentpix(n,6))];
            if isempty(data_temp)
                stats(:,n,k)=[0,0,0,0,0];
            else
                stats(:,n,k)=[mean(abs(double(data_temp))),mean(double(data_temp)),std(double(data_temp)),median(double(data_temp)),length(double(data_temp))];
            end
        end
    end
end
end
%%



