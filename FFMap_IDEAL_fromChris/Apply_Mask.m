%% Creat Mask
% Creat mask

if ~exist([MRdat_path,'Manual_Mask.mat'])% if original mask doesn't exist
    mag_all = sum(squeeze(abs(AllPhasemap(n).compleximg(:,:,:,:,end))),4).*AllPhasemap(1).Mask;
      RD_FOV=[(round(size(mag_all,1)/2)+1-round(size(mag_all,1)/4))...
        :(round(size(mag_all,1)/2)+1+round(size(mag_all,1)/4))];
    mag_all=permute(mag_all, [3 1 2]);% get the side view and readout to be dim 2

    centMask=zeros(size(mag_all));
    centMask(1:end, RD_FOV,1:end,1:end)=1;
    mag_all=mag_all.*centMask;
   % AllPhasemap(n).Mask=ones(size(mag_all));
    AllPhasemap(n).Mask=permute(AllPhasemap(n).Mask,[3 1 2]);
    Generatecontour=1;
    for ms=1:2:(size(mag_all,3)-1);
        mag=mag_all(:,:,ms);
        if FlagMask
            if (IsManual_Mask && AllPhasemap(n).ManualMask==0)
                if n==1;
                    %Charles 27Jun2019, add to show slice# as well
                    disp(['Mask ', num2str(ms/size(AllPhasemap(n).Mask,3)*100),' %,', num2str(ms), '/', num2str(size(AllPhasemap(n).Mask,3)), 'slice'])% to display in a better orientation
                    %                     figure(1);imagesc(mag'); axis equal ;
                    %                     colormap gray;
                    %                     Mask_m=roipoly';
                    %                     AllPhasemap(n).Mask(:,:,ms)=AllPhasemap(n).Mask(:,:,ms).*Mask_m;
                    %                     AllPhasemap(n).Mask(:,:,ms+1)=AllPhasemap(n).Mask(:,:,ms+1).*Mask_m;
                    nmag=mag/max(mag(:));%for imshow normalize the figure
                    if ms==1
                        imshow(nmag);
                        set(gcf, 'units', 'normalized', 'Position', [0.045, 0.3, 0.4, 0.5]) %Charles 30Jan2019, show the figure at a better location
                        %Orginally: set(gcf, 'Position', [100, 100, 1000, 1000])
                        [Mask_m tmp_xi tmp_yi]=roipoly;
                        
                    else
                        imshow(nmag); 
                        set(gcf, 'units', 'normalized', 'Position', [0.045, 0.3, 0.4, 0.5]) %Charles 30Jan2019, show the figure at a better location
                        %Orginally: set(gcf, 'Position', [100, 100, 1000, 1000])
                        h = impoly(gca,[tmp_xi tmp_yi]); 
                        Generatecontour=input('Re-contour?');
                        Mask_m = createMask(h); %BW contains the mask which you just altered.
                        pos = getPosition(h); 
                        tmp_xi=pos(:,1);
                        tmp_yi=pos(:,2);
                        if Generatecontour % overwrite mask if need more points
                            [Mask_m tmp_xi tmp_yi]=roipoly(nmag);
                        end
                        
                    end
                    AllPhasemap(n).Mask(:,:,ms)=AllPhasemap(n).Mask(:,:,ms).*Mask_m;
                    AllPhasemap(n).Mask(:,:,ms+1)=AllPhasemap(n).Mask(:,:,ms+1).*Mask_m;
                    if ms==(size(AllPhasemap(n).Mask,3)-1);
                        AllPhasemap(n).ManualMask=1;
                    end
                else
                    AllPhasemap(n).Mask(:,:,ms)= AllPhasemap(1).Mask(:,:,ms);
                    if ms==(size(AllPhasemap(n).Mask,3)-1);
                        AllPhasemap(n).ManualMask=1;
                    end
                end
                
            end
            % mask = repmat(AllPhasemap(n).Mask(:,:,ms),[1,1,size(phase,3)]);
            % phase(~mask)=nan;
            
        end
    end
    
    
    %fill holes
    for s=1:size(AllPhasemap(n).Mask,3)
        AllPhasemap(n).Mask(:,:,s)= imfill(AllPhasemap(n).Mask(:,:,s),'holes');
    end
    
    
    AllPhasemap(n).Mask=permute(AllPhasemap(n).Mask,[2 3 1]);
    Mask=AllPhasemap(n).Mask;
    save([MRdat_path,'Manual_Mask.mat'],'Mask');
else
    load([MRdat_path,'Manual_Mask.mat']);
    AllPhasemap(n).Mask=Mask;
end
    %clear Mask