if ~strcmp(ScanType,'Cine')
    if ~exist('card_roi','var')
        h=figure;
        imshow(dispim(abs(composite_fbp)),[]);
        title('Select Cardiac ROI')
        card_roi=imrect;
        card_roi=padarray(createMask(card_roi),[ovs/2 ovs/2]);
        close(h)
        drawnow
    end
    Nseg = params.lSegments;
    temp=reshape(U_init,[],L_init);
    temp=Phi_rt_init.'*temp(card_roi,:).';
    %   temp=navdata_temp;
    frame_rate = 1/(params.lEchoSpacing*2);
    fs=linspace(-frame_rate/2,frame_rate/2,floor(size(temp,1)/2)*2+1);
    if mod(size(temp,1),2)==0 %iseven
        fs=fs(1:(end-1));
    end
    notchfilt = mod(fs,1/(Nseg*params.lEchoSpacing));
    notchfilt = notchfilt < (fs(2)-fs(1)) | notchfilt > (1/(Nseg*params.lEchoSpacing)-(fs(2)-fs(1)));
    notchfilt = ifftshift(~notchfilt);
    temp=ifft(fft(temp).*repmat(notchfilt(:),[1 size(temp,2)]));
    [~,~,Segidx,Hidx,~,wall_clock]=self_gate(temp,[55 109]/60,cbins,[10 49]/60,1,frame_rate,Nseg/2);
else
    Nseg = 2;
    navdata_temp=nav_data(:,:);
%     navdata_temp=reshape(nav_data(:,:,:,1),size(nav_data,1),[]); %assumes 1st "coil" is virtual coil covering cardiac region
    [~,~,Segidx,Hidx,~,wall_clock]=self_gate(navdata_temp,[55 109]/60,cbins,[10 49]/60,1,1/(params.lEchoSpacing*2),Nseg/2);
end