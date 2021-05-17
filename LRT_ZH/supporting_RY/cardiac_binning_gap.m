if ~exist('card_roi','var')
    h=figure;
    imshow(dispim(abs(composite_fbp)),[]);
    title('Select Cardiac ROI')
    card_roi=imrect;
    card_roi=padarray(createMask(card_roi),[ovs/2 ovs/2]);
    close(h)
    drawnow
end

missing=round((params.alTR_seconds-params.lSegments*params.lEchoSpacing)/params.lEchoSpacing/2)
Phi_rt=reshape(Phi_rt_init,size(Phi_rt_init,1),Nseg/2,[]);
Phi_rt(:,end+(1:missing),:)=repmat(Phi_rt(:,end,:),[1 missing 1]);
Phi_rt=Phi_rt(:,:);

Nseg2 = params.lSegments+missing*2;
temp=reshape(U_init,[],L_init);
temp=Phi_rt.'*temp(card_roi,:).';
%   temp=navdata_temp;
frame_rate = 1/(params.lEchoSpacing*2);
fs=linspace(-frame_rate/2,frame_rate/2,floor(size(temp,1)/2)*2+1);
if mod(size(temp,1),2)==0 %iseven
    fs=fs(1:(end-1));
end
notchfilt = mod(fs,1/(Nseg2*params.lEchoSpacing));
notchfilt = notchfilt < (fs(2)-fs(1)) | notchfilt > (1/(Nseg2*params.lEchoSpacing)-(fs(2)-fs(1)));
notchfilt = ifftshift(~notchfilt);
temp=ifft(fft(temp).*repmat(notchfilt(:),[1 size(temp,2)]));
[~,~,Segidx,Hidx,~,wall_clock]=self_gate(temp,[55 109]/60,cbins,[10 49]/60,1,frame_rate,Nseg2/2);
Hidx(Segidx>Nseg/2)=[];
wall_clock(Segidx>Nseg/2)=[];
Segidx(Segidx>Nseg/2)=[];